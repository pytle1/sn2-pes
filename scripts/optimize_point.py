from pathlib import Path
import math
import re
import shutil
import subprocess

PROJECT = Path(__file__).resolve().parents[1]

SOURCE_XYZ = (
        PROJECT
        / "calculations"
        / "relaxed_point_cn1.600_ccl2.800" 
        / "job.xyz"
)

WORKDIR = PROJECT / "calculations" / "relaxed_point_cn1.800_ccl2.600"

# ORCA uses 0-based atom indices.
IPSO_C = 9
NITRO_N = 2
CL = 24

TARGET_CN = 1.80
TARGET_CCL = 2.60
NCORES = 4


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    natoms = int(lines[0])
    rows = []
    for line in lines[2 : 2 + natoms]:
        fields = line.split()
        rows.append((fields[0], tuple(map(float, fields[1:4]))))
    return rows


def distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def report_distances(xyz_path: Path):
    atoms = read_xyz(xyz_path)
    cn = distance(atoms[IPSO_C][1], atoms[NITRO_N][1])
    ccl = distance(atoms[IPSO_C][1], atoms[CL][1])
    print(f"C-N  = {cn:.6f} Å  (target {TARGET_CN:.3f})")
    print(f"C-Cl = {ccl:.6f} Å  (target {TARGET_CCL:.3f})")


def main():
    orca = shutil.which("orca")
    if orca is None:
        raise RuntimeError("Could not find ORCA executable in PATH")

    WORKDIR.mkdir(parents=True, exist_ok=True)
    start_xyz = WORKDIR / "start.xyz"
    shutil.copy2(SOURCE_XYZ, start_xyz)

    inp = WORKDIR / "job.inp"
    out = WORKDIR / "job.out"

    inp.write_text(
        f"""! r2SCAN-3c Opt TightSCF

%pal
  nprocs {NCORES}
end

%geom
  Constraints
    {{ B {IPSO_C} {NITRO_N} {TARGET_CN:.6f} C }}
    {{ B {IPSO_C} {CL} {TARGET_CCL:.6f} C }}
  end
end

* xyzfile -1 1 start.xyz
"""
    )

    print("Starting geometry:")
    report_distances(start_xyz)
    print(f"\nRunning ORCA in {WORKDIR}")

    with out.open("w") as fh:
        result = subprocess.run(
            [orca, inp.name],
            cwd=WORKDIR,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )

    text = out.read_text(errors="replace")
    if result.returncode != 0 or "ORCA TERMINATED NORMALLY" not in text:
        raise RuntimeError(f"ORCA failed; inspect {out}")

    if "THE OPTIMIZATION HAS CONVERGED" not in text:
        raise RuntimeError(f"Geometry optimization did not converge; inspect {out}")

    energies = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text
    )
    if not energies:
        raise RuntimeError(f"Could not find final energy in {out}")

    optimized_xyz = WORKDIR / "job.xyz"
    if not optimized_xyz.exists():
        raise RuntimeError(
            f"ORCA terminated normally, but {optimized_xyz} was not found"
        )

    print(f"\nFinal energy: {float(energies[-1]):.12f} Eh")
    print("Optimized geometry:")
    report_distances(optimized_xyz)
    print(f"\nXYZ: {optimized_xyz}")
    print(f"Output: {out}")


if __name__ == "__main__":
    main()

