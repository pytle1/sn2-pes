from pathlib import Path
import math
import re
import shutil
import subprocess
import argparse

PROJECT = Path(__file__).resolve().parents[1]

source_xyz = (
        PROJECT
        / "calculations"
        / "relaxed_point_cn1.600_ccl2.800" 
        / "job.xyz"
)

# ORCA uses 0-based atom indices.
IPSO_C = 9
NITRO_N = 2
CL = 24

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


def report_distances(
        xyz_path: Path,
):
    atoms = read_xyz(xyz_path)
    cn = distance(atoms[IPSO_C][1], atoms[NITRO_N][1])
    ccl = distance(atoms[IPSO_C][1], atoms[CL][1])
    print(f"C-N  = {cn:.6f} Å  ")
    print(f"C-Cl = {ccl:.6f} Å  ")


def optimize_point(
        source_xyz: Path,
):
    workdir = (
            PROJECT
            / "calculations"
            / "radical"
            / f"D_unconstrained_opt"
    )

    orca = shutil.which("orca")
    if orca is None:
        raise RuntimeError("Could not find ORCA executable in PATH")

    workdir.mkdir(parents=True, exist_ok=True)
    start_xyz = workdir / "start.xyz"
    shutil.copy2(source_xyz, start_xyz)

    inp = workdir / "job.inp"
    out = workdir / "job.out"

    inp.write_text(
        f"""! r2SCAN-3c Opt TightSCF

%pal
  nprocs {NCORES}
end

* xyzfile 0 2 start.xyz
"""
    )
    
    print("Radical pathway: charge 0, multiplicity 2")
    print("Starting geometry:")
    report_distances(start_xyz)
    print(f"\nRunning ORCA in {workdir}")

    with out.open("w") as fh:
        result = subprocess.run(
            [orca, inp.name],
            cwd=workdir,
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

    optimized_xyz = workdir / "job.xyz"
    if not optimized_xyz.exists():
        raise RuntimeError(
            f"ORCA terminated normally, but {optimized_xyz} was not found"
        )

    print(f"\nFinal energy: {float(energies[-1]):.12f} Eh")
    print("Optimized geometry:")
    report_distances(optimized_xyz)
    print(f"\nXYZ: {optimized_xyz}")
    print(f"Output: {out}")

    return {
            "energy": float(energies[-1]),
            "xyz": str(optimized_xyz),
    }


def parse_args():
    parser = argparse.ArgumentParser(
            description="Run one constrained C-N / C-Cl optimization."
    )

    parser.add_argument(
            "--xyz",
            type=Path,
            required=True,
            help="Starting XYZ geometry",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    optimize_point(
            source_xyz=args.xyz,
    )
