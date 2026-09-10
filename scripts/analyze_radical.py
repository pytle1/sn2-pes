from pathlib import Path
import argparse
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description="Frequency and spin-density analysis of a radical geometry."
    )

    parser.add_argument(
        "--xyz",
        type=Path,
        required=True,
        help="Optimized radical XYZ geometry",
    )

    parser.add_argument(
        "--name",
        required=True,
        help="Name of the species, e.g. D_candidate",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    project = Path(__file__).resolve().parents[1]

    workdir = (
        project
        / "calculations"
        / "radical"
        / args.name
        / "freq_spin"
    )

    workdir.mkdir(parents=True, exist_ok=True)

    start_xyz = workdir / "geometry.xyz"
    shutil.copy2(args.xyz, start_xyz)

    inp = workdir / "analysis.inp"
    out = workdir / "analysis.out"

    inp.write_text(
        """! r2SCAN-3c Freq TightSCF

%pal
  nprocs 4
end

%output
  Print[P_Hirshfeld] 1
  Print[P_Mayer] 1
end

%plots
  Format Gaussian_Cube
  SpinDens("spin_density.cube");
end

* xyzfile 0 2 geometry.xyz
"""
    )

    orca = shutil.which("orca")

    if orca is None:
        raise RuntimeError("ORCA not found in PATH")

    orca = Path(orca).resolve()

    print(f"Geometry: {start_xyz}")
    print("Electronic state: charge 0, multiplicity 2")
    print(f"Running analysis in {workdir}")

    with out.open("w") as f:
        result = subprocess.run(
            [str(orca), inp.name],
            cwd=workdir,
            stdout=f,
            stderr=subprocess.STDOUT,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"ORCA failed. See {out}"
        )

    text = out.read_text()

    if "ORCA TERMINATED NORMALLY" not in text:
        raise RuntimeError(
            f"ORCA did not terminate normally. See {out}"
        )

    print("Analysis completed.")
    print(f"Output:       {out}")
    print(f"Spin density: {workdir / 'spin_density.cube'}")


if __name__ == "__main__":
    main()
