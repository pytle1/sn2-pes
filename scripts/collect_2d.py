from pathlib import Path
import csv
import re


PROJECT = Path(__file__).resolve().parents[1]
CALCULATIONS = PROJECT / "calculations"

HARTREE_TO_KJMOL = 2625.49962

CN_VALUES = [round(1.9 + 0.1 * i, 1) for i in range(7)]
CCL_VALUES = [round(1.9 + 0.1 * i, 1) for i in range(7)]

POINT_RE = re.compile(r"relaxed_point_cn([0-9.]+)_ccl([0-9.]+)")


def get_energy(output_file: Path) -> float:
    text = output_file.read_text()

    if "THE OPTIMIZATION HAS CONVERGED" not in text:
        text = output_file.read_text()

    matches = re.findall(
            r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
            text,
    )

    if not matches:
        raise ValueError(f"No energy found in {output_file}")

    return float(matches[-1])


def collect_points():
    points = []

    allowed_cn = set(CN_VALUES)
    allowed_ccl = set(CCL_VALUES)

    for directory in CALCULATIONS.glob("relaxed_point_cn*_ccl*"):
        match = POINT_RE.fullmatch(directory.name)

        if match is None:
            continue

        cn = round(float(match.group(1)), 1)
        ccl = round(float(match.group(2)), 1)

        if cn not in allowed_cn or ccl not in allowed_ccl:
            continue

        output_file = directory / "job.out"
        xyz_file = directory / "job.xyz"

        if not output_file.exists() or not xyz_file.exists():
            continue

        energy = get_energy(output_file)

        points.append(
                {
                    "cn": cn,
                    "ccl": ccl,
                    "energy": energy,
                    "xyz": str(xyz_file),
                }
        )

    return points


def main():
    points = collect_points()

    if len(points) != 49:
        raise ValueError(
                f"Expected 49 grid points, found {len(points)}"
        )

    points.sort(key=lambda point: (point["cn"], point["ccl"]))

    min_energy = min(point["energy"] for point in points)

    for point in points:
        point["relative_energy"] = (
                point["energy"] - min_energy
        ) * HARTREE_TO_KJMOL

    output_csv = CALCULATIONS / "scan_2d.csv"

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
                f,
                fieldnames=[
                    "point",
                    "cn",
                    "ccl",
                    "energy",
                    "relative_energy",
                    "xyz",
                ],
        )
        writer.writeheader()
        writer.writerows(points)

    print(f"Saved {output_csv}")
    print()

    for point in points:
        print(
                f"CN={point['cn']:.1f}  "
                f"CCl={point['ccl']:.1f}  "
                f"dE={point['relative_energy']:8.2f} kJ/mol"
        )


if __name__ == "__main__":
    main()
