from pathlib import Path
import csv
import re

PROJECT = Path(__file__).resolve().parents[1]
CALCULATIONS = PROJECT / "calculations"

HARTREE_TO_KJMOL = 2625.49962


def get_energy(output_file: Path) -> float:
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

    pattern = re.compile(
            r"relaxed_point_cn([0-9.]+)_ccl([0-9.]+)"
    )

    for directory in CALCULATIONS.glob("relaxed_point_cn*_ccl*"):
        match = pattern.fullmatch(directory.name)

        if match is None:
            continue

        cn = float(match.group(1))
        ccl = float(match.group(2))

        output_file = directory / "job.out"

        if not output_file.exists():
            continue

        energy = get_energy(output_file)

        points.append(
                {
                    "cn": cn,
                    "ccl": ccl,
                    "energy": energy,
                    "xyz": str(directory / "job.xyz"),
                }
        )

    return points


def main():
    points = collect_points()

    points.sort(key=lambda point: point["cn"])

    reference_energy = points[0]["energy"]

    for i, point in enumerate(points):
        point["point"] = i
        point["relative_energy"] = (
                point["energy"] - reference_energy
        ) * HARTREE_TO_KJMOL

    output_csv = CALCULATIONS / "relaxed_scan.csv"

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

    for point in points:
        print(
                f"{point['point']:2d}  "
                f"CN={point['cn']:.2f}  "
                f"CCl={point['ccl']:.2f}  "
                f"dE={point['relative_energy']:8.2f} kJ/mol"
        )


if __name__ == "__main__":
    main()
