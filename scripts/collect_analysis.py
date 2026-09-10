from pathlib import Path
import csv
import re
import math


PROJECT = Path(__file__).resolve().parents[1]
CALCULATIONS = PROJECT / "calculations"

ENERGY_CSV = CALCULATIONS / "scan_2d.csv"
OUTPUT_CSV = CALCULATIONS / "scan_2d_analysis.csv"


def find_bond_order(text, atom1, elem1, atom2, elem2):
    patterns = [
            rf"B\(\s*{atom1}-{elem1}\s*,\s*{atom2}-{elem2}\s*\)\s*:\s*([+-]?\d+\.\d+)",
            rf"B\(\s*{atom2}-{elem2}\s*,\s*{atom1}-{elem1}\s*\)\s*:\s*([+-]?\d+\.\d+)",
    ]

    for pattern in patterns:
        match = re.search(pattern, text)

        if match:
            return float(match.group(1))

    return math.nan


def hirshfeld_charge(text, atom, element):
    if "HIRSHFELD ANALYSIS" not in text:
        return math.nan

    hirshfeld = text.split("HIRSHFELD ANALYSIS", 1)[1]

    pattern = rf"^\s*{atom}\s+{element}\s+([+-]?\d+\.\d+)"

    match = re.search(
            pattern,
            hirshfeld,
            flags=re.MULTILINE,
    )

    if match is None:
        return math.nan

    return float(match.group(1))


def main():
    rows = []

    with ENERGY_CSV.open() as f:
        reader = csv.DictReader(f)

        for row in reader:
            cn = float(row["cn"])
            ccl = float(row["ccl"])

            directory = (
                    CALCULATIONS
                    / f"relaxed_point_cn{cn:.3f}_ccl{ccl:.3f}"
            )

            analysis_file = directory / "analysis.out"

            if not analysis_file.exists():
                raise FileNotFoundError(analysis_file)

            text = analysis_file.read_text()

            bo_cn = find_bond_order(
                    text,
                    2, "N",
                    9, "C",
            )

            bo_ccl = find_bond_order(
                    text,
                    9, "C",
                    24, "Cl",
            )

            bo_difference = bo_ccl - bo_cn

            q_o0 = hirshfeld_charge(text, 0, "O")
            q_o1 = hirshfeld_charge(text, 1, "O")
            q_n = hirshfeld_charge(text, 2, "N")
            q_cipso = hirshfeld_charge(text, 9, "C")
            q_cl = hirshfeld_charge(text, 24, "Cl")

            q_no2 = q_o0 + q_o1 + q_n
            q_difference = q_no2 - q_cl

            rows.append(
                    {
                        "cn": cn,
                        "ccl": ccl,
                        "energy": float(row["energy"]),
                        "relative_energy": float(row["relative_energy"]),
                        "bo_cn": bo_cn,
                        "bo_ccl": bo_ccl,
                        "bo_difference": bo_difference,
                        "q_cl": q_cl,
                        "q_cipso": q_cipso,
                        "q_n": q_n,
                        "q_no2": q_no2,
                        "q_difference": q_difference,
                    }
            )

    with OUTPUT_CSV.open("w", newline="")as f:
        writer = csv.DictWriter(
                f,
                fieldnames=rows[0].keys(),
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {OUTPUT_CSV}")

    for row in rows:
        print(
                f"CN={row['cn']:.1f}  "
                f"CCl={row['ccl']:.1f}  "
                f"BO_CN={row['bo_cn']:.3f}  "
                f"BO_CCl={row['bo_ccl']:.3f}  "
                f"BO_diff={row['bo_difference']:+.3f}  "
                f"qCl={row['q_cl']:+.3f}  "
                f"qNO2={row['q_no2']:+.3f}  "
                f"q_diff={row['q_difference']:+.3f}"
        )


if __name__ == "__main__":
    main()



