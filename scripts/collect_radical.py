from pathlib import Path
import csv
import re


PROJECT = Path(__file__).resolve().parents[1]
RADICAL = PROJECT / "calculations" / "radical"

OUTPUT = RADICAL / "radical_summary.csv"


def get_energy(text):
    matches = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
        text,
    )

    if not matches:
        return None

    return float(matches[-1])


def get_frequencies(text):
    pattern = re.compile(
        r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1",
        re.MULTILINE,
    )

    return [
        float(value)
        for value in pattern.findall(text)
    ]


def get_hirshfeld(text):
    marker = "HIRSHFELD ANALYSIS"

    if marker not in text:
        return {}

    section = text.split(marker, 1)[1]

    atoms = {}

    pattern = re.compile(
        r"^\s*(\d+)\s+([A-Za-z]+)\s+"
        r"([+-]?\d+\.\d+)\s+"
        r"([+-]?\d+\.\d+)",
        re.MULTILINE,
    )

    for match in pattern.finditer(section):
        index = int(match.group(1))

        atoms[index] = {
            "element": match.group(2),
            "charge": float(match.group(3)),
            "spin": float(match.group(4)),
        }

    return atoms


def xyz_distance(xyz_file, i, j):
    lines = xyz_file.read_text().splitlines()[2:]

    atoms = []

    for line in lines:
        fields = line.split()

        atoms.append(
            (
                fields[0],
                float(fields[1]),
                float(fields[2]),
                float(fields[3]),
            )
        )

    _, xi, yi, zi = atoms[i]
    _, xj, yj, zj = atoms[j]

    return (
        (xi - xj) ** 2
        + (yi - yj) ** 2
        + (zi - zj) ** 2
    ) ** 0.5


def main():
    rows = []

    for directory in sorted(RADICAL.iterdir()):
        analysis = directory / "freq_spin" / "analysis.out"
        xyz = directory / "freq_spin" / "geometry.xyz"

        if not analysis.exists():
            continue

        text = analysis.read_text()

        frequencies = get_frequencies(text)

        # Ignore projected translational/rotational zeros.
        vibrational = [
            f for f in frequencies
            if abs(f) > 1.0
        ]

        negative = [
            f for f in vibrational
            if f < 0.0
        ]

        hirshfeld = get_hirshfeld(text)

        # Atom numbering from our structure:
        # O = 0,1
        # N = 2
        # ipso C = 9
        # Cl = 24

        q_no2 = sum(
            hirshfeld[i]["charge"]
            for i in (0, 1, 2)
        )

        spin_no2 = sum(
            hirshfeld[i]["spin"]
            for i in (0, 1, 2)
        )

        rows.append(
            {
                "species": directory.name,
                "energy": get_energy(text),

                "cn": xyz_distance(xyz, 2, 9),
                "ccl": xyz_distance(xyz, 9, 24),

                "negative_frequencies": len(negative),
                "lowest_frequency": min(vibrational)
                if vibrational else None,

                "q_cl": hirshfeld[24]["charge"],
                "spin_cl": hirshfeld[24]["spin"],

                "q_cipso": hirshfeld[9]["charge"],
                "spin_cipso": hirshfeld[9]["spin"],

                "q_no2": q_no2,
                "spin_no2": spin_no2,
            }
        )

    if not rows:
        raise RuntimeError("No radical analyses found")

    with OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=rows[0].keys(),
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {OUTPUT}\n")

    for row in rows:
        print(
            f"{row['species']}: "
            f"CN={row['cn']:.3f} "
            f"CCl={row['ccl']:.3f}  "
            f"imag={row['negative_frequencies']}  "
            f"lowest={row['lowest_frequency']:.1f} cm^-1  "
            f"spin(Cl)={row['spin_cl']:+.3f}  "
            f"spin(C)={row['spin_cipso']:+.3f}  "
            f"spin(NO2)={row['spin_no2']:+.3f}"
        )


if __name__ == "__main__":
    main()
