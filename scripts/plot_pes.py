from pathlib import Path
import csv

import matplotlib.pyplot as plt


PROJECT = Path(__file__).resolve().parents[1]

CSV_FILE = PROJECT / "calculations" / "relaxed_scan.csv"
FIGURES = PROJECT / "figures"


def main():
    q = []
    energies = []
    labels = []
    
    with CSV_FILE.open() as f:
        reader = csv.DictReader(f)

        for row in reader:
            cn = float(row["cn"])
            ccl = float(row["ccl"])
            energy = float(row["relative_energy"])

            q.append(cn - ccl)
            energies.append(energy)
            labels.append(row["point"])

    fig, ax = plt.subplots()

    ax.plot(q, energies, marker="o")

    for x, y, label in zip(q, energies, labels):
        ax.annotate(
                label,
                (x, y),
                xytext=(4, 4),
                textcoords="offset points",
        )

    ax.axhline(0, linewidth=0.8)

    ax.set_xlabel(r"$r_{\mathrm{C-N}} - r_{\mathrm{C-Cl}}$ / Å")
    ax.set_ylabel(r"Relative energy / kJ mol$^{-1}$")

    FIGURES.mkdir(exist_ok=True)

    output = FIGURES / "relaxed_scan_1d.png"

    fig.tight_layout()
    fig.savefig(output, dpi=200)

    print(f"Saved {output}")

    plt.show()


if __name__ == "__main__":
    main()
