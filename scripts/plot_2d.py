from pathlib import Path
import csv

import numpy as np
import matplotlib.pyplot as plt


PROJECT = Path(__file__).resolve().parents[1]
CSV_FILE = PROJECT / "calculations" / "scan_2d.csv"
FIGURES = PROJECT / "figures"


def main():
    rows = []

    with CSV_FILE.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                    {
                        "cn": float(row["cn"]),
                        "ccl": float(row["ccl"]),
                        "relative_energy": float(row["relative_energy"]),
                    }
            )

    cn_values = sorted({row["cn"] for row in rows})
    ccl_values = sorted({row["ccl"] for row in rows})

    cn_index = {value: i for i, value in enumerate(cn_values)}
    ccl_index = {value: i for i, value in enumerate(ccl_values)}

    Z = np.full((len(ccl_values), len(cn_values)), np.nan)

    for row in rows:
        i = ccl_index[row["ccl"]]
        j = cn_index[row["cn"]]
        Z[i, j] = row["relative_energy"]


    print(rows)

    X, Y = np.meshgrid(cn_values, ccl_values)

    print("rows:", len(rows))
    print("cn_values:", cn_values)
    print("ccl_values:", ccl_values)
    print("Z:")
    print(Z)
    print("NaNs:", np.isnan(Z).sum())
    print("min/max:", np.nanmin(Z), np.nanmax(Z))

    fig, ax = plt.subplots()

    filled = ax.contourf(X, Y, Z, levels=15)
    contours = ax.contour(X, Y, Z, levels=15, linewidths=0.8)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.0f")

    fig.colorbar(filled, ax=ax, label="Relative energy / kJ mol$^{-1}$")

    ax.set_xlabel(r"$r_{\mathrm{C-N}}$ / Å")
    ax.set_ylabel(r"$r_{\mathrm{C-Cl}}$ / Å")
    ax.set_title("2D relaxed PES around the barrier region")

    FIGURES.mkdir(exist_ok=True)

    output = FIGURES / "scan_2d_contour.png"
    fig.tight_layout()
    fig.savefig(output, dpi=200)

    print(f"Saved {output}")

    plt.show()


if __name__ == "__main__":
    main()
