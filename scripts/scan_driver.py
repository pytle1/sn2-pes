from pathlib import Path
import csv

from optimize_point import optimize_point

PROJECT = Path(__file__).resolve().parents[1]

def main():
    start_xyz = (
            PROJECT
            / "calculations"
            / "relaxed_point_cn2.000_ccl2.400"
            / "job.xyz"
    )
    
    points = [
            (2.60, 1.90),
            (2.80, 1.80),
            (3.00, 1.74),
    ]
    
    previous_xyz = start_xyz
    results = []

    for cn, ccl in points:
        print(f"\n=== C-N {cn:.2f} / C-Cl {ccl:.2f} ===")

        try:
            result = optimize_point(
                    source_xyz=previous_xyz,
                    cn_target=cn,
                    ccl_target=ccl,
            )
        except Exception as error:
            print(f"FAILED at C-N={cn:.2f}, C-Cl={ccl:.2f}")
            print(error)
            break

        results.append(result)
        previous_xyz = Path(result["xyz"])

    results_file = PROJECT / "calculations" / "relaxed_scan.csv"

    with results_file.open("w", newline="") as f:
        writer = csv.DictWriter(
                f,
                fieldnames=["cn", "ccl", "energy", "xyz"],
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved results to {results_file}")
    print("\nScan results:")

    for result in results:
        print(
                f"C-N={result['cn']:.2f}  "
                f"C-Cl={result['ccl']:.2f}  "
                f"E={result['energy']:.12f} Eh"
        )


if __name__ == "__main__":
    main()
