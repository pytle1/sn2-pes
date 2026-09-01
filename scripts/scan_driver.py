from pathlib import Path

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
            (2.20, 2.20),
            (2.40, 2.00),
    ]
    
    previous_xyz = start_xyz
    results = []

    for cn, ccl in points:
        result = optimize_point(
                source_xyz=previous_xyz,
                cn_target=cn,
                ccl_target=ccl,
        )
    
        results.append(result)
        previous_xyz = result["xyz"]

    print("\nScan results:")

for result in results:
    print(
            f"C-N={result['cn']:.2f}  "
            f"C-Cl={result['ccl']:.2f}  "
            f"E={result['energy']:.12f} Eh"
    )


if __name__ == "__main__":
    main()
