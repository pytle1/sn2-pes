from pathlib import Path
import re

from optimize_point import optimize_point


PROJECT = Path(__file__).resolve().parents[1]
CALCULATIONS = PROJECT / "calculations"

CN_VALUES = [round(1.9 + 0.1 * i, 1) for i in range(7)]
CCL_VALUES = [round(1.9 + 0.1 * i, 1) for i in range(7)]

POINT_RE = re.compile(
        r"relaxed_point_cn([0-9.]+)_ccl([0-9.]+)"
)


def existing_points():
    points = {}

    for directory in CALCULATIONS.glob(
            "relaxed_point_cn*_ccl*"
    ):
        match = POINT_RE.fullmatch(directory.name)

        if match is None:
            continue

        xyz = directory / "job.xyz"
        output = directory / "job.out"

        if not xyz.exists() or not output.exists():
            continue
        
        cn = round(float(match.group(1)), 3)
        ccl = round(float(match.group(2)), 3)

        points[(cn, ccl)] = xyz

    return points


def coordinate_distance(a, b):
    cn1, ccl1 = a
    cn2, ccl2 = b

    return (
            (cn1 - cn2) ** 2
            + (ccl1 - ccl2) ** 2
    ) ** 0.5


def nearest_point(target, points):
    return min(
            points,
            key=lambda point: coordinate_distance(
                target,
                point,
            ),
    )


def main():
    calculated = existing_points()

    targets = [
            (cn, ccl)
            for cn in CN_VALUES
            for ccl in CCL_VALUES
    ]
    
    remaining = [
            point
            for point in targets
            if point not in calculated
    ]

    while remaining:
        target = min(
                remaining,
                key=lambda point: min(
                    coordinate_distance(point, done)
                    for done in calculated
                ),
        )
    
        seed = nearest_point(target, calculated)

        source_xyz = calculated[seed]

        cn, ccl = target

        print()
        print(
                f"=== C-N {cn:.2f} / C-Cl {ccl:.2f} ==="
        )
        print(
                f"Starting from "
                f"C-N {seed[0]:.2f} / "
                f"C-Cl {seed[1]:.2f}"
        )

        result = optimize_point(
                source_xyz=source_xyz,
                cn_target=cn,
                ccl_target=ccl,
        )

        calculated[target] = Path(result["xyz"])
        remaining.remove(target)


if __name__ == "__main__":
    main()
