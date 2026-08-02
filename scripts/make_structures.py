from pathlib import Path
import numpy as np


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    atoms = []
    coords = []

    for line in lines[2:]:
        parts = line.split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])

    return atoms, np.array(coords)


def write_xyz(path, atoms, coords, comment="generated"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = [str(len(atoms)), comment]
    for atom, (x, y, z) in zip(atoms, coords):
        lines.append(f"{atom:2s} {x:12.6f} {y:12.6f} {z:12.6f}")

    path.write_text("\n".join(lines) + "\n")


def make_reactant_complex(
    xyz_path,
    out_path,
    ipso_index=10,
    distance=3.0,
    side="above",
):
    atoms, coords = read_xyz(xyz_path)

    c = coords[ipso_index - 1]

    if side == "above":
        cl = c + np.array([0.0, 0.0, distance])
    elif side == "below":
        cl = c + np.array([0.0, 0.0, -distance])
    else:
        raise ValueError("side must be 'above' or 'below'")

    atoms2 = atoms + ["Cl"]
    coords2 = np.vstack([coords, cl])

    write_xyz(out_path, atoms2, coords2, "4-nitrobiphenyl + Cl- reactant complex")

def add_floating_no2(
    xyz_path,
    out_path,
    ipso_index=7,
    distance=3.0,
    no_bond=1.25,
    ono_angle_deg=120.0,
    side="above",
):
    atoms, coords = read_xyz(xyz_path)

    c = coords[ipso_index - 1]

    # molecule is almost in xy plane, so z is orthogonal
    if side == "above":
        direction = np.array([0.0, 0.0, 1.0])
    elif side == "below":
        direction = np.array([0.0, 0.0, -1.0])
    else:
        raise ValueError("side must be 'above' or 'below'")

    # N position
    n = c + direction * distance

    # define NO2 plane: use x/y directions around N
    u = np.array([0.0, 0.0, 1.0])
    v = np.array([0.0, 1.0, 0.0])

    half_angle = np.deg2rad(ono_angle_deg / 2)

    o1 = n + no_bond * (np.cos(half_angle) * u + np.sin(half_angle) * v)
    o2 = n + no_bond * (np.cos(half_angle) * u - np.sin(half_angle) * v)

    atoms2 = atoms + ["N", "O", "O"]
    coords2 = np.vstack([coords, n, o1, o2])

    write_xyz(
        out_path,
        atoms2,
        coords2,
        "4-chlorobiphenyl + floating NO2- product complex",
    )
