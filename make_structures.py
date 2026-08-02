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


def make_product_complex(
    xyz_path,
    out_path,
    ipso_index=10,
    nitro_n_index=3,
    c_cl_distance=1.75,
    no2_distance=3.2,
):
    atoms, coords = read_xyz(xyz_path)

    c = coords[ipso_index - 1]
    n = coords[nitro_n_index - 1]

    direction = n - c
    direction = direction / np.linalg.norm(direction)

    # remove NO2 atoms: O1, O2, N3
    remove = {0, 1, 2}
    atoms_core = [a for i, a in enumerate(atoms) if i not in remove]
    coords_core = np.array([r for i, r in enumerate(coords) if i not in remove])

    cl = c + direction * c_cl_distance

    # keep NO2- as a group translated away
    no2_atoms = [atoms[i] for i in [0, 1, 2]]
    no2_coords = coords[[0, 1, 2]]

    old_n = coords[2]
    new_n = c + direction * no2_distance
    shift = new_n - old_n
    no2_coords = no2_coords + shift

    atoms2 = atoms_core + ["Cl"] + no2_atoms
    coords2 = np.vstack([coords_core, cl, no2_coords])

    write_xyz(out_path, atoms2, coords2, "4-chlorobiphenyl + NO2- product complex")
