import numpy as np
from pathlib import Path


def set_atom_distance(coords, anchor_i, moving_i, distance):
    """
    Move atom `moving_i` along the line from anchor_i to moving_i.
    Indices are 1-based.
    """
    coords = coords.copy()

    a = coords[anchor_i - 1]
    b = coords[moving_i - 1]

    direction = b - a
    direction = direction / np.linalg.norm(direction)

    coords[moving_i - 1] = a + direction * distance

    return coords

def move_group_to_distance(coords, anchor_i, group_indices, reference_i, distance):
    """
    Move a whole group so that reference_i is at chosen distance from anchor_i.
    Indices are 1-based.
    """
    coords = coords.copy()

    anchor = coords[anchor_i - 1]
    ref = coords[reference_i - 1]

    direction = ref - anchor
    direction = direction / np.linalg.norm(direction)

    new_ref = anchor + direction * distance
    shift = new_ref - ref

    for i in group_indices:
        coords[i - 1] += shift

    return coords

def generate_cn_ccl_scan(
    start_xyz,
    out_dir,
    ipso_i,
    nitro_group=(1, 2, 3),
    nitro_n_i=3,
    cl_i=25,
    cn_values=None,
    ccl_values=None,
):
    atoms, coords0 = read_xyz(start_xyz)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if cn_values is None:
        cn_values = np.linspace(1.45, 3.20, 10)

    if ccl_values is None:
        ccl_values = np.linspace(3.20, 1.80, 10)

    rows = []

    for k, (cn, ccl) in enumerate(zip(cn_values, ccl_values)):
        coords = coords0.copy()

        # move whole NO2 group by N atom
        coords = move_group_to_distance(
            coords,
            anchor_i=ipso_i,
            group_indices=nitro_group,
            reference_i=nitro_n_i,
            distance=cn,
        )

        # move Cl atom
        coords = set_atom_distance(
            coords,
            anchor_i=ipso_i,
            moving_i=cl_i,
            distance=ccl,
        )

        xyz_path = out_dir / f"scan_{k:02d}.xyz"

        write_xyz(
            xyz_path,
            atoms,
            coords,
            comment=f"C-N={cn:.3f} C-Cl={ccl:.3f}",
        )

        rows.append({
            "frame": k,
            "xyz": str(xyz_path),
            "r_CN": cn,
            "r_CCl": ccl,
        })

    return pd.DataFrame(rows)


