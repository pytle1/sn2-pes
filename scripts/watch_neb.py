from paths import VIEW_DIR, CALC_DIR
from visualization import save_molecule_html

NEB_VIEW_DIR = VIEW_DIR / "neb_reaction"
NEB_IMAGE_VIEW_DIR = NEB_VIEW_DIR / "images"
NEB_TRAJ_VIEW_DIR = NEB_VIEW_DIR / "trajectory"

for d in [NEB_VIEW_DIR, NEB_IMAGE_VIEW_DIR, NEB_TRAJ_VIEW_DIR]:
    d.mkdir(parents=True, exist_ok=True)

neb_dir = CALC_DIR / "neb_reaction_ordered" / "nudged_elastic_band"
print(neb_dir)

for xyz_path in sorted(neb_dir.glob("ts_im*.xyz")):
    out_file = NEB_IMAGE_VIEW_DIR / f"{xyz_path.stem}.html"

    save_molecule_html(
            xyz_path.read_text(),
            out_file,
            style={'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}},
            add_labels=True
    )

def save_py3dmol_trajectory_html(trajectory_file, out_file, width=800, height=600):
    from pathlib import Path
    import py3Dmol

    trajectory_file = Path(trajectory_file)
    out_file = Path(out_file)
    print("traj file ", trajectory_file)
    print(out_file)

    xyz_data = trajectory_file.read_text()
    print(len(xyz_data))
    print(xyz_data)

    view = py3Dmol.view(width=width, height=height)
    view.addModelsAsFrames(xyz_data, "xyz")
    view.setStyle({}, {"stick": {"radius": 0.15}, "sphere": {"scale": 0.25}})
    view.animate({"interval": 300, "loop": "forward"})
    view.zoomTo()

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(view._make_html(), encoding="utf-8")

    return out_file

save_py3dmol_trajectory_html(
    neb_dir / "ts_MEP_trj.xyz",
    NEB_TRAJ_VIEW_DIR / "ts_MEP_trj.html",
)
