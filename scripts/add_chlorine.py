import paths
import make_structures
import visualization

make_structures.make_reactant_complex(
        paths.XYZ_DIR / "4nitrobiphenyl.xyz",
        paths.XYZ_DIR / "initial_complex.xyz",
        side = "above",
)

'''
def save_py3dmol_html(xyz_data, out_file, width=800, height=600):
    from pathlib import Path
    import py3Dmol

    if isinstance(xyz_data, Path):
        xyz_data = xyz_data.read_text()

    view = py3Dmol.view(width=width, height=height)
    view.addModel(xyz_data, "xyz")

    view.setStyle({"stick": {}})
    view.setStyle(
        {"elem": "Cl"},
        {
            "sphere": {"scale": 0.45, "color": "green"},
            "stick": {"radius": 0.2, "color": "green"},
        },
    )

    view.zoomTo()

    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(view._make_html(), encoding="utf-8")

    return out_file
'''
# visualization.save_py3dmol_html(
 #       (paths.XYZ_DIR / "initial_complex.xyz").read_text(),
 #       paths.VIEW_DIR / "initial_complex.html",
#)

visualization.save_molecule_html(
        (paths.XYZ_DIR / "initial_complex.xyz").read_text(),
        paths.VIEW_DIR / "initial_complex.html",
        style={'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}},
        add_labels=True
)
