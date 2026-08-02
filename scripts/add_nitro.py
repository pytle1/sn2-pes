import paths
import make_structures
import visualization

visualization.save_molecule_html(
        (paths.XYZ_DIR / "4chlorobiphenyl.xyz").read_text(),
        paths.VIEW_DIR / "4chlorobiphenyl.html",
        style={'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}},
        add_labels=True
)

make_structures.add_floating_no2(
        paths.XYZ_DIR / "4chlorobiphenyl.xyz",
        paths.XYZ_DIR / "final_complex.xyz",
)

visualization.save_molecule_html(
        (paths.XYZ_DIR / "final_complex.xyz").read_text(),
        paths.VIEW_DIR / "final_complex.html",
        style={'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}},
        add_labels=True
)
