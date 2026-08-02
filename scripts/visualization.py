from pathlib import Path
import py3Dmol


def save_molecule_html(
    xyz_data,
    out_file,
    width=800,
    height=600,
    style=None,
    add_labels=False
):
    if style is None:
        style = {"stick": {}}

    view = py3Dmol.view(width=width, height=height)
    view.addModel(xyz_data, "xyz")
    view.setStyle(style)

    if (add_labels==True):
        view.addPropertyLabels(
            "index",
            "",
            {"fontSize": 14}
        )

    view.zoomTo()

    html = view._make_html()

    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(html, encoding="utf-8")

    return out_file
