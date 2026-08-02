from pathlib import Path
import os
from opi.core import Calculator
import py3Dmol

print(Path.cwd())

# PROJECT = Path(__file__).resolve().parents[1]
PROJECT = Path.cwd()

SAVE_DIR = PROJECT / "results"
FIG_DIR = PROJECT / "figures"
CALC_DIR = PROJECT / "calculations"
VIEW_DIR = FIG_DIR / "molecules"
XYZ_DIR = PROJECT / "molecules"

for d in [SAVE_DIR, FIG_DIR, CALC_DIR, VIEW_DIR, XYZ_DIR]:
    d.mkdir(parents=True, exist_ok=True)


folder = CALC_DIR / "easy_water2"
print(folder)

calc = Calculator(basename="job",
                  working_dir=Path(folder))

output = calc.get_output()
output.parse()

energy = output.get_final_energy()

print(energy)

graph = output.print_graph()

pointgroup = output.gbw_results[0].molecule.pointgroup
print(pointgroup)

cube_data = output.plot_mo(mo_index).cube

# > py3Dmol viewer
view = py3Dmol.view(width=500, height=400)
view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})

# > add cube data to py3Dmol viewer
settings = {"isoval": 0.005, "color": "blue", "opacity":0.8}
view.AddVolumetricData(cube_data, "cube", settings)
settings = {"isoval": 0.005, "color": "red", "opacity":0.8}
view.AddVolumetricData(cube_data, "cube", settings)

# > zoom in and save
view.zoomTo()

name = "easy_water2_orbitals“
html = view._make_html()
out = VIEW_DIR / name3
out.write_text(html, encoding="utf-8")

view.save(
