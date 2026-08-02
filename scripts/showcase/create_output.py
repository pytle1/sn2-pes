# >> get output object for current calculation
output = calc.get_output()

# > parse results from ORCA's JSON files
output.parse()

# > point group of the molecule
output.gbw_results[0].molecule.pointgroup

# > total run time
output.results_properties.calculation_timings.sum

# > graph
output.print_graph()

# > getter to get final energy
output.get_final_energy()

# > Grepper

from opi.output.grepper.core import Grepper

# > get output
output = calc.get_output()
