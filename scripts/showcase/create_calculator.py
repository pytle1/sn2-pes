from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, SimpleKeyword
from opi.input.blocks import BlockGeom

# >> initializing calculator
calc = Calculator(basename="job")
# >> addition of simple keyword
calc.input.add_simple_keywords(BasisSet.DEF2_SVP)
# >> addition of blocks with a respective block option
calc.input.add_blocks(BlockGeom(maxiter=500))
# >> memory per core
calc.input.maxcore = 1000
# >> number of CPU cores
calc.input.ncores=4
# >> initial molecular orbital guess
calc.input.moinp = "initial_guess.gbw"
# >> addition of arbitrary strings
calc.input.add_arbitrary_string("# arbitrary string", pos="top")
# >> write ORCA input file
calc.write_input()
# >> run the ORCA calculation
calc.run()
