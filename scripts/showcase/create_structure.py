from opi.input.structures import Structure
from rdkit import Chem
from rdkit.Chem import AllChem

# from XYZ file
water_xyz = Structure.from_xyz("water_xyz.xyz", charge=0, multiplicity=0)

# from SMILES string
water_smiles = Structure.from_smiles("O", charge=0, multiplicity=0)

# HF from lists
HF_list = Structure.from_lists(["F", "H"],
                               [(0, 0, 0), (0, 0, 0.95)],
                               charge=0, multiplicity=1)

# from RDKit's ROMol
mol = Chem.MolFromSmiles("CCO")
# add hydrogens
mol = Chem.AddHs(mol)
# generate 3D coordinates
AllChem.EmbedMMolecule(mol)
cco_rdkit = Structure.from_rdkitmol(mol, charge=0, multiplicity=1)

