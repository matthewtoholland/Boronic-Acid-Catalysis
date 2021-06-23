import autode as ade
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from autode import Molecule
from autode.methods import XTB, ORCA
import sys

row_index = sys.argv[1]

#set path to file
path_to_database = r'/u/fd/spet5093/ML_boronic_acids/Boronic Acid Database.xlsx'

data = pd.read_excel(path_to_database)


amine = Molecule(smiles=data['Amine'][row_index])
acid = Molecule(smiles=data['Acid'][row_index])
product = Molecule(smiles=data['Product'][row_index])

for mol in (amine, acid, product):
    mol.optimise(method=ORCA())
    mol.single_point(method=ORCA())

amine.print_xyz_file(filename='/u/fd/spet5093/ML_boronic_acids/Energies/xyz_files/amine_{}.xyz'.format{row_index})
acid.print_xyz_file(filename='/u/fd/spet5093/ML_boronic_acids/Energies/xyz_files/acid_{}.xyz'.format{row_index})
product.print_xyz_file(filename='/u/fd/spet5093/ML_boronic_acids/Energies/xyz_files/product_{}.xyz'.format{row_index})

#Calculate delta g
delta_g = (product.energy) - (amine.energy + acid.energy)

energy = pd.DataFrame([amine.energy, acid.energy, product.energy, delta_g], 
                        columns=['amine_single_point', 'acid_single_point', 'product_single_point', 'delta_g_reaction'])

energies.to_csv(r'/u/fd/spet5093/ML_boronic_acids/Energies/single_points{}.csv'.format(row_index))
