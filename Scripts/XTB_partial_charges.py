import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
import os
import pandas as pd
from rdkit.Chem import AllChem
import autode
from autode import Molecule

#set path variable to database
path_to_database = r'/u/fd/spet5093/ML_boronic_acids/Boronic Acid Database.xlsx'
#path_to_database = r'C:\Users\mtoho\OneDrive\Oxford\Amide Bond Formation\Data\Boronic Acid Database.xlsx'

#import data
raw_data = pd.read_excel(path_to_database)

#Generate list of catalysts, amines and acids
catalysts = [Molecule(smiles=smiles) for smiles in raw_data['Catalyst']]
amines = [Molecule(smiles=smiles) for smiles in raw_data['Amine']]
acids = [Molecule(smiles=smiles) for smiles in raw_data['Acid']]

boron_charges = []
carbonyl_charges = []
nitrogen_charges = []

#Calculate partial charges for boron of catalyst

for catalyst in catalysts:
    calc = autode.Calculation(name=catalyst.name, method=autode.methods.XTB(), keywords=autode.Config.XTB.keywords.sp, molecule = catalyst)
    calc.run()
    atomic_charges= calc.get_atomic_charges()
    boron_charge = next(charge for i, charge in enumerate(atomic_charges) if catalyst.atoms[i].label == 'B')
    boron_charges.append(boron_charge)

#Calculate partial charges for nitrogen of amine

boron_charges = np.asarray(boron_charges)
np.savetxt('boron_partial_charges_xtb.csv', boron_charges, delimiter=",")