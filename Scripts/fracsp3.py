import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
import os
import pandas as pd
from rdkit.Chem import AllChem

#set path variable to database
path_to_database = r'C:\Users\mtoho\OneDrive\Oxford\Amide Bond Formation\Data\Boronic Acid Database.xlsx'

#import data
raw_data = pd.read_excel(path_to_database)

#generate amine list
amine_list = [Chem.MolFromSmiles(smiles) for smiles in raw_data['Amine']]
amine_frac_sp3 = [Chem.rdMolDescriptors.CalcFractionCSP3(amine) for amine in amine_list]

#generate acid list
acid_list = [Chem.MolFromSmiles(smiles) for smiles in raw_data['Acid']]
acid_frac_sp3 = [Chem.rdMolDescriptors.CalcFractionCSP3(acid) for acid in acid_list]

#Create plot
fig = plt.figure()

plt.scatter(amine_frac_sp3, acid_frac_sp3)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('Amine Fraction sp3')
plt.ylabel('Acid Fraction sp3')
plt.title('Proportion of sp3 atoms in couplings')
plt.savefig('Figures/frac_sp3')

