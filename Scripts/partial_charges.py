import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
import os
import pandas as pd
from rdkit.Chem import AllChem
from brokenaxes import brokenaxes

#set path variable to database
path_to_database = r'C:\Users\mtoho\OneDrive\Oxford\Amide Bond Formation\Data\Boronic Acid Database.xlsx'

#import data
raw_data = pd.read_excel(path_to_database)

#generate catalyst list of rd_mol objects
catalyst_list = [Chem.MolFromSmiles(smiles) for smiles in raw_data['Catalyst']]

#Compute Gasteiger partial charges
for catalyst in catalyst_list:
    AllChem.ComputeGasteigerCharges(catalyst)

#Generate list of Boron Gasteiger charges
boron_partial_charges = []
for catalyst in catalyst_list:
    for atom in catalyst.GetAtoms():
        if atom.GetSymbol() == 'B':
            boron_partial_charges.append(atom.GetDoubleProp('_GasteigerCharge'))



#create histogram of partial charges
fig = plt.figure()
ax1, ax2 = fig.add_subplot(2, 1)

ax1.hist(boron_partial_charges)
ax2.hist(boron_partial_charges)

plt.show()