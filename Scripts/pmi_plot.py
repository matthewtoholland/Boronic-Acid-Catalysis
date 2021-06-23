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
path_to_database = '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Boronic Acid Database.xlsx'

#import data
raw_data = pd.read_excel(path_to_database)

def conformer(mol):

    """
    Takes a rdkit mol object and generates 300 conformers

    mol -- rdkit.Chem.Mol object

    """
    m_string = 'ETKDGv3' if hasattr(AllChem, 'ETKDGv3') else 'ETKDGv2'
    method_class = getattr(AllChem, m_string)
    method = method_class()
    method.pruneRmsThresh = 0.1
    method.numThreads = 6
    method.useSmallRingTorsion = True
    conf_ids = list(AllChem.EmbedMultipleConfs(mol,
                                                numConfs=300,
                                                params=method))

#generate amine list of rdkit mol objects and generate conformers
amine_list = [Chem.MolFromSmiles(smiles) for smiles in raw_data['Amine']]

for amine in amine_list:
   conformer(amine)
    
#generate list of NPR1 and NPR2 for the amines
npr1_amine = [Descriptors3D.NPR1(mol_object) for mol_object in amine_list]
npr2_amine = [Descriptors3D.NPR2(mol_object) for mol_object in amine_list]

#The same as above, for the acids
acid_list = [Chem.MolFromSmiles(smiles) for smiles in raw_data['Acid']]

for acid in acid_list:
    conformer(acid)

npr1_acid = [Descriptors3D.NPR1(mol_object) for mol_object in acid_list]
npr2_acid = [Descriptors3D.NPR2(mol_object) for mol_object in acid_list]

#generate the PMI plot
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(npr1_amine, npr2_amine, c='b', marker ='o', label='Amines', s=5)
ax1.scatter(npr1_acid, npr2_acid, c='r', marker='o', label='Acids', s=5)

#add the triangle to the plot
ax1.plot([0.5, 0], [0.5, 1], c='grey')
ax1.plot([0.5, 1], [0.5, 1], c='grey')
ax1.spines['right'].set_visible(False)

plt.xlabel('NPR1')
plt.ylabel('NPR2')
plt.legend(loc='upper right')
#plt.title('PMI Plot of Acid and Amine Chemical Space Distribution')
plt.xlim([0, 1.0])
plt.ylim([0.5, 1.0])
plt.savefig('Figures/pmi_plot.png')