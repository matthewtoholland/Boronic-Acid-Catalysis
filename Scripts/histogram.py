import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
import os
import pandas as pd
from rdkit.Chem import AllChem
import seaborn as sns

#set path variable to database
path_to_database = '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Boronic Acid Database.xlsx'
#import data
raw_data = pd.read_excel(path_to_database)

#Generate histogram of yields
fig = plt.figure(1)
ax1 = fig.add_subplot(111)

sns.set()
sns.set_palette('pastel')
sns.set_style(style='white')
ax1.hist(raw_data['Yield /%'], bins=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], alpha=0.5, edgecolor='black')
plt.ylim([0, 70])
plt.xlim([0,100])
plt.xlabel('Yield /%', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
#plt.title('Distribution of Yields')
plt.savefig('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Visualisation/Figures/yield_distribution.png')
#plt.show()
#Now for Solvents
solvent_counts = raw_data['Organic Solvent'].value_counts()

solvents = list(solvent_counts.index.values)
count = list(solvent_counts)
solvents = [solvent.replace('fluorobenzene', 'PhF') for solvent in solvents]

fig = plt.figure(2)

plt.bar(solvents, count)
plt.xlabel('Solvent')
plt.ylabel('Count')
#plt.title('Distribution of Solvents')
plt.savefig('Figures/Solvent_distribution.png')
