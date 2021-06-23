import pandas as pd
import numpy as np


def find_homo_lumo(molecule, no_instances=298, homo=True):
    molecule_homo_lumo = []
    molecule_fmo = []
    molecule_dipole = []
    for i in range(0, no_instances):
        with open(f'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/second_attempt/{molecule}{i}_sp_orca.out') as f:
            contents = f.readlines()

            energies_flag = False
            for idx, line in enumerate(contents):
                if 'Magnitude (Debye)' in line:
                    dipole_line = line.split()
                    dipole = float(dipole_line[3])
                    molecule_dipole.append(dipole)
                if line == '  NO   OCC          E(Eh)            E(eV) \n':
                    energies_flag = True
                    continue
                if energies_flag == True:
                    entry = list(map(float, line.split()))

                    if entry[1] == 0 and list(map(float, contents[idx-1].split()))[1] == 2:
                        molecule_homo_lumo.append(entry[3]-list(map(float, contents[idx-1].split()))[3])
                        if homo == True:
                            molecule_fmo.append(list(map(float, contents[idx-1].split()))[3])
                            energies_flag=False
                        elif homo == False:
                            molecule_fmo.append(entry[3])
                            energies_flag=False
                        #break
    return molecule_homo_lumo, molecule_fmo, molecule_dipole

amine_homo_lumo, amine_homo, amine_dipole = find_homo_lumo('Amine', homo=True)
acid_homo_lumo, acid_lumo, acid_dipole = find_homo_lumo('Acid', homo=False)
product_homo_lumo, product_homo, product_dipole = find_homo_lumo('Product')

acid_amine_gap = [acid_lumo[i] - amine_homo[i] for i in range(0, 298)]

data = pd.DataFrame(list(zip(amine_homo_lumo, amine_homo, amine_dipole, acid_homo_lumo, acid_lumo, acid_dipole, product_homo_lumo, product_dipole, acid_amine_gap)),
          columns= ['amine_homo_lumo', 'amine_homo', 'amine_dipole', 'acid_homo_lumo', 'acid_lumo','acid_dipole', 'product_homo_lumo','product_dipole', 'acid_amine_gap'])

data.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/all_homo_lumos.csv', index=False)

