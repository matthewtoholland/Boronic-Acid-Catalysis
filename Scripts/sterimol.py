from rdkit import Chem
from rdkit.Chem import AllChem
from morfeus import Sterimol, read_xyz
import pandas as pd

product_data = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/boronic_acid_with_products.csv')

amine_atom_index = []
acid_atom_index = []
catalyst_atom_index = []

def get_partner(atom, molecule, neighbour='C'):
  for atom in molecule.GetAtomWithIdx(atom).GetNeighbors():
    if atom.GetSymbol() == neighbour:
      partner = atom.GetIdx()
      return partner

for idx in product_data.index:
  amine = Chem.AddHs(Chem.MolFromSmiles(product_data['Amine'][idx]))
  acid = Chem.AddHs(Chem.MolFromSmiles(product_data['Acid'][idx]))
  catalyst = Chem.AddHs(Chem.MolFromSmiles(product_data['Catalyst'][idx]))

  amine_partners = []
  #Check to see if primary amine
  if len(amine.GetSubstructMatch(Chem.MolFromSmarts('[N;H2]'))) != 0:
    nitrogen = amine.GetSubstructMatch(Chem.MolFromSmarts('[N;H2]'))[0]
    amine_partners.append(nitrogen)

    #Determine the index of the connection to the rest of the molecule
    for atom in amine.GetAtomWithIdx(nitrogen).GetNeighbors():
      amine_partners.append(atom.GetIdx())
    amine_atom_index.append(amine_partners)

  #If not primary amine, check for secondary amine (not considering 3ry here)
  elif len(amine.GetSubstructMatch(Chem.MolFromSmarts('[N;H1]'))) != 0:
    #now check that there aren't any amides - these cause problems
    amides = Chem.MolFromSmarts('[CX3]([NX3])=[OX1]')
    if len(amine.GetSubstructMatch(amides)) !=0:
      nitrogen = amine.GetSubstructMatch(Chem.MolFromSmarts('[n;H1]'))[0]
      amine_partners.append(nitrogen)
    else:
      nitrogen = amine.GetSubstructMatch(Chem.MolFromSmarts('[N;H1]'))[0]
      amine_partners.append(nitrogen)

    for atom in amine.GetAtomWithIdx(nitrogen).GetNeighbors():
      amine_partners.append(atom.GetIdx())
    amine_atom_index.append(amine_partners)

  #Determine the indices of the acid moiety
  substruct = acid.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'))
  acid_partners = []
  for atom in substruct:
    if acid.GetAtomWithIdx(atom).GetSymbol() == 'C':
      acid_partners.append(atom)
      acid_partners.append(get_partner(atom,acid, 'C'))
      acid_atom_index.append(acid_partners)

  #Now for the Boron catalysts
  boron_partners = []
  oxygens = []
  boronic_acid_motif = catalyst.GetSubstructMatch(Chem.MolFromSmarts('B(O)(O)'))
  for atom in boronic_acid_motif:
    if catalyst.GetAtomWithIdx(atom).GetSymbol() == 'B':
      boron_partners.append(atom)
  for atom in catalyst.GetAtomWithIdx(boron_partners[0]).GetNeighbors():
    boron_partners.append(atom.GetIdx())
  catalyst_atom_index.append(boron_partners)

sterimols = pd.DataFrame(columns=['amine_L_1', 'amine_B1_1', 'amine_B5_1', 'amine_L_2', 'amine_B1_2', 'amine_B5_2', 'amine_L_3', 'amine_B1_3', 'amine_B5_3', 'acid_L', 'acid_B1', 'acid_B5','catalyst_L_1', 'catalyst_B1_1', 'catalyst_B5_1',
                                    'catalyst_L_2', 'catalyst_B1_2', 'catalyst_B5_2', 'catalyst_L_3', 'catalyst_B1_3', 'catalyst_B5_3'])


def calculate_sterimol(index, molecule,atom1, atom2):
    elements, coordinates = read_xyz(f'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/xyz_files/{molecule}_{index}.xyz')
    sterimol = Sterimol(elements, coordinates, atom1, atom2)
    return sterimol.L_value_uncorrected, sterimol.B_1_value, sterimol.B_5_value


for idx in product_data.index:
    sterimol_dict = {}
    for i in range(1, len(amine_atom_index[idx])):
        sterimol_dict[f'amine_L_{i}'], sterimol_dict[f'amine_B1_{i}'], sterimol_dict[f'amine_B5_{i}'] = calculate_sterimol(idx, 'amine', amine_atom_index[idx][0] +1, amine_atom_index[idx][i]+1)

    for i in range(1, len(catalyst_atom_index[idx])):
        sterimol_dict[f'catalyst_L_{i}'], sterimol_dict[f'catalyst_B1_{i}'], sterimol_dict[f'catalyst_B5_{i}'] = calculate_sterimol(idx, 'catalyst', catalyst_atom_index[idx][0]+1, catalyst_atom_index[idx][i]+1)
    
    sterimol_dict['acid_L'], sterimol_dict['acid_B1'], sterimol_dict['acid_B5'] = calculate_sterimol(idx, 'acid', acid_atom_index[idx][0]+1, acid_atom_index[idx][1]+1)

    sterimols = sterimols.append(sterimol_dict, ignore_index=True)

sterimols.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/DFT_features/sterimols.csv', index=False)
