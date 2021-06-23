from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

path_to_data = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Boronic Acid Database.xlsx'

data = pd.read_excel(path_to_data)

#Blank column for products
data['Product'] = ''

#Create products
for ind in data.index:
  amine = Chem.MolFromSmiles(data['Amine'][ind])
  acid = Chem.MolFromSmiles(data['Acid'][ind])

  #test to see if primary amine
  if amine.HasSubstructMatch(Chem.MolFromSmarts('[N;H2]')) == True:
    rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N;H2:3]>>[C:1](=[O:2])[N:3]')
    ps = rxn.RunReactants((acid,amine))
    data['Product'][ind] = Chem.MolToSmiles(ps[0][0])
    #if no primary amine, react at secondary amine
  elif amine.HasSubstructMatch(Chem.MolFromSmarts('[N;H1]')) == True:
    rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N;H1:3]>>[C:1](=[O:2])[N:3]')
    ps = rxn.RunReactants((acid,amine))
    data['Product'][ind] = Chem.MolToSmiles(ps[0][0])
  else:
      #if there's no primary or secondary amine, you've probably made a mistake
    print(error)
    
data.to_csv('Data/boronic_acid_with_products.csv')

