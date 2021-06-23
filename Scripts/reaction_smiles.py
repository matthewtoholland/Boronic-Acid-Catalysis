"""
Generates reaction smiles string from reactants, catalysts and products

Returns: csv of reaction smiles with yield labels
"""

import pandas as pd

database = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/boronic_acid_with_products.csv', index_col=0)

reagent_smiles = [amine + '.' + acid for amine, acid in zip(database['Amine'], database['Acid'])]
reaction_smiles = [reagent_smiles + '>' + catalyst + '>' + product for reagent_smiles, catalyst, product in zip(reagent_smiles,
                                                    database['Catalyst'], database['Product'])]

reaction_smiles = pd.DataFrame(reaction_smiles, columns=['Reaction Smiles'])
reaction_smiles = pd.concat([reaction_smiles, database['Yield /%']], axis=1)
reaction_smiles.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/reaction_smiles.csv', index=False)

