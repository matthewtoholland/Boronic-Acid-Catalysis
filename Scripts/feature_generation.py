import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
import os
import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def remove_same_values(dataframe):
    """
    identify and remove columns that only contain a single value

    :param dataframe: dataframe of descriptors
    :type dataframe: pandas dataframe

    :return: pandas dataframe of descriptors with no columns of single value
    :type: pandas dataframe

    """

    same_values = [column for column in dataframe if dataframe[column].nunique() == 1]
    output = dataframe.drop(columns = same_values)
    return output

def generate_features(molecules, output_file=None):
    """
    Computes RD-kit 2D features for given list of smiles strings, and removes columns that only contain a single value

    :param molecules: a dataframe column of smiles strings of input molecules
    :type molecules: pd.dataframe column

    :param output_file: name for csv file containing features
    :type output_file: string

    :return: csv of descriptor name and values for each inputted molecule
    :rtype: csv file

    """
    molecule_list = molecules.to_list()

    features = {}
    descriptors = {d[0]:d[1] for d in Descriptors.descList}

    rdmol_list = [Chem.MolFromSmiles(smiles) for smiles in molecules]

    for rdmol in rdmol_list:
        features[rdmol] = {d: descriptors[d](rdmol) for d in descriptors}

    features = pd.DataFrame.from_dict(features).T
    features = features.rename(index=lambda x: Chem.MolToSmiles(x))

    features = remove_same_values(features)
    
    if output_file !=None: 
        features.to_csv('Features/{}.csv'.format(output_file))
        return 

    else:
        return features

#set path variable to database
path_to_database = r'C:\Users\mtoho\OneDrive\Oxford\Amide Bond Formation\Data\Boronic Acid Database.xlsx'
path_to_partial_charges = r'C:\Users\mtoho\OneDrive\Oxford\Amide Bond Formation\Data\boron_partial_charges_xtb.csv'

#import data
raw_data = pd.read_excel(path_to_database)
boron_partial_charges = pd.read_csv(path_to_partial_charges, header=None, names=['partial_charge'])

#Combine raw data with catalyst partial charges
partial_charges_list = [partial_charge for partial_charge in boron_partial_charges['partial_charge']]


#Generate CSV of RDkit features for amines and acids
generate_features(raw_data['Amine'], 'amine_features')
generate_features(raw_data['Acid'], 'acid_features')

#Generate features for Catalyst, and merge 
catalyst_features = generate_features(raw_data['Catalyst'])
catalyst_features['Boron Partial Charge'] = partial_charges_list
catalyst_features.to_csv('Features/catalyst_features.csv')
