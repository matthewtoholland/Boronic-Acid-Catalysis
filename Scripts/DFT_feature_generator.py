import pandas as pd
import numpy as np

#import the various datasets
sterimols = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/DFT_features/sterimols.csv')
single_points = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/second_attempt/single_points/single_points/single_point_energies.csv', usecols=[1, 2, 3,4], header=None, names=['amine_single_point', 'acid_single_point', 'product_single_point', 'delta_g_reaction'])
homo_lumo = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/all_homo_lumos.csv')
boron_partial_charges = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/boron_partial_charges_xtb.csv', header=None, names=['boron_partial_charge'])

homo_lumo = homo_lumo.drop(['amine_homo', 'acid_lumo'], axis=1)
#print(boron_partial_charges)
#Create one megadatabase comprised of all the other datasets
all_features = pd.concat([sterimols, single_points, homo_lumo, boron_partial_charges], axis = 1)

all_features.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/DFT_features/all_dft_features.csv', index=False)