"""
Generates a set of 20 random features for the training instances

Returns .csv file of 30 random features + labels for the training set instances
"""

import numpy as np
import pandas as pd

#Import data
training_instances = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/training_set.csv', usecols=[0])
training_labels = pd.read_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/training_set.csv', index_col=0, usecols=[0, 4])

random_features = pd.DataFrame(np.random.normal(loc=0, scale=1, size=(238, 30)))

training_data = pd.concat([training_instances, random_features], axis=1)
training_data.set_index('Unnamed: 0', inplace=True)
training_data = pd.concat([training_data, training_labels], axis=1)

training_data.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/random_features/random_features.csv', index_label='Index')