import sklearn
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit

#set path to necessary data
acid_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/acid_features.csv'
amine_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/amine_features.csv'
catalyst_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Features/catalyst_features.csv'
data_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Boronic Acid Database.xlsx'
partial_charge_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/boron_partial_charges_xtb.csv'

#Import Data
acid_features = pd.read_csv(acid_path)
amine_features = pd.read_csv(amine_path)
catalyst_features = pd.read_csv(catalyst_path)
partial_charges = pd.read_csv(partial_charge_path, header=None)
database = pd.read_excel(data_path)

#Rename column headers to identify which features belong to which reagent
for column in acid_features:
  acid_features = acid_features.rename(columns={column: column+'_acid'} )

for column in amine_features:
  amine_features = amine_features.rename(columns={column: column+'_amine'})

for column in catalyst_features:
  catalyst_features = catalyst_features.rename(columns={column: column+'_catalyst'})

#Merge partial charges into main database
#database.insert(5, 'Boron_partial_charge', partial_charges)

#Remove non-numerical data from the database dataframe
numeric_database = database.select_dtypes(include='number')

#remove incomplete data - temperature and solvent ratio
numeric_database = numeric_database.drop(columns=['Temperature', 'Ratio (v/v)'])

#combine all numerical features into one enormous database
all_numerical_features = pd.concat([numeric_database, acid_features.iloc[:, 1:], amine_features.iloc[:, 1:], catalyst_features.iloc[:, 1:]], axis=1)

#Generate bins for stratifying by boron partial charge
all_numerical_features['boron_charge_cat'] = pd.cut(all_numerical_features['Boron Partial Charge_catalyst'], bins=[0, 0.3, 0.43, 0.47, 0.5, np.inf], labels=[1, 2, 3, 4, 5])

#Generate a training and a test set stratified by boron partial charge as a proxy for a range of catalysts.
split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
for train_index, test_index in split.split(all_numerical_features, all_numerical_features['boron_charge_cat']):
  strat_train_set = all_numerical_features.loc[train_index]
  strat_test_set = all_numerical_features.loc[test_index]

#Remove the boron charge category labels from the data
for set_ in (strat_train_set, strat_test_set):
  set_.drop('boron_charge_cat', axis=1, inplace=True)

#Save test set to csv and then leave well alone
strat_test_set.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Stratified_test_set.csv')

#Save training set to csv
strat_train_set.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/training_set.csv')
