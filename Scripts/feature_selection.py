import sklearn
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV

#Set paths to data
strat_train_set_path = r'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/training_set.csv'

#Import data
strat_train_set = pd.read_csv(strat_train_set_path)

#Remove the yield labels from the datasets and create a separate csv
train_features = strat_train_set.drop('Yield /%', axis =1)
train_labels = strat_train_set['Yield /%'].copy()
print(list(train_features.index))
train_labels.reindex(list(train_features.index))
train_labels.to_csv('training_labels.csv')

#Standardise the training features 
scaler = StandardScaler()
train_features_scaled = scaler.fit_transform(train_features)

#convert back to pandas dataframe - just generally more convenient to visualise and work with
column_titles = list(train_features.columns.values)
row_titles = list(train_features.index)
train_features_scaled = pd.DataFrame(train_features_scaled, index=row_titles, columns=column_titles)

#Use Lasso CV to select important features
reg = LassoCV(max_iter=10000)
reg.fit(train_features_scaled, train_labels)
feature_coefficients = pd.Series(reg.coef_, index = train_features_scaled.columns)
reg.get_params()

#Generate a list of the selected features
selected_features =[index for index, value in feature_coefficients.iteritems() if value !=0]
selected_features_values = [feature for feature in feature_coefficients if feature !=0]

#Plot the non-zero co-efficient features on a horizontal bar
imp_coef = feature_coefficients
imp_coef = imp_coef[selected_features].sort_values()
matplotlib.rcParams['figure.figsize'] = (20.0, 20.0)
imp_coef.plot(kind = "barh")
#plt.title("Feature importance using Lasso Model")
plt.savefig('../Visualisation/Figures/lasso_nonzero_figures.png')

#Further refine features by removing features that are below Abs() 0.25
selected_features =[index for index, value in feature_coefficients.iteritems() if value>0.5 or value<-0.5]
imp_coef = feature_coefficients
imp_coef = imp_coef[selected_features].sort_values()

matplotlib.rcParams['figure.figsize'] = (20.0, 20.0)
imp_coef.plot(kind = "barh")
#plt.title("Feature importance using Lasso Model")
plt.xticks(fontsize=14)
plt.xlabel('Lasso Co-Efficient', fontsize=18)
plt.savefig('../Visualisation/Figures/Final_lasso_features.png')

#Generate list of scaled, lasso important features
final_features = train_features_scaled[selected_features]
print(final_features)
final_features.to_csv('scaled_lasso_features.csv')
