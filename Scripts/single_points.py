import os
import glob
import pandas as pd


all_filenames = [i for i in glob.glob('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/single_points/*')]

combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])

combined_csv.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/all_single_points.csv', index=False, header=['index', 'amine_single_point','acid_single_point','product_single_point','delta_g_reaction'])
