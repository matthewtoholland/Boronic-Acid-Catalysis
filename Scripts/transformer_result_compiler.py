#Compiles the results from each fold validation from the yield-Bert transformer

import os
import pandas as pd
import numpy as np


def compile_results(path_to_evaluation, path_to_training_data, no_of_folds=5):
    """
    extracts the evaluation loss, rmse, and r2 values from the final validation step of each fold in the yield-BERT transformer

    Parameters:
    path_to_evaluation (string): path to evaluation file outputted by transformer - iterable (usually fold) must be replaced with {}
    path_to_training_data (string): path to dataset used to train the fold (not validation data) - iterable replaced as above
    no_of_folds (int): number of folds used in the cross-validation

    Returns:
    pandas DataFrame containing the root mean squared error, r2 correlation and the evaluation loss for each fold

    """
    eval_loss = []
    rmse = []
    r2 = []


    for fold in range(1,no_of_folds + 1):
        eval_file = path_to_evaluation.format(fold)
        training_data = pd.read_csv(path_to_training_data.format(fold), index_col=0)
        standard_deviation = np.std(training_data['Yield /%'])

        with open(eval_file, 'r') as file:
            for line in file:
                if 'eval_loss' in line:
                    eval_loss.append(float(line[12:]))
                elif 'r2' in line:
                    r2.append(float(line[5:]))
                elif 'mse' in line:
                    mse = float(line[6:])
                    rmse_one = np.sqrt(mse)
                    rmse_one = rmse_one * standard_deviation
                    rmse.append(rmse_one)

    return pd.DataFrame(list(zip(rmse, r2, eval_loss)), columns=['rmse', 'r2', 'eval_loss'])

pretrained_output = compile_results('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/amidation_boronic_acid_run_1_ft_cv_fold_{}_eval_val_df/eval_results.txt', 
                                    '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_inputs/fold_{}_train', 5)
pretrained_output.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/pretrained_transformer_output.csv', index=False)

finetuned_output = compile_results('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/amidation_boronic_acid_run_1_ft_cv_fold_{}_eval_val_df/eval_results.txt',
                                    '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_inputs/fold_{}_train', 5)
finetuned_output.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/finetuned_transformer_outputs.csv', index=False)
