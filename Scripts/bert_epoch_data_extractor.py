import os
import pandas as pd
import numpy as np

def epoch_data_collater(data_directory, train_data_directory, no_of_folds=5):
    """
    Extracts and compiles all the training metrics for each epoch for each fold, returning a csv file with all the metrics for each fold

    Parameters:
    data_directory (string): an f string of the data directory path with the iterable as a {fold} brace
    train_data_directory (string): an f string of the training data path for each fold, with the iterable as a {fold} brace
    no_of_folds (int): number of folds the model was trained on

    Returns:
    csv file of all metrics for each fold (one per fold)

    """
    for fold in range(1, no_of_folds + 1):

        directory = data_directory + str(fold)
        fold_train_data = train_data_directory + str(fold) + '_train'
        training_data = pd.read_csv(fold_train_data)
        eval_loss = []
        epoch = []
        mse = []
        r2 = []

        for sub_directory in os.listdir(directory):
            if sub_directory.startswith('.'):
                continue
            #print(sub_directory)
            sub_directory_path = os.path.join(directory, sub_directory)
            #print(int(sub_directory_path[-1]))
            epoch.append(int(sub_directory_path[-1])+1)
            for filename in os.listdir(sub_directory_path):
                if filename == 'eval_results.txt':
                    f = os.path.join(sub_directory_path, filename)
                    with open(f, 'r') as file:
                        for line in file:
                            if "eval_loss" in line:
                                eval_loss.append(float(line[12:]))
                            elif "mse" in line:
                                mse.append(float(line[6:]))
                            elif "r2" in line:
                                r2.append(float(line[5:]))
        #print(epoch)
        std_dev = np.std(training_data['Yield /%'])

        rmse = []
        for error in mse:
            rms_error = np.sqrt(error) * std_dev
            rmse.append(rms_error)

        fold_data = pd.DataFrame(list(zip(epoch, eval_loss, rmse, r2)), columns=['epoch', 'eval_loss', 'rmse', 'r2'])
        fold_data = fold_data.sort_values(by='epoch')
        print(fold_data)

        if 'pretrained' in data_directory:
            fold_data.to_csv(f'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/run_1_pretrained_fold_{fold}_training_result', index=False)
        elif 'ft_cv' in data_directory:
            fold_data.to_csv(f'/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/run_1_ft_fold_{fold}_training_result', index=False)
        
#for pretrained model
epoch_data_collater('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/amidation_boronic_acid_run_1_pretrained_cv_fold_',
                    '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_inputs/fold_', 5)

#for finetuned model
epoch_data_collater('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/amidation_boronic_acid_run_1_ft_cv_fold_',
                    '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_inputs/fold_', 5)