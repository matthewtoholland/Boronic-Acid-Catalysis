import os
import pandas as pd
import numpy as np
import scipy
from scipy import stats

eval_loss_all_folds = []
r2_all_folds = []
rmse_all_folds = []

def metric_list_appender(epoch, training_output_df, eval_loss_list, r2_list, rmse_list):
        eval_loss_list.append(training_output_df["eval_loss"][epoch])
        r2_list.append(training_output_df["r2"][epoch])
        rmse_list.append(training_output_df["rmse"][epoch])

def averager(metric_list):
    mean_metric_per_epoc = list([np.mean(x) for x in zip(*metric_list)])
    sem_per_epoc = list([scipy.stats.sem(x) for x in zip(*metric_list)])
    return mean_metric_per_epoc, sem_per_epoc

def iterator(fold, parent_path):
    eval_loss_per_fold = []
    r2_per_fold = []
    rmse_per_fold = []

    path = parent_path + str(fold) + '_training_result'
    training_output_df = pd.read_csv(path)
    epochs = training_output_df.epoch.to_list()
    training_output_df.index = epochs

    for epoch in epochs:
        metric_list_appender(epoch, training_output_df, eval_loss_per_fold, r2_per_fold, rmse_per_fold)

    eval_loss_all_folds.append(eval_loss_per_fold)
    r2_all_folds.append(r2_per_fold)
    rmse_all_folds.append(rmse_per_fold)


#Pretrained Model
for i in range(1,6):
    iterator(i, '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/run_1_pretrained_fold_')

mean_eval_loss_per_epoch, sem_eval_loss_per_epoch = averager(eval_loss_all_folds)
mean_r2_per_epoch, sem_r2_per_epoch = averager(r2_all_folds)
mean_rmse_per_epoch, sem_rmse_per_epoch = averager(rmse_all_folds)

results_df = pd.DataFrame([mean_eval_loss_per_epoch, sem_eval_loss_per_epoch, mean_r2_per_epoch, sem_r2_per_epoch, mean_rmse_per_epoch, sem_rmse_per_epoch], columns=["epoch_1", "epoch_2", "epoch_3", "epoch_4", "epoch_5", "epoch_6", "epoch_7", "epoch_8", "epoch_9", "epoch_10"], index = ["mean eval_loss", "eval_loss sterr", "mean r2", "r2 sterr", "mean RMSE/ %", "RMSE sterr / %"])
results_df.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/mean_results_per_epoch_pretrained.csv')

#finetuned model
for i in range(1, 6):
    iterator(i, '/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/run_1_ft_fold_')

mean_eval_loss_per_epoch, sem_eval_loss_per_epoch = averager(eval_loss_all_folds)
mean_r2_per_epoch, sem_r2_per_epoch = averager(r2_all_folds)
mean_rmse_per_epoch, sem_rmse_per_epoch = averager(rmse_all_folds)

results_df = pd.DataFrame([mean_eval_loss_per_epoch, sem_eval_loss_per_epoch, mean_r2_per_epoch, sem_r2_per_epoch, mean_rmse_per_epoch, sem_rmse_per_epoch], columns=["epoch_1", "epoch_2", "epoch_3", "epoch_4", "epoch_5", "epoch_6", "epoch_7", "epoch_8", "epoch_9", "epoch_10"], index = ["mean eval_loss", "eval_loss sterr", "mean r2", "r2 sterr", "mean RMSE/ %", "RMSE sterr / %"])
results_df.to_csv('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/mean_results_per_epoch_finetuned.csv')