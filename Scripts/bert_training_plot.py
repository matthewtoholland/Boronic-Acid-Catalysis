import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def get_plot_params(path_to_data):
    """
    Retrieves data to plot from a datafile

    Parameters:
    path_to_data (string): root path to data file

    Returns:
    lists of evaluation parameters (eval_loss, r2 and rmse values), and lists of their standard errors
    """
    data = pd.read_csv(path_to_data, index_col=0)
    y_eval_loss = list(data.loc["mean eval_loss"])
    y_r2 = list(data.loc["mean r2"])
    y_rmse = list(data.loc["mean RMSE/ %"])
    sterr_eval_loss = list(data.loc["eval_loss sterr"])
    sterr_r2 = list(data.loc["r2 sterr"])
    sterr_rmse = list(data.loc["RMSE sterr / %"])
    return y_eval_loss, y_r2, y_rmse, sterr_eval_loss, sterr_r2, sterr_rmse


fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111)
ax1 = fig.add_subplot(121)
ax2= fig.add_subplot(122)

x = list(range(1,11))

sns.set()
sns.set_style(style='white')
sns.set_palette('deep')

#set parameters for main axes
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_xlabel("Epoch")

#Get data
pt_eval_loss, pt_r2, pt_rmse, pt_sterr_eval_loss, pt_sterr_r2, pt_sterr_rmse = get_plot_params('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/mean_results_per_epoch_pretrained.csv')
ft_eval_loss, ft_r2, ft_rmse, ft_sterr_eval_loss, ft_sterr_r2, ft_sterr_rmse = get_plot_params('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Data/Transformer_outputs/mean_results_per_epoch_finetuned.csv')

#RMSE plots
ax1.errorbar(x, pt_rmse, fmt="o",yerr = pt_sterr_rmse, label = "Pre-trained", markersize=6)
ax1.errorbar(x, ft_rmse, fmt="o", yerr = ft_sterr_rmse, label = "Finetuned", markersize=6)
ax1.set_ylabel('RMSE %')
ax1.set_ylim([0, 40])

#R2 plots
ax2.errorbar(x, pt_r2, fmt='o', yerr=pt_sterr_r2, label='Pre-trained', markersize=6)
ax2.errorbar(x, ft_r2, fmt='o', yerr=ft_sterr_r2, label='Finetuned', markersize=6)
ax2.set_ylim([-0.5, 1.0])
ax2.set_ylabel('$R^{2}$')
ax2.legend()


plt.tight_layout()
plt.savefig('/Users/matthewholland/OneDrive/Oxford/Amide Bond Formation/Visualisation/Figures/yield_bert_training_performance.pdf')
