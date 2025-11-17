import matplotlib.pyplot as plt
import math
from matplotlib.ticker import StrMethodFormatter



def plot_metrics_bi(result_df):
    plt.figure(figsize=(10,8)) #,layout='constrained'
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    metrics = ['loss', 'auc','acc'] ##
    for n, metric in enumerate(metrics):
        name = metric.replace("_"," ").capitalize()
        plt.subplot(2,2,n+1)
        plt.plot(result_df.index, result_df[metric], color=colors[0], label='Train')
        plt.plot(result_df.index, result_df['dev_'+metric],
                    color=colors[3], linestyle="--", label='Dev')
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        plt.legend()
        plt.xlabel('Epoch')
        plt.ylabel(name)
