# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 21:48:52 2024

@author: adolf01p
"""

#----------------------GOF plots------------------------------------
# loading important packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pharmpy.modeling import *
from pharmpy.tools import fit
import pharmpy
from pharmpy.tools import read_modelfit_results
import re
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error
from statsmodels.nonparametric.smoothers_lowess import lowess
import seaborn as sns

# Choosing working directory
os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final_model_dfs")

#----------------------------------Efficiently extracting GOF plots and residual plots from output table -------------------------------


# Define a mapping of abbreviated names to full names
name_mapping = {
    'DV': 'Measured conc. (mg/L)',
    'PRED': 'Population predicted conc. (mg/L)',
    'WRES': 'Weighted residual',
    'RES': 'Residual',
    'IWRES': 'Individual weighted residual',
    'CWRES': 'Conditional weighted residual',
    'IPRED': 'Individual predicted conc. (mg/L)',
    'TSPD': 'Time since previous dose (hours)',
    'TAD' : 'Time since previous dose (days)',
    'ETA1': 'Eta on CL',
    'ETA2': 'Eta on V1',
    'ETA3': 'Eta on V2',
    'ETA4': 'Eta on VM'
}

# MAPE_dict = {}

def gof_scatter(x_values, y_values, model_name, x_label, y_label, prediction_type):
    # Calculating MAE
    #MAPE = mean_absolute_error(y_true=x_values, y_pred=y_values)
    MAPE = round(np.mean(np.abs((x_values-y_values) / x_values))*100,2)
    # MAPE_dict[model_name, prediction_type]=MAPE
    
    plt.figure(figsize=(8, 8))
    # Add a small constant to avoid zero values
    x_values = x_values + 1e-10
    y_values = y_values + 1e-10

    plt.scatter(x_values, y_values, color='blue', alpha=0.3, s=10, label=f'Observed vs. {y_label}')
    plt.plot([min(x_values.min(), y_values.min()), max(x_values.max(), y_values.max())], 
             [min(x_values.min(), y_values.min()), max(x_values.max(), y_values.max())], 
             color='red', linestyle='--', label='Perfect Fit')
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal')  # Keep aspect ratio equal

    # Add grey lines
    plt.plot([min(x_values.min(), y_values.min()), max(x_values.max(), y_values.max())], 
             [0.5 * min(x_values.min(), y_values.min()), 0.5 * max(x_values.max(), y_values.max())], 
             color='gray', alpha=0.3, linestyle='--', label='Two times more/less than Perfect Fit')
    
    plt.plot([min(x_values.min(), y_values.min()), max(x_values.max(), y_values.max())], 
             [2 * min(x_values.min(), y_values.min()), 2 * max(x_values.max(), y_values.max())], 
             color='gray', alpha=0.3, linestyle='--')

    # Fill the space between the grey lines with light grey
    plt.fill_between([min(x_values.min(), y_values.min()), max(x_values.max(), y_values.max())], 
                     [0.5 * min(x_values.min(), y_values.min()), 0.5 * max(x_values.max(), y_values.max())], 
                     [2 * min(x_values.min(), y_values.min()), 2 * max(x_values.max(), y_values.max())], 
                     color='lightgray', alpha=0.3)

    max_limit = max(x_values.max(), y_values.max())
    min_limit = min(x_values.min(), y_values.min())

    # Set equal range for both axes
    plt.xlim(min_limit, max_limit)
    plt.ylim(min_limit, max_limit)
    
    # Plotting legend and MAE
    plt.legend()
    plt.figtext(0.5, 0.01, 'MAPE = %0.4f %%' % MAPE, wrap=True, horizontalalignment='center', fontsize=12)
    
    # Save the plot as PNG
    file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_thesis\\GOF_{model_name.replace('.', '')}_{prediction_type}.png"
    plt.savefig(file_name)
    
    plt.show()


# Define the directory containing the model files
directory = "C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final_model_dfs"

# Iterate over all files in the directory
for file_name in os.listdir(directory):
    if file_name.endswith(".csv"):
        # Load the output table
        file_path = os.path.join(directory, file_name)
        output_table = pd.read_csv(file_path)
        output_table = output_table[output_table["TAD"] != 0]

        # Extract the model number and name from the file name
        match = re.search(r'final_output_table_(\d+-\d+-\d+)', file_name) #pas dit aan voor andere patronen in naam!!
        if match:
            model_number = match.group(1)
            model_name = f"run{model_number}.mod"

            # Calculate the absolute difference in TIME since the previous dose (hours)
            output_table['TSPD'] = output_table.groupby('ID')['TIME'].diff().fillna(0)

            # Generate GOF Scatter plot of Observations versus Population Predictions
            x_values = output_table['DV']
            y_values = output_table['PRED']
            W = output_table['RES'] / output_table['WRES']

            gof_scatter(x_values, y_values, model_name, name_mapping['DV'], name_mapping['PRED'], 'PRED')

            # Generate GOF Scatter plot of Observations versus Individual Predictions
            x_values = output_table['DV']
            y_values = output_table['IPRED']
            gof_scatter(x_values, y_values, model_name, name_mapping['DV'], name_mapping['IPRED'], 'IPRED')

# MAPE_list = [MAPE_dict.keys(), MAPE_dict.values()]
# MAPE_output = pd.DataFrame(MAPE_list)
# MAPE_output.to_csv("MAPE.csv")
