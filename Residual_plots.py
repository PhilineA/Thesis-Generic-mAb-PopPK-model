# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 11:33:08 2023

@author: adolf01p
"""
# loading important packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
os.chdir("C:\\Users\\adolf01p\\__pycache__")
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


def eta_cov(x_values, y_values, model_name, x_label, y_label):
    # Calculating correlation coefficient and p-value
    corr = pearsonr(x_values, y_values)
    
    # Calculate LOESS fit
    loess_smoothed = lowess(y_values, x_values, frac=0.9)
    
    # Plotting figure
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, color='blue', alpha=0.3)
    plt.plot(loess_smoothed[:, 0], loess_smoothed[:, 1], color='red', linewidth=2, label='LOESS')
    plt.axhline(y=0, color='grey', linestyle='--')
    plt.axvline(x=0, color='grey', linestyle='--')
    
    # Check if the labels are in the name_mapping dictionary
    x_label_full = name_mapping.get(x_label, x_label)
    y_label_full = name_mapping.get(y_label, y_label)
    
    plt.xlabel(x_label_full)  # Use the full name for the x-axis label
    plt.ylabel(y_label_full)  # Use the full name for the y-axis label
    
    # Adding correlation coefficient and corresponding p-value to the plot
    plt.figtext(0.5, 0.01, 'r = %0.5f p-value = %0.4f' % corr, wrap=True, horizontalalignment='center', fontsize=12)
    
    # Save the plot as PNG using abbreviations in the file name
    file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\eta_cov_{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
    plt.savefig(file_name)
    
    plt.legend()  # Add legend for LOESS curve
    plt.show()

def eta_biplot(x_values, y_values, model_name, x_label, y_label):
    # Calculating correlation coefficient and p-value
    corr = pearsonr(x_values, y_values)
    
    # Plotting figuer
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, color='blue', alpha=0.3)
    plt.axhline(y=0, color='grey', linestyle='--')
    plt.axvline(x=0, color='grey', linestyle='--')
    # Check if the labels are in the name_mapping dictionary
    x_label_full = name_mapping.get(x_label, x_label)
    y_label_full = name_mapping.get(y_label, y_label)
    
    plt.xlabel(x_label_full)  # Use the full name for the x-axis label
    plt.ylabel(y_label_full)  # Use the full name for the y-axis label
    
    # Adding correlation coefficient and corresponding p-value to the plot
    plt.figtext(0.5, 0.01, 'r = %0.5f p-value = %0.4f' % corr, wrap=True, horizontalalignment='center', fontsize=12)
    
    # Save the plot as PNG using abbreviations in the file name
    file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\eta_biplot_{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
    plt.savefig(file_name)
    
    plt.show()
    
def gof_scatter(x_values, y_values, model_name, x_label, y_label, prediction_type):
    # Calculating MAE
    MAE = mean_absolute_error(y_true=x_values, y_pred=y_values)

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
    plt.figtext(0.5, 0.01, 'MAE = %0.4f' % MAE, wrap=True, horizontalalignment='center', fontsize=12)
    
    # Save the plot as PNG
    file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_thesis\\GOF_{model_name.replace('.', '')}_{prediction_type}.png"
    plt.savefig(file_name)
    
    plt.show()


# Residual plot function
def residual_plot(x_values, y_values, model_name, x_label, y_label):
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, alpha=0.3)
    plt.axhline(y=0, color='grey', linestyle='--')
    
    # Check if the labels are in the name_mapping dictionary
    x_label_full = name_mapping.get(x_label, x_label)
    y_label_full = name_mapping.get(y_label, y_label)
    
    plt.xlabel(x_label_full)  # Use the full name for the x-axis label
    plt.ylabel(y_label_full)  # Use the full name for the y-axis label
    
    # Save the plot as PNG using abbreviations in the file name
    file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\residual_{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
    plt.savefig(file_name)
    
    plt.show()

# Load the output table
file_name = "final_output_table_202_15_35_32_20231221.csv"
output_table = pd.read_csv(file_name)
# Deleting drug administration data rows
output_table = output_table[output_table["TAD"] != 0]

# Extract the model number and name from the file name
model_number = re.search(r'final_output_table_(\d{3})-(\d{1})', file_name).group(1)
model_name = f"run{model_number}.mod"

#------------------------GOF plots-----------------------------------------------

# Generate GOF Scatter plot of Observations versus Population Predictions
x_values = output_table['DV']
y_values = output_table['PRED']
W = output_table['RES']/output_table['WRES']

gof_scatter(x_values, y_values, model_name, name_mapping['DV'], name_mapping['PRED'], 'PRED')

# Generate GOF Scatter plot of Observations versus Individual Predictions
x_values = output_table['DV']
y_values = output_table['IPRED']
gof_scatter(x_values, y_values, model_name, name_mapping['DV'], name_mapping['IPRED'], 'IPRED')

#----------------------Residual plots-------------------------------------------
# Plotting CWRES against PRED
x_values = output_table['PRED']
y_values = output_table['CWRES']
residual_plot(x_values, y_values, model_name, 'PRED', 'CWRES')

# Plotting RES against PRED
x_values = output_table['PRED']
y_values = output_table['RES']
residual_plot(x_values, y_values, model_name, 'PRED', 'RES')

# Plotting WRES against PRED
x_values = output_table['PRED']
y_values = output_table['WRES']
residual_plot(x_values, y_values, model_name,'PRED', 'WRES')

# Plotting WRES against Time Since Previous Dose
x_values = output_table['TAD']
y_values = output_table['WRES']
residual_plot(x_values, y_values, model_name,'TAD', 'WRES')

# -----------------------Plotting eta bi-plots--------------------------------
# Iff r>= 0.8 --> shared-eta approach

#Eta CL vs Eta V1
x_values = output_table['ETA2']
y_values = output_table['ETA1']
eta_biplot(x_values, y_values, model_name, 'ETA2', 'ETA1')

# Eta CL vs Eta V2
x_values = output_table['ETA3']
y_values = output_table['ETA1']
eta_biplot(x_values, y_values, model_name, 'ETA3', 'ETA1')

# Eta CL vs Eta Vm
x_values = output_table['ETA4']
y_values = output_table['ETA1']
eta_biplot(x_values, y_values, model_name, 'ETA4', 'ETA1')

# Eta V1 vs Eta V2
x_values = output_table['ETA3']
y_values = output_table['ETA2']
eta_biplot(x_values, y_values, model_name, 'ETA3', 'ETA2')

# Eta V1 vs Eta Vm
x_values = output_table['ETA4']
y_values = output_table['ETA2']
eta_biplot(x_values, y_values, model_name, 'ETA4', 'ETA2')

# Eta Vm vs Eta V2
x_values = output_table['ETA3']
y_values = output_table['ETA4']
eta_biplot(x_values, y_values, model_name, 'ETA3', 'ETA4')


#_-----------------Eta histogram----------------
# Extracting 'ETA2' values
eta2_values = output_table['ETA2']

# Creating a histogram for ETA2 values
plt.figure(figsize=(8, 6))
plt.hist(eta2_values, bins=20, edgecolor='black')  # Adjust the number of bins as needed

# Adding labels and title
plt.xlabel(name_mapping.get('ETA2', 'ETA2'))
plt.ylabel('Frequency')

# Show the plot
plt.grid(True)
plt.tight_layout()
plt.show()

#--------------Covariance Eta testing------------------
#Eta CL vs WT
x_values = output_table['WT']
y_values = output_table['ETA1']
eta_cov(x_values, y_values, model_name, 'WT', 'ETA1')

#Eta CL vs NWT
output_table["NWT"] = (output_table["WT"]-output_table["WT"].min()) / (output_table["WT"].max() - output_table["WT"].min())
x_values = output_table['NWT']
y_values = output_table['ETA1']
eta_biplot(x_values, y_values, model_name, 'NWT', 'ETA1')

# --------------- Plotting Eta distributions------------------------

# Get a list of all columns with names starting with 'ETA'
eta_columns = [col for col in output_table.columns if col.startswith('ETA')]

# Create subplots
fig, axes = plt.subplots(nrows=len(eta_columns), figsize=(10, 6 * len(eta_columns)))

# Define a colormap
colors = plt.cm.viridis.colors

# Loop through each ETA column
for i, eta_col in enumerate(eta_columns):
    groups_eta = output_table.groupby('ID')[eta_col].first()
    groups_eta = pd.DataFrame(groups_eta)
    
    # Plot density with separate color
    groups_eta.plot.density(color=colors[i], ax=axes[i])
    axes[i].set_xlabel(eta_col)

plt.tight_layout()
plt.show()

# Save the plot as PNG using abbreviations in the file name
file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\eta_distr_{model_name.replace('.', '').replace('run','')}.png"
plt.savefig(file_name)

# -------------------------- PRED/IPRED vs Time---------------------------------------
# PRED vs Time
plt.figure(figsize=(10, 6))

x_values = output_table['TAD']
y_values = output_table['PRED']

plt.scatter(x_values, y_values, alpha=0.3)
# Check if the labels are in the name_mapping dictionary
x_label = 'TAD'
y_label = 'PRED'
x_label_full = name_mapping.get(x_label, x_label)
y_label_full = name_mapping.get(y_label, y_label)

plt.xlabel(x_label_full)  # Use the full name for the x-axis label
plt.ylabel(y_label_full)  # Use the full name for the y-axis label

# Save the plot as PNG using abbreviations in the file name
file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
plt.savefig(file_name)
 
plt.show()


# IPRED vs Time
plt.figure(figsize=(10, 6))

x_values = output_table['TAD']
y_values = output_table['IPRED']

plt.scatter(x_values, y_values, alpha=0.3)
# Check if the labels are in the name_mapping dictionary
x_label = 'TAD'
y_label = 'IPRED'
x_label_full = name_mapping.get(x_label, x_label)
y_label_full = name_mapping.get(y_label, y_label)

plt.xlabel(x_label_full)  # Use the full name for the x-axis label
plt.ylabel(y_label_full)  # Use the full name for the y-axis label

# Save the plot as PNG using abbreviations in the file name
file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
plt.savefig(file_name)
 
plt.show()

# # ----------Plotting DV against TAD--------------------------------------------------------------
plt.figure(figsize=(10, 6))
sub_data = output_table[output_table["TAD"]!= 0]
x_values = sub_data['TAD']
y_values = sub_data['DV']
plt.scatter(x_values, y_values, alpha=0.3)
# Check if the labels are in the name_mapping dictionary
x_label = 'TAD'
y_label = 'DV'
x_label_full = name_mapping.get(x_label, x_label)
y_label_full = name_mapping.get(y_label, y_label)

plt.xlabel(x_label_full)  # Use the full name for the x-axis label
plt.ylabel(y_label_full)  # Use the full name for the y-axis label

# Save the plot as PNG using abbreviations in the file name
file_name = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\residual_{y_label}_vs_{x_label}_{model_name.replace('.', '').replace('run','')}.png"
plt.savefig(file_name)
 
plt.show()

# ----------------- Comparing Eta in different datasets ---------------------------------
file_name1 = "final_output_table_062_11_16_19_20231020.csv"
file_name2 = "final_output_table_063_16_47_32_20231020.csv"
basemod_fix_nata = "final_output_table_064_11_16_19_20231020.csv"
basemod_fix_toci = "final_output_table_065_16_47_32_20231020.csv"
basemod_nata = "final_output_table_060_11_16_19_20231020.csv"
basemod_toci = "final_output_table_061_16_47_32_20231020.csv"

Natalizumab_cov = pd.read_csv(file_name1)
Tocilizumab_cov = pd.read_csv(file_name2)
Nata_base_fix = pd.read_csv(basemod_fix_nata)
Toci_base_fix = pd.read_csv(basemod_fix_toci)
Nata_base = pd.read_csv(basemod_nata)
Toci_base = pd.read_csv(basemod_toci)

Natalizumab = Natalizumab_cov[Natalizumab_cov["TAD"] != 0]
Tocilizumab = Tocilizumab_cov[Tocilizumab_cov["TAD"] != 0]
Nata_base_fix = Nata_base_fix[Nata_base_fix["TAD"] != 0]
Toci_base_fix = Toci_base_fix[Toci_base_fix["TAD"] != 0]
Nata_base = Nata_base[Nata_base["TAD"] != 0]
Toci_base = Toci_base[Toci_base["TAD"] != 0]

# Combining the data
combined_data = pd.concat([Natalizumab[['ETA1']], Nata_base_fix[['ETA1']], Nata_base[['ETA1']], Tocilizumab[['ETA1']], Toci_base_fix[['ETA1']], Toci_base[['ETA1']]], axis=1)
combined_data.columns = ['Natalizumab + cov WT_CL', 'Natalizumab + fixed parameters', 'Natalizumab base', 'Tocilizumab + cov WT_CL', 'Tocilizumab + fixed parameters', 'Tocilizumab base']

# Create a figure
plt.figure(figsize=(10, 6))

# Set the color palette to 'flare'
sns.set_palette('flare')

# Create a box plot
sns.boxplot(data=combined_data)

# Add strip plot with adjusted transparency and color
sns.stripplot(data=combined_data, jitter=True, edgecolor='grey', alpha=0.3)  # Adjusted color and alpha values

# Set axis labels
plt.ylabel(r'$\eta$CL')
plt.xlabel('mAbs dataset')

# Rotate x-axis ticks
plt.xticks(rotation=45)  # Adjust the rotation angle as needed

# Show the plot
plt.tight_layout()  # Ensures all elements fit within the figure area

file_name = "C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\GOF_residual_plots\\ETACL_across_different_datasets.png"
plt.savefig(file_name)

plt.show()

# -------------Exporting descriptions of combined data----------
# making the descriptions table
combined_data = combined_data.rename(columns={'Natalizumab + cov WT_CL': 'NTLZ + FIX + Cov', 'Natalizumab + fixed parameters': 'NTLZ + FIX',
       'Natalizumab base': 'NTLZ base', 'Tocilizumab + cov WT_CL': 'TCL + FIX + Cov',
       'Tocilizumab + fixed parameters': 'TCL + FIX', 'Tocilizumab base': 'TCL base'})
description = pd.DataFrame(combined_data.describe())
description.reset_index(inplace=True)
description = description.rename(columns={'index': 'Statistics'})

def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    data_rounded = data.apply(lambda x: np.round(x, 4) if np.issubdtype(x.dtype, np.number) else x)
    mpl_table = ax.table(cellText=data_rounded.values, bbox=bbox, colLabels=data.columns, **kwargs)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
            cell.set_text_props(wrap=True, va='center')  # Wrap text for header cells
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
            cell.set_text_props(wrap=True, va='top')  # Wrap text for non-header cells
    return ax.get_figure(), ax

fig, ax = render_mpl_table(description, header_columns=0, col_width=4.0)  # Adjusted col_width

plt.show()

