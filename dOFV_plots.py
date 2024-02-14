# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:05:51 2023

@author: adolf01p
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import statistics as st
os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Thesis")


#-----------------Parameter base models-----------------------------------------------

df = pd.read_excel("dOFV_param.xlsx")

# Calculate dOFV for each model
base_model = 'B0'
for model in range(1, 9):
    df[f'dOFV_Model{model}'] = df[f'B{model}'] - df[base_model]

# Create the plot
plt.figure(figsize=(10, 6))

# Plot dOFV for Natalizumab in first row
plt.plot(df.columns[2:10], df.loc[0, 'dOFV_Model1':'dOFV_Model8'], label=df.loc[0, 'mAb'], marker='x')

# Plot dOFV for Tocilizumab in second row
plt.plot(df.columns[2:10], df.loc[1, 'dOFV_Model1':'dOFV_Model8'], label=df.loc[1, 'mAb'], marker='x')

# Plot dOFV for Guselkumab in third row
plt.plot(df.columns[2:10], df.loc[2, 'dOFV_Model1':'dOFV_Model8'], label=df.loc[2, 'mAb'], marker='x')

# # Plot the sum of dOFV for each model
# plt.plot(df.columns[2:10], df.iloc[:, 10:].sum(axis=0), label='Sum of dOFV', linestyle='--')

# Add a dashed zero line
plt.axhline(0, color='gray', linestyle='--')

# Add labels and title
plt.ylabel('dOFV')

# Add legend
plt.legend()

# Show the plot
plt.show()




#---------------------------------Covariate models------------------------------------------


df = pd.read_excel("dOFV_cov.xlsx")

# Calculate dOFV for each model
base_model = 'C0'
for model in range(1, 10):
    df[f'dOFV_Model{model}'] = df[f'C{model}'] - df[base_model]

# Create the plot
plt.figure(figsize=(10, 6))

# Plot dOFV for Natalizumab
plt.plot(df.columns[2:11], df.loc[0, 'dOFV_Model1':'dOFV_Model9'], label=df.loc[0, 'mAb'], marker='x')

# Plot dOFV for Tocilizumab
plt.plot(df.columns[2:11], df.loc[1, 'dOFV_Model1':'dOFV_Model9'], label=df.loc[1, 'mAb'], marker='x')

# Plot dOFV for Guselkumab
plt.plot(df.columns[2:11], df.loc[2, 'dOFV_Model1':'dOFV_Model9'], label=df.loc[2, 'mAb'], marker='x')

# Plot the sum of dOFV for each model
plt.plot(df.columns[2:11], df.iloc[:, 11:].sum(axis=0), label='Sum of dOFV', linestyle='--')

# Add a dashed zero line
plt.axhline(0, color='gray', linestyle='--')

# Add labels and title
plt.ylabel('dOFV')

# Add legend
plt.legend()

# Show the plot
plt.show()

