# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 21:53:41 2024

@author: adolf01p
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

#---------------------- MAPE for covariate model analysis in PRED (population predictions)--------------------------------

# Change the working directory
os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Thesis")

# Read the data from Excel file
MAPE_file = pd.read_excel("Covariate_MAPE_pop.xlsx")

# Plot MAPE values for each mAb
plt.figure(figsize=(10, 6))
plt.plot(MAPE_file['Covariate Model'], MAPE_file['Natalizumab'], label='Natalizumab', marker='o')
plt.plot(MAPE_file['Covariate Model'], MAPE_file['Tocilizumab'], label='Tocilizumab', marker='o')
plt.plot(MAPE_file['Covariate Model'], MAPE_file['Guselkumab'], label='Guselkumab', marker='o')

# Plot the sum of MAPE for all three mAbs with a dashed red line
total_mape = (MAPE_file['Natalizumab'] + MAPE_file['Tocilizumab'] + MAPE_file['Guselkumab']) / 3
plt.plot(MAPE_file['Covariate Model'], total_mape, label='Average MAPE', linestyle='--', color='red')

# Change x-axis ticks and labels
plt.xticks(MAPE_file['Covariate Model'], ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'])

# Add horizontal grid lines in the background
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Add labels and legend
plt.xlabel('Covariate Model Number')
plt.ylabel('MAPE (%)')
plt.legend()

# Show the plot
plt.show()