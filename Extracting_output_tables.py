# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:16:30 2023

@author: adolf01p
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
os.chdir("C:\\Users\\adolf01p\\__pycache__")
from pharmpy.modeling import *
from pharmpy.tools import fit
import pharmpy
from pharmpy.tools import read_modelfit_results

os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\pirana_model_develop")

# Set this variable based on how you're running the script
Ran_from_PIRANA = True  # Set to True if running from Pirana, False if running from Python

# Run NONMEM model from python 
model_name = "run300-7-0.mod" # Can be altered manually to specify model name
model_number = int(model_name[3:6]) #Depends on model name style

# Defining path to output tables
if not Ran_from_PIRANA:
    # Picking latest ran model from python
    modelfit_dirs = [int(d.split("modelfit_dir")[-1]) for d in os.listdir() if d.startswith("modelfit_dir")]
    latest_modelfit_dir = max(modelfit_dirs, default=0)
    
    # Define path to output tables
    sdtab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\modelfit_dir{latest_modelfit_dir}\\models\\run{model_number:03d}\\sdtab{model_number:03d}"
    cotab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\modelfit_dir{latest_modelfit_dir}\\models\\run{model_number:03d}\\cotab{model_number:03d}"
    patab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\modelfit_dir{latest_modelfit_dir}\\models\\run{model_number:03d}\\patab{model_number:03d}"

else:
    # Or picking from Pirana ran models
    latest_modelfit_dir = None
    
    # Define path to output tables
    sdtab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\sdtab{model_number:03d}"
    cotab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\cotab{model_number:03d}"
    patab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final Pirana\\patab{model_number:03d}"

# Function to read table files
def read_table_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    column_names = lines[1].strip().split()
    data = [line.strip().split() for line in lines[2:]]
    return pd.DataFrame(data, columns=column_names)

# Read sdtab, patab, and cotab files into DataFrames
sdtab_df = read_table_file(sdtab_path)
patab_df = read_table_file(patab_path)
cotab_df = read_table_file(cotab_path)

# Concatenate cotab_df and patab_df to sdtab_df
merged_df = pd.concat([sdtab_df, cotab_df, patab_df], axis=1)

# Drop duplicate columns
merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]

# Change LENGTH column name to HEIGHT
merged_df = merged_df.rename(columns={'LENGTH': 'HEIGHT'})

# Define columns to change to int and float
int_columns = ['ID', 'TIME', 'TAD', 'MDV', 'EVID', 'AGE', 'SEX']
float_columns = ['DV', 'IPRED', 'IPREDH', 'IPREDL', 'IWRES', 'CWRES', 
                 'SNR', 'PRED', 'RES', 'WRES', 'WT', 'HEIGHT', 'CL', 'V1', 'V2', 'VM', 
                 'ETA1', 'ETA2', 'ETA3', 'ETA4']

# Convert scientific notation to float
for col in float_columns:
    merged_df[col] = merged_df[col].apply(lambda x: float(x))

# Convert to int
for col in int_columns:
    merged_df[col] = merged_df[col].apply(lambda x: int(float(x)))

#----------------Exporting the table-----------------------------
# Define the path for the output directory
output_dir = "C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final_model_dfs"


def export_to_csv(merged_df, model_number, output_dir):
    
    # Define the output file path with model number and timestamp
    output_file = os.path.join(output_dir, f"final_output_table_{model_number:03d}.csv")
    
    # Export the merged_df to the CSV file
    merged_df.to_csv(output_file, index=False)

# Export merged_df
export_to_csv(merged_df, model_number, output_dir)


#-----------------------Exporting all PIRANA ran models from the past----------------------------------------------------------------

import os
import glob
import re
import pandas as pd

# Define the pattern to match the model files
pattern = os.path.join("run[0-9][0-9][0-9]-[0-9]-[0-9].mod") #alter for different file formatting

# Get a list of all model files
model_files = glob.glob(pattern)
model_files = "run300-7-0.mod"
# Extract the model number and process each model
for model_name in model_files:
    # Extract model number from the filename using regular expression
    match = re.search(r'run(\d+-\d+-\d+)\.mod', model_name) # Alter for different file formatting
    if match:
        model_number = match.group(1)
        
        # Define paths to output tables
        sdtab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\pirana_model_develop\\sdtab{model_number}"
        cotab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\pirana_model_develop\\cotab{model_number}"
        patab_path = f"C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\pirana_model_develop\\patab{model_number}"

        # Initialize DataFrames as None
        sdtab_df = None
        cotab_df = None
        patab_df = None

        # Function to read table files
        def read_table_file(file_path):
            if os.path.exists(file_path):
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                column_names = lines[1].strip().split()
                data = [line.strip().split() for line in lines[2:]]
                return pd.DataFrame(data, columns=column_names)
            else:
                return None

        # Read sdtab, patab, and cotab files into DataFrames
        sdtab_df = read_table_file(sdtab_path)
        cotab_df = read_table_file(cotab_path)
        patab_df = read_table_file(patab_path)

        # Create a list of DataFrames that are not None
        available_dfs = [df for df in [sdtab_df, cotab_df, patab_df] if df is not None]

        if available_dfs:
            # Concatenate DataFrames
            merged_df = pd.concat(available_dfs, axis=1)

            # Drop duplicate columns
            merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]

            # Change LENGTH column name to HEIGHT
            merged_df = merged_df.rename(columns={'LENGTH': 'HEIGHT'})

            # Define columns to change to int and float
            int_columns = ['ID', 'TIME', 'TAD', 'MDV', 'EVID', 'SMPL', 'AGE', 'SEX']
            float_columns = ['TSFD', 'DV', 'IPRED', 'IPREDH', 'IPREDL', 'IWRES', 'CWRES', 
                            'SNR', 'PRED', 'RES', 'WRES', 'WT', 'HEIGHT', 'CL', 'V1', 'V2', 'VM', 
                            'ETA1', 'ETA2', 'ETA3', 'ETA4']

            # Convert to int where applicable
            for col in int_columns:
                if col in merged_df.columns:
                    merged_df[col] = merged_df[col].apply(lambda x: int(float(x)))

            # Convert to float where applicable
            for col in float_columns:
                if col in merged_df.columns:
                    merged_df[col] = merged_df[col].apply(lambda x: float(x))

            # Export the merged_df to the CSV file
            output_dir = "C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\Final_model_dfs"
            output_file = os.path.join(output_dir, f"final_output_table_{model_number}.csv")
            merged_df.to_csv(output_file, index=False)
