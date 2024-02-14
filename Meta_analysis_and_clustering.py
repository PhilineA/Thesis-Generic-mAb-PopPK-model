# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:53:05 2023

@author: adolf01p
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import statistics as st
import numpy as np
import pandas as pd
from itertools import combinations
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.hierarchy import dendrogram
from sklearn.preprocessing import StandardScaler

os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\mAbs datasets")

df = pd.read_excel("Theo_mAbs_PK(1).xlsx", sheet_name = "Other")

# -------------------- Data preprocessing--------------------------------------

selected_columns= ['Pdf', 'mAb', 'Route', 'ka', 'V1', 'CL', 'Q',
       'V2', 'V_total', 'F', 'V1_BSV %', 'V2_BSV %', 'Q_BSV %', 'CL_BSV %',
       'VM_BSV %', 'KM_BSV %', 'Ka_BSV %', 'F_BSV %', 'Prop_res %',
       'CL/F_BSV %', 'V1/F_BSV %', 'Q/F_BSV %', 'V2/F_BSV %', 'IG type',
       'Patients', 'Special patient group', 'Illness Classification',
       'Autoimmune Disorder', 'Healthy patients included','Model',
       'NonLinCL', 'V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'WT', 'N', 'V1_Sx',
       'CL_Sx', 'V2_Sx', 'SEX', 'V1_Alb', 'CL_Alb', 'V2_Alb', 'ALB', 'CL_ADA',
       'Cov_other', 'V1n', 'CLn', 'Qn', 'V2n', 'Mcount']

data = df[selected_columns]

# Selecting only 2 compartment models
data = data[data["Model"] == 2]

#---------------------------following part is excluded after clustering-------------------------
# Excluding models for children
data = data[data["Special patient group"]!= "Children"]
# Only including models with parallel linear and nonlinear clearance
data = data[data['NonLinCL'] == 1]
#-----------------------------------------------------------------------------------------------

# Only including IV models
desired_routes = ["IV", "IV/SC", "SC/IV"]
IV_SC= data[data["Route"].isin(desired_routes)].reset_index(drop=True)

# Removing models with Nan values for V1n, V2n, Qn and Cln
IV_SC = IV_SC[IV_SC['V2n'].notna()]

# Correcting the Sex covariates for reference Female, which results in Male dependent covariates:
for index, value in enumerate(IV_SC["SEX"]):
    if value == "M":
        theta_old_V1Sx = IV_SC["V1_Sx"].iloc[index]
        theta_old_CLSx = IV_SC["CL_Sx"].iloc[index]
        IV_SC["V1_Sx"].iloc[index] = (1/(1+theta_old_V1Sx)) - 1
        IV_SC["CL_Sx"].iloc[index] = (1/(1+theta_old_CLSx)) - 1

#----------------------Data analysis-------------------------------------------
# ------parameter values-------------
# Group by "mAb" column and calculate mean for each group
mean_V1n_grouped = IV_SC.groupby("mAb")["V1n"].mean()
mean_V2n_grouped = IV_SC.groupby("mAb")["V2n"].mean()
mean_Qn_grouped = IV_SC.groupby("mAb")["Qn"].mean()
mean_CLn_grouped = IV_SC.groupby("mAb")["CLn"].mean()


# Take the overall median of the grouped means
median_V1n = st.median(mean_V1n_grouped.dropna())
median_V2n = st.median(mean_V2n_grouped.dropna())
median_Qn = st.median(mean_Qn_grouped.dropna())
median_CLn = st.median(mean_CLn_grouped.dropna())
median_V1n_V2n = st.median(mean_V1n_grouped + mean_V2n_grouped)

# Plotting distributions of V1 and V2 and Q and CL
plt.figure(figsize=(10, 6))
IV_SC.V1n.plot.density(color='blue', label='V1')
IV_SC.V2n.plot.density(color='orange', label='V2')
IV_SC.Qn.plot.density(color='green', label='Q')
IV_SC.CLn.plot.density(color='red', label='CL')
# IV_SC.V_total.plot.density(color='purple', label='V1 + V2')

plt.legend()
plt.show()

#--------covariate values-----------------------

columns_to_convert = ['V1_Sx', 'CL_Sx', 'ka','F',
        'Q_BSV %', 'VM_BSV %', 'KM_BSV %',
       'Ka_BSV %', 'CL/F_BSV %', 'V1/F_BSV %',
       'Q/F_BSV %', 'V2/F_BSV %']

# Loop through each column and convert to numeric values
for column in columns_to_convert:
    IV_SC[column] = pd.to_numeric(IV_SC[column], errors='coerce')


# Plotting distributions of covariates
# Filter out NA and zero values
IV_SC_CL_WT_filtered = IV_SC.CL_WT.dropna()[IV_SC.CL_WT.dropna() != 0]
IV_SC_V1_WT_filtered = IV_SC.V1_WT.dropna()[IV_SC.V1_WT.dropna() != 0]
IV_SC_V1_Sx_filtered = IV_SC.V1_Sx.dropna()[IV_SC.V1_Sx.dropna() != 0]
IV_SC_CL_Alb_filtered = IV_SC.CL_Alb.dropna()[IV_SC.CL_Alb.dropna() != 0]
IV_SC_CL_Sx_filtered = IV_SC.CL_Sx.dropna()[IV_SC.CL_Sx.dropna() != 0]
IV_SC_V2_WT_filtered = IV_SC.V2_WT.dropna()[IV_SC.V2_WT.dropna() != 0]
IV_SC_Q_WT_filtered = IV_SC.Q_WT.dropna()[IV_SC.Q_WT.dropna() != 0]


# Create the density plots
plt.figure(figsize=(10, 6))
IV_SC_CL_WT_filtered.plot.density(color='orange', label='CL_WT')
IV_SC_V1_WT_filtered.plot.density(color='blue', label='V1_WT')
IV_SC_V1_Sx_filtered.plot.density(color='green', label='V1_Sx')
IV_SC_CL_Alb_filtered.plot.density(color='yellow', label='CL_Alb')
IV_SC_CL_Sx_filtered.plot.density(color='purple', label='CL_Sx')
IV_SC_V2_WT_filtered.plot.density(color='lightblue', label='V2_WT')
IV_SC_Q_WT_filtered.plot.density(color='red', label='Q_WT')

# Add legend and show the plot
plt.legend()
plt.show()

# Define a function to calculate mean excluding zero and NaN values
def mean_without_zeros_and_nans(group):
    return group.replace(0, np.nan).dropna().mean()


mean_CL_WT_grouped = IV_SC.groupby("mAb")["CL_WT"].apply(mean_without_zeros_and_nans)
mean_V1_WT_grouped = IV_SC.groupby("mAb")["V1_WT"].apply(mean_without_zeros_and_nans)
mean_V1_Sx_grouped = IV_SC.groupby("mAb")["V1_Sx"].apply(mean_without_zeros_and_nans)
mean_CL_Alb_grouped = IV_SC.groupby("mAb")["CL_Alb"].apply(mean_without_zeros_and_nans)
mean_CL_Sx_grouped = IV_SC.groupby("mAb")["CL_Sx"].apply(mean_without_zeros_and_nans)
mean_V2_WT_grouped = IV_SC.groupby("mAb")["V2_WT"].apply(mean_without_zeros_and_nans)
mean_Q_WT_grouped = IV_SC.groupby("mAb")["Q_WT"].apply(mean_without_zeros_and_nans)

median_CL_WT= st.median(mean_CL_WT_grouped.dropna())
median_V1_WT = st.median(mean_V1_WT_grouped.dropna())
median_V1_Sx = st.median(mean_V1_Sx_grouped.dropna())
median_CL_Alb = st.median(mean_CL_Alb_grouped.dropna())
median_CL_Sx = st.median(mean_CL_Sx_grouped.dropna())
median_V2_WT = st.median(mean_V2_WT_grouped.dropna())
median_Q_WT = st.median(mean_Q_WT_grouped.dropna())


#-----------showing median values in a df----------
median_values = pd.DataFrame({
    'Parameter': ['V1n', 'V2n', 'Qn', 'CLn', 'CL_WT', 'V1_WT', 'V1_Sx', 
                  'CL_Alb', 'CL_Sx', 'V2_WT', 'Q_WT'],
    'Median': [median_V1n, median_V2n, median_Qn, median_CLn, median_CL_WT, 
               median_V1_WT, median_V1_Sx, median_CL_Alb, median_CL_Sx, median_V2_WT, median_Q_WT]
})

# Display the DataFrame
print(median_values)


#--------------------------------Covariate counts---------------------------------------------

# Replace NaN values with 0 for accurate calculation
covariate_data = IV_SC.copy()
covariate_data.fillna(0, inplace=True)

# List of columns to consider
columns_to_consider = ['V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'V1_Sx', 'CL_Sx', 
                       'V2_Sx', 'V1_Alb', 'CL_Alb', 'V2_Alb', 'CL_ADA', 'Cov_other']

# Define the condition to check for 0 or NaN
condition = (covariate_data[columns_to_consider] == 0) | (covariate_data[columns_to_consider].isna())

# Check if ALL specified columns meet the condition
rows_with_zeros_or_nan = condition.all(axis=1)

# Calculate the number of rows that meet the condition
num_rows_with_zeros_or_nan = rows_with_zeros_or_nan.sum()

print("Number of rows with 0 or NaN in ALL specified columns:", num_rows_with_zeros_or_nan)

# Calculate the number of non-zero and non-NaN values for each column
non_zero_non_nan_counts = covariate_data[columns_to_consider].astype(bool).sum(axis=0)

print(non_zero_non_nan_counts)


# ----- CHecking which combination of covariates is the most present in the models---

#All possible combinations of columns
combinations_list = []
for r in range(1, len(columns_to_consider)+1):
    combinations_list.extend(combinations(columns_to_consider, r))

# Initialize a dictionary to store the counts of each combination
combination_counts = {}

# Iterate through each combination
for combination in combinations_list:
    # Create a subset of the DataFrame with only the selected columns
    selected_columns = list(combination)
    other_columns = [col for col in columns_to_consider if col not in selected_columns]
    subset_df = IV_SC[selected_columns + other_columns]

    # Check if selected columns in the combination have non-zero or non-NaN values,
    # and other columns outside the combination contain ONLY 0 or NaN values
    condition = (
        (subset_df[selected_columns] != 0).all(axis=1) & 
        subset_df[selected_columns].notna().all(axis=1) &
        (subset_df[other_columns].fillna(0) == 0).all(axis=1)
    )

    # Filter rows based on the condition
    filtered_df = subset_df[condition]

    # Count the number of occurrences
    count = len(filtered_df)

    # Store the count in the dictionary
    combination_counts[combination] = count

# Find the top 10 most prevalent combinations
top_15_combinations = sorted(combination_counts.items(), key=lambda x: x[1], reverse=True)[:15]

# Print the top 10 combinations and their counts
for combination, count in top_15_combinations:
    print(f"Combination {combination} has a count of {count}.")

# Create a DataFrame from the top 10 combinations and their counts
top_15_df = pd.DataFrame(top_15_combinations, columns=['Combination', 'Count'])

# ------At least present:---------------------
# find the "At least present combination count" of the corresponding strictly present combinations
# Initialize a dictionary to store the counts of each combination
combination_counts = {}

# Iterate through each combination
for combination in combinations_list:
    # Create a subset of the DataFrame with only the selected columns
    selected_columns = list(combination)
    other_columns = [col for col in covariate_data.columns if col not in selected_columns]
    subset_df = covariate_data[selected_columns + other_columns]

    # Apply the filtering conditions for the selected combination
    condition = (
        (subset_df[selected_columns] != 0).all(axis=1) & 
        subset_df[selected_columns].notna().all(axis=1) 
    )

    # Filter rows based on the condition
    filtered_df = subset_df[condition]

    # Count the number of occurrences
    count = len(filtered_df)

    # Store the count in the dictionary
    combination_counts[combination] = count


# Key combinations to search for in the combination_counts dictionary
key_combinations_to_search = [
    ('V1_WT', 'CL_WT'),
    ('V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'CL_Alb'),
    ('V1_WT', 'CL_WT', 'V1_Sx', 'CL_Sx'),
    ('V1_WT', 'CL_WT', 'V1_Sx', 'CL_Sx', 'V1_Alb', 'CL_Alb'),
    ('CL_WT',),
    ('V1_Sx',),
    ('CL_WT', 'CL_Alb'),
    ('V1_WT', 'CL_WT', 'CL_Alb'),
    ('V1_WT', 'CL_WT', 'V1_Sx', 'CL_Alb'),
    ('V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'V1_Sx',),
    ('V1_WT', 'CL_WT', 'V1_Sx', 'V2_Sx', 'V1_Alb', 'CL_Alb'),
    ('V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'V1_Sx', 'CL_Sx', 'CL_Alb')
    ]

# Search for values of specified key combinations in the combination_counts dictionary
for key_combination in key_combinations_to_search:
    if key_combination in combination_counts:
        print(f"Key Combination {key_combination}: {combination_counts[key_combination]}")
    else:
        print(f"Key Combination {key_combination} not found in combination_counts")

#-------------------------- Group differences------------------------------------

#---------Illness Classification---------------------------------------------------------------------
# -----Calculating the number of Nan count in V1n, V2n, Qn and CLn for each Ilness Class--

# Initialize dictionaries to store counts for each unique illness classification and for NaN values
classification_counts = {}
nan_counts = {'V1n': {}, 'V2n': {}, 'Qn': {}, 'CLn': {}}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    # Handle multiple classifications in a single row
    classifications = row['Illness Classification'].split(', ')
    for classification in classifications:
        classification_cleaned = classification.strip()  # Remove extra spaces
        
        # Add count for individual classifications
        if classification_cleaned in classification_counts:
            classification_counts[classification_cleaned] += 1
        else:
            classification_counts[classification_cleaned] = 1
        
        # Calculate NaN counts for each 'Illness Classification' group and column
        classification_key = '_'.join(classification_cleaned.split())
        for col in ['V1n', 'V2n', 'Qn', 'CLn']:
            if pd.isnull(row[col]):
                if classification_key in nan_counts[col]:
                    nan_counts[col][classification_key] += 1
                else:
                    nan_counts[col][classification_key] = 1

# Display NaN counts for each 'Illness Classification' group and column
for col, counts in nan_counts.items():
    print(f"NaN counts for column '{col}':")
    for classification, nan_count in counts.items():
        print(f"Illness Classification: {classification}, NaN count: {nan_count}")

# If you want to store the counts in a DataFrame:
nan_counts_df = pd.DataFrame(nan_counts)


#-----------Calculating median values for each Illness class--------------

# Initialize dictionaries to store counts for each unique illness classification and for mean/median values
classification_counts = {}
mean_values = {}
median_values = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    # Handle multiple classifications in a single row
    classifications = row['Illness Classification'].split(', ')
    for classification in classifications:
        classification_cleaned = classification.strip()  # Remove extra spaces
        
        # Add count for individual classifications
        if classification_cleaned in classification_counts:
            classification_counts[classification_cleaned] += 1
        else:
            classification_counts[classification_cleaned] = 1
        
        # Calculate means for each classification based on 'mAb' groups
        classification_key = '_'.join(classification_cleaned.split())
        mAb_group = row['mAb']  # Assuming the column name is 'mAb'
        
        if not pd.isnull(row['V1n']) and not pd.isnull(row['V2n']) and not pd.isnull(row['Qn']) and not pd.isnull(row['CLn']):
            if classification_key in mean_values:
                if mAb_group in mean_values[classification_key]:
                    mean_values[classification_key][mAb_group]['V1n'].append(row['V1n'])
                    mean_values[classification_key][mAb_group]['V2n'].append(row['V2n'])
                    mean_values[classification_key][mAb_group]['Qn'].append(row['Qn'])
                    mean_values[classification_key][mAb_group]['CLn'].append(row['CLn'])
                else:
                    mean_values[classification_key][mAb_group] = {
                        'V1n': [row['V1n']],
                        'V2n': [row['V2n']],
                        'Qn': [row['Qn']],
                        'CLn': [row['CLn']]
                    }
            else:
                mean_values[classification_key] = {
                    mAb_group: {
                        'V1n': [row['V1n']],
                        'V2n': [row['V2n']],
                        'Qn': [row['Qn']],
                        'CLn': [row['CLn']]
                    }
                }

# Calculate median of means for each classification and iqr bounds
iqr_bounds = {}
variables = ['V1n', 'V2n', 'Qn', 'CLn']
for classification, mAb_groups in mean_values.items():
    means = {
        'V1n': [],
        'V2n': [],
        'Qn': [],
        'CLn': []
    }
    for mAb_group, values in mAb_groups.items():
        if values['V1n'] and values['V2n'] and values['Qn'] and values['CLn']:
            means['V1n'].append(st.mean(values['V1n']))
            means['V2n'].append(st.mean(values['V2n']))
            means['Qn'].append(st.mean(values['Qn']))
            means['CLn'].append(st.mean(values['CLn']))
    if means['V1n'] and means['V2n'] and means['Qn'] and means['CLn']:
        median_values[classification + '_V1n'] = st.median(means['V1n'])
        median_values[classification + '_V2n'] = st.median(means['V2n'])
        median_values[classification + '_Qn'] = st.median(means['Qn'])
        median_values[classification + '_CLn'] = st.median(means['CLn'])

        #Calculating IQR bounds
        iqr_bounds[classification] = {}
        for variable in variables:  
            lower_bound = np.percentile(means[variable], 25)
            upper_bound = np.percentile(means[variable], 75)
            iqr_bounds[classification][variable] = (lower_bound, upper_bound)
        
# Convert the dictionaries to DataFrames with proper indexing
median_df = pd.DataFrame.from_dict(median_values, orient='index')


# Display counts for each unique value in 'Illness Classification'
print("Number of occurrences for each unique value in 'Illness Classification':")
print(classification_counts)

# Filter the DataFrame for each illness classification and calculate means
illness_classifications = [
    'Cardiovascular', 'Gastroenterology', 'Neurology', 'Oncology', 'Other', 'Rheumatology', 'Transplantation', 
    'Respiratory']

# Data for each illness classification

data_illness = {
    'Illness Classification': illness_classifications,
    'Median V1': [median_values.get(cls + '_V1n', None) for cls in illness_classifications],
    'Median V2': [median_values.get(cls + '_V2n', None) for cls in illness_classifications],
    'Median Q': [median_values.get(cls + '_Qn', None) for cls in illness_classifications],
    'Median CL': [median_values.get(cls + '_CLn', None) for cls in illness_classifications]
}

# Create DataFrame for illness classification median values
median_df_illness = pd.DataFrame(data_illness)

# Melt the DataFrame for Seaborn
melted_df = pd.melt(median_df_illness, id_vars='Illness Classification', var_name='Variable', value_name='Median Value')


# Extract counts for each illness classification
illness_counts = [f"{illness} (n={len(mean_values[illness])})" for illness in illness_classifications]


# # Calculate IQR bounds for each variable within each illness classification
illness_classes = ['Cardiovascular', 'Gastroenterology', 'Neurology', 'Oncology', 'Other',
                   'Rheumatology', 'Transplantation', 'Respiratory']


# Set the figure size
plt.figure(figsize=(12, 8))

# Plotting the medians for different variables using Seaborn
ax = sns.barplot(x='Illness Classification', y='Median Value',  hue='Variable', data=melted_df)

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting error bars for IQR bounds
for i, illness_class in enumerate(illness_classes):
    for j, variable in enumerate(variables):
        lower, upper = iqr_bounds[illness_class][variable]  # Retrieve bounds from the dictionary
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = median_df_illness[median_column_name].iloc[i]  # Retrieve the median value
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables) + 1) + (j + 1)) * bar_width - 0.5
        
        # Plot error bars
        ax.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)

# Annotate bars with their median values
for p in ax.patches:
    ax.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight='bold', color='black', xytext=(0, 5), 
                textcoords='offset points')


# Customize the x-axis labels with counts
plt.xlabel('Medical branch', fontsize = 14)
plt.ylabel('Median Values', fontsize = 14)

# Update x-axis ticks with number of models included
plt.xticks(range(len(illness_counts)), illness_counts, rotation=45, ha='right', fontsize=14)

# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Create the legend for both Seaborn plot and the dotted lines
handles_sns, labels_sns = ax.get_legend_handles_labels()
plt.legend(handles_sns, labels_sns, bbox_to_anchor=(1.02, 1), loc='upper left')

# Rotate x-axis ticks for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()
plt.show()


# signifiance testing


def kruskal_wallis_test(df, columns, group_column):
    results = {}
    for col in columns:
        groups = [df[df[group_column] == group][col] for group in df[group_column].unique()]
        statistic, p_value = kruskal(*groups)
        results[col] = {'Kruskal-Wallis Statistic': statistic, 'p-value': p_value}
    return results

columns_to_compare = ['V1n', 'V2n', 'Qn', 'CLn']
results = kruskal_wallis_test(IV_SC, columns_to_compare, 'Illness Classification')

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Kruskal-Wallis Statistic: {result['Kruskal-Wallis Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print("There is a significant difference in medians among medical branch types for this parameter.")
    else:
        print("There is no significant difference in medians among Imedical branch types for this parameter.")
    print("--------------")


# comparing pairwise groups

# def mann_whitney_u_test(df, columns, group_column, group_values):
#     results = {}
#     for col in columns:
#         group1 = df[df[group_column] == group_values[0]][col]
#         group2 = df[df[group_column] == group_values[1]][col]

#         statistic, p_value = mannwhitneyu(group1, group2)

#         results[col] = {'Mann-Whitney U Statistic': statistic, 'p-value': p_value}
#     return results

# columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  
# # ['Cardiovascular', 'Gastroenterology', 'Neurology', 'Oncology', 'Other', 'Rheumatology', 'Transplantation', 'Respiratory']
# group_values = ['Gastroenterology', 'Neurology']  # Groups to compare

# results_mannwhitney = mann_whitney_u_test(IV_SC, columns_to_test, 'mAb_Group', group_values)

# # Display results
# for col, result in results_mannwhitney.items():
#     print(f"Column: {col}")
#     print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
#     print(f"p-value: {result['p-value']}")
#     if result['p-value'] < 0.05:
#         print(f"There is a significant difference in '{col}' between 'Gastroenterology' and 'Neurology' groups.")
#     else:
#         print(f"There is no significant difference in '{col}' between 'Gastroenterology' and 'Neurology' groups.")
#     print("--------------")


#-------------------------------------- IgG type-----------------------------------------------------------

# Initialize dictionaries to store counts for each unique IG type and for mean/median values
ig_type_counts = {}
mean_values_ig = {}
median_values_ig = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    ig_type = row['IG type']
    if ig_type in ig_type_counts:
        ig_type_counts[ig_type] += 1
    else:
        ig_type_counts[ig_type] = 1
        
# Calculate medians using grouped means for each 'IgG' type

# For IgG1
IGG1_data = IV_SC[IV_SC['IG type'] == 'IgG1']
mean_V1n_IGG1 = IGG1_data.groupby("mAb")["V1n"].mean()
mean_V2n_IGG1 = IGG1_data.groupby("mAb")["V2n"].mean()
mean_Qn_IGG1 = IGG1_data.groupby("mAb")["Qn"].mean()
mean_CLn_IGG1 = IGG1_data.groupby("mAb")["CLn"].mean()

median_V1n_IGG1 = st.median(mean_V1n_IGG1.dropna())
median_V2n_IGG1 = st.median(mean_V2n_IGG1.dropna())
median_Qn_IGG1 = st.median(mean_Qn_IGG1.dropna())
median_CLn_IGG1 = st.median(mean_CLn_IGG1.dropna())

# For IgG2
IGG2_data = IV_SC[IV_SC['IG type'] == 'IgG2']
mean_V1n_IGG2 = IGG2_data.groupby("mAb")["V1n"].mean()
mean_V2n_IGG2 = IGG2_data.groupby("mAb")["V2n"].mean()
mean_Qn_IGG2 = IGG2_data.groupby("mAb")["Qn"].mean()
mean_CLn_IGG2 = IGG2_data.groupby("mAb")["CLn"].mean()

median_V1n_IGG2 = st.median(mean_V1n_IGG2.dropna())
median_V2n_IGG2 = st.median(mean_V2n_IGG2.dropna())
median_Qn_IGG2 = st.median(mean_Qn_IGG2.dropna())
median_CLn_IGG2 = st.median(mean_CLn_IGG2.dropna())

# For IgG4
IGG4_data = IV_SC[IV_SC['IG type'] == 'IgG4']
mean_V1n_IGG4 = IGG4_data.groupby("mAb")["V1n"].mean()
mean_V2n_IGG4 = IGG4_data.groupby("mAb")["V2n"].mean()
mean_Qn_IGG4 = IGG4_data.groupby("mAb")["Qn"].mean()
mean_CLn_IGG4 = IGG4_data.groupby("mAb")["CLn"].mean()

median_V1n_IGG4 = st.median(mean_V1n_IGG4.dropna())
median_V2n_IGG4 = st.median(mean_V2n_IGG4.dropna())
median_Qn_IGG4 = st.median(mean_Qn_IGG4.dropna())
median_CLn_IGG4 = st.median(mean_CLn_IGG4.dropna())

# Data for each IG type
data_ig = {
    'IG Type': ['IgG1', 'IgG2', 'IgG4'],
    'Median V1': [median_V1n_IGG1, median_V1n_IGG2, median_V1n_IGG4],
    'Median V2': [median_V2n_IGG1, median_V2n_IGG2, median_V2n_IGG4],
    'Median Q': [median_Qn_IGG1, median_Qn_IGG2, median_Qn_IGG4],
    'Median CL': [median_CLn_IGG1, median_CLn_IGG2, median_CLn_IGG4]
}

# Create a DataFrame from the data
df = pd.DataFrame(data_ig)

# Function to calculate lower and upper bounds of IQR for a variable within an IG type
def calculate_iqr_bounds(data, ig_type, variable):
    ig_data = data[data['IG type'] == ig_type][variable]
    lower_bound = np.percentile(ig_data.dropna(), 25)
    upper_bound = np.percentile(ig_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each IG type
ig_types = ['IgG1', 'IgG2', 'IgG4']
variables = ['V1n', 'V2n', 'Qn', 'CLn']

iqr_bounds = {}
for ig_type in ig_types:
    iqr_bounds[ig_type] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, ig_type, variable)
        iqr_bounds[ig_type][variable] = (lower, upper)



# Set the figure size
plt.figure(figsize=(12, 8))

# Plotting the medians for different variables using Seaborn
ax_ig = sns.barplot(x='IG Type', y='value', hue='variable', data=pd.melt(df, id_vars='IG Type', var_name='variable', value_name='value'))

# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting error bars for IQR bounds
for i, ig_type in enumerate(ig_types):
    for j, variable in enumerate(variables):
        lower, upper = iqr_bounds[ig_type][variable]  # Retrieve bounds from the dictionary
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df[df['IG Type'] == ig_type][median_column_name].values[0]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables)+1) + (j+1)) * bar_width - 0.5
        
        # Plot error bars
        ax_ig.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5 , linewidth=0.5)

# Annotate bars with their median values
for p in ax_ig.patches:
    ax_ig.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize = 14, fontweight='bold', color='black', xytext=(0, 5), 
                textcoords='offset points')


# Extract counts for each IG type
ig_counts = [f"{ig} (n={ig_type_counts[ig]})" for ig in ig_types]

# Customize the x-axis labels with IG types and counts
plt.xlabel('IgG class', fontsize = 14)
plt.ylabel('Median Values', fontsize = 14)

# Customize the x-axis ticks with IG types and counts
plt.xticks(range(len(ig_types)), ig_counts, rotation=0,fontsize = 14)

# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

# Show the plot
plt.tight_layout()
plt.show()


# signifiance testing

results = kruskal_wallis_test(IV_SC, columns_to_compare, 'IG type')

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Kruskal-Wallis Statistic: {result['Kruskal-Wallis Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print("There is a significant difference in medians among IG types for this parameter.")
    else:
        print("There is no significant difference in medians among IG types for this parameter.")
    print("--------------")


#----------------------Nomenclature-------------------------------------------------
def categorize_mAb(mAb):
    if mAb.endswith('zumab'):
        return 'Humanized'
    elif mAb.endswith('umab'):
        return 'Human'
    elif mAb.endswith('ximab'):
        return 'Chimeric'
    elif mAb.endswith('omab'):
        return 'Murine'
    else: # the single one remaining is a fully human mAb
        return 'Human'

# Apply the categorization function to create a new column 'mAb_Group'
IV_SC['mAb_Group'] = IV_SC['mAb'].apply(categorize_mAb)


# Initialize dictionaries to store counts for each unique 'mAb_Group' and for mean/median values
mab_group_counts = {}
mean_values_mab = {}
median_values_mab = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    mab_group = row['mAb_Group']
    if mab_group in mab_group_counts:
        mab_group_counts[mab_group] += 1
    else:
        mab_group_counts[mab_group] = 1
        
# Calculate medians using grouped means for each 'mAb_Group'

# For 'mAb_Group1'
mAb_Group1_data = IV_SC[IV_SC['mAb_Group'] == 'Human']
mean_V1n_mAb_Group1 = mAb_Group1_data.groupby("mAb")["V1n"].mean()
mean_V2n_mAb_Group1 = mAb_Group1_data.groupby("mAb")["V2n"].mean()
mean_Qn_mAb_Group1 = mAb_Group1_data.groupby("mAb")["Qn"].mean()
mean_CLn_mAb_Group1 = mAb_Group1_data.groupby("mAb")["CLn"].mean()

median_V1n_mAb_Group1 = st.median(mean_V1n_mAb_Group1.dropna())
median_V2n_mAb_Group1 = st.median(mean_V2n_mAb_Group1.dropna())
median_Qn_mAb_Group1 = st.median(mean_Qn_mAb_Group1.dropna())
median_CLn_mAb_Group1 = st.median(mean_CLn_mAb_Group1.dropna())

# For 'mAb_Group2'
mAb_Group2_data = IV_SC[IV_SC['mAb_Group'] == 'Humanized']
mean_V1n_mAb_Group2 = mAb_Group2_data.groupby("mAb")["V1n"].mean()
mean_V2n_mAb_Group2 = mAb_Group2_data.groupby("mAb")["V2n"].mean()
mean_Qn_mAb_Group2 = mAb_Group2_data.groupby("mAb")["Qn"].mean()
mean_CLn_mAb_Group2 = mAb_Group2_data.groupby("mAb")["CLn"].mean()

median_V1n_mAb_Group2 = st.median(mean_V1n_mAb_Group2.dropna())
median_V2n_mAb_Group2 = st.median(mean_V2n_mAb_Group2.dropna())
median_Qn_mAb_Group2 = st.median(mean_Qn_mAb_Group2.dropna())
median_CLn_mAb_Group2 = st.median(mean_CLn_mAb_Group2.dropna())

# # For 'mAb_Group3'
# mAb_Group3_data = IV_SC[IV_SC['mAb_Group'] == 'Other']
# mean_V1n_mAb_Group3 = mAb_Group3_data.groupby("mAb")["V1n"].mean()
# mean_V2n_mAb_Group3 = mAb_Group3_data.groupby("mAb")["V2n"].mean()
# mean_Qn_mAb_Group3 = mAb_Group3_data.groupby("mAb")["Qn"].mean()
# mean_CLn_mAb_Group3 = mAb_Group3_data.groupby("mAb")["CLn"].mean()

# median_V1n_mAb_Group3 = st.median(mean_V1n_mAb_Group3.dropna())
# median_V2n_mAb_Group3 = st.median(mean_V2n_mAb_Group3.dropna())
# median_Qn_mAb_Group3 = st.median(mean_Qn_mAb_Group3.dropna())
# median_CLn_mAb_Group3 = st.median(mean_CLn_mAb_Group3.dropna())

# For 'mAb_Group4'
mAb_Group4_data = IV_SC[IV_SC['mAb_Group'] == 'Chimeric']
mean_V1n_mAb_Group4 = mAb_Group4_data.groupby("mAb")["V1n"].mean()
mean_V2n_mAb_Group4 = mAb_Group4_data.groupby("mAb")["V2n"].mean()
mean_Qn_mAb_Group4 = mAb_Group4_data.groupby("mAb")["Qn"].mean()
mean_CLn_mAb_Group4 = mAb_Group4_data.groupby("mAb")["CLn"].mean()

median_V1n_mAb_Group4 = st.median(mean_V1n_mAb_Group4.dropna())
median_V2n_mAb_Group4 = st.median(mean_V2n_mAb_Group4.dropna())
median_Qn_mAb_Group4 = st.median(mean_Qn_mAb_Group4.dropna())
median_CLn_mAb_Group4 = st.median(mean_CLn_mAb_Group4.dropna())

# Data for each 'mAb_Group'
data_mab = {
    'mAb_Group': ['Human','Humanized', 'Chimeric'],
    'Median V1': [median_V1n_mAb_Group1, median_V1n_mAb_Group2, median_V1n_mAb_Group4],
    'Median V2': [median_V2n_mAb_Group1, median_V2n_mAb_Group2, median_V2n_mAb_Group4],
    'Median Q': [median_Qn_mAb_Group1, median_Qn_mAb_Group2, median_Qn_mAb_Group4],
    'Median CL': [median_CLn_mAb_Group1, median_CLn_mAb_Group2, median_CLn_mAb_Group4]
}

# Create a DataFrame from the data
df_mab = pd.DataFrame(data_mab)

# Function to calculate lower and upper bounds of IQR for a variable within an 'mAb_Group'
def calculate_iqr_bounds(data, mab_group, variable):
    mab_data = data[data['mAb_Group'] == mab_group][variable]
    lower_bound = np.percentile(mab_data.dropna(), 25)
    upper_bound = np.percentile(mab_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each 'mAb_Group'
mab_groups = ['Human', 'Humanized', 'Chimeric']
variables = ['V1n', 'V2n', 'Qn', 'CLn']

iqr_bounds_mab = {}
for mab_group in mab_groups:
    iqr_bounds_mab[mab_group] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, mab_group, variable)
        iqr_bounds_mab[mab_group][variable] = (lower, upper)

# Set the figure size
plt.figure(figsize=(12, 8))

# Plotting the medians for different variables using Seaborn
ax_mab = sns.barplot(x='mAb_Group', y='value', hue='variable', data=pd.melt(df_mab, id_vars='mAb_Group', var_name='variable', value_name='value'))

# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting error bars for IQR bounds
for i, mab_group in enumerate(mab_groups):
    for j, variable in enumerate(variables):
        lower, upper = iqr_bounds_mab[mab_group][variable]  # Retrieve bounds from the dictionary
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df_mab[df_mab['mAb_Group'] == mab_group][median_column_name].values[0]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables)+1) + (j+1)) * bar_width - 0.5
        
        # Plot error bars
        ax_mab.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)

# Annotate bars with their median values
for p in ax_mab.patches:
    ax_mab.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight= "bold", color='black', xytext=(0, 5), 
                textcoords='offset points')

# Extract counts for each 'mAb_Group'
mab_counts = [f"{mab} (n={mab_group_counts[mab]})" for mab in mab_groups]

# Customize the x-axis labels with 'mAb_Group' and counts
plt.xlabel('mAb_Group', fontsize=14)
plt.ylabel('Median Values', fontsize=14)

# Customize the x-axis ticks with 'mAb_Group' and counts
plt.xticks(range(len(mab_groups)), mab_counts, rotation=0, fontsize=14)

# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()

# Show the plot
plt.show()

# significance testing

results_mab = kruskal_wallis_test(IV_SC, columns_to_compare, 'mAb_Group')

# Display results
for col, result in results_mab.items():
    print(f"Column: {col}")
    print(f"Kruskal-Wallis Statistic: {result['Kruskal-Wallis Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print("There is a significant difference in medians among 'mAb_Groups' for this parameter.")
    else:
        print("There is no significant difference in medians among 'mAb_Groups' for this parameter.")
    print("--------------")

def mann_whitney_u_test(df, columns, group_column, group_values):
    results = {}
    for col in columns:
        group1 = df[df[group_column] == group_values[0]][col]
        group2 = df[df[group_column] == group_values[1]][col]

        statistic, p_value = mannwhitneyu(group1, group2)

        results[col] = {'Mann-Whitney U Statistic': statistic, 'p-value': p_value}
    return results

# comparing pairwise groups
columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  # Adjust as needed
group_values = ['Humanized', 'Human']  # Groups to compare

results_mannwhitney = mann_whitney_u_test(IV_SC, columns_to_test, 'mAb_Group', group_values)

# Display results
for col, result in results_mannwhitney.items():
    print(f"Column: {col}")
    print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print(f"There is a significant difference in '{col}' between 'Human' and 'Humanized' groups.")
    else:
        print(f"There is no significant difference in '{col}' between 'Human' and 'Humanized' groups.")
    print("--------------")

#------------------------Children--------------------------------------------------
# Initialize dictionaries to store counts for each unique column value and for mean/median values
pediatric_counts = {}
mean_values_child = {}
median_values_child = {}

# making extra column for binary values for pediatric model or not
IV_SC['pediatric model'] = IV_SC['Special patient group'].apply(lambda x: 1 if x == 'Children' else 0)
# or x == ">= 12 year" 

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    child = row['pediatric model']
    if child in pediatric_counts:
        pediatric_counts[child] += 1
    else:
        pediatric_counts[child] = 1

# changing keys of dictionary
pediatric_counts['No pediatric model'] = pediatric_counts.pop(0)
pediatric_counts['Pediatric model'] = pediatric_counts.pop(1)

# Filter the DataFrame for each Health type and calculate means or medians
pediatric_types = ['No pediatric model', 'Pediatric model']

# For Healthy patients included
pediatric_data = IV_SC[IV_SC['pediatric model'] == 1]
mean_V1n_child = pediatric_data.groupby("mAb")["V1n"].mean()
mean_V2n_child = pediatric_data.groupby("mAb")["V2n"].mean()
mean_Qn_child = pediatric_data.groupby("mAb")["Qn"].mean()
mean_CLn_child = pediatric_data.groupby("mAb")["CLn"].mean()

median_V1n_child = st.median(mean_V1n_child.dropna())
median_V2n_child = st.median(mean_V2n_child.dropna())
median_Qn_child = st.median(mean_Qn_child.dropna())
median_CLn_child = st.median(mean_CLn_child.dropna())

# For No healthy patients included
no_pediatric_data = IV_SC[IV_SC['pediatric model'] == 0]
mean_V1n_childno = no_pediatric_data.groupby("mAb")["V1n"].mean()
mean_V2n_childno = no_pediatric_data.groupby("mAb")["V2n"].mean()
mean_Qn_childno = no_pediatric_data.groupby("mAb")["Qn"].mean()
mean_CLn_childno = no_pediatric_data.groupby("mAb")["CLn"].mean()

median_V1n_childno = st.median(mean_V1n_childno.dropna())
median_V2n_childno = st.median(mean_V2n_childno.dropna())
median_Qn_childno = st.median(mean_Qn_childno.dropna())
median_CLn_childno = st.median(mean_CLn_childno.dropna())

# Data for each type
median_values_child = {
    'Patient population': ['No pediatric model', 'Pediatric model'],
    'Median V1': [median_V1n_childno, median_V1n_child],
    'Median V2': [median_V2n_childno, median_V2n_child],
    'Median Q': [median_Qn_childno, median_Qn_child],
    'Median CL': [median_CLn_childno, median_CLn_child]
}


# Create a DataFrame from the median values
df_childmod = pd.DataFrame(median_values_child)

# Function to calculate lower and upper bounds of IQR for a variable within a group
def calculate_iqr_bounds(data, group, variable):
    group_data = data[data['pediatric model'] == group][variable]
    lower_bound = np.percentile(group_data.dropna(), 25)
    upper_bound = np.percentile(group_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each group (Healthy and No healthy patients)
groups = [1, 0]  # Assuming '0' represents 'No pediatric model' and '1' represents 'pediatric model'
variables = ['V1n', 'V2n', 'Qn', 'CLn']

# Calculating IQR bounds for health data
iqr_bounds_pediatric = {}
for group in groups:
    iqr_bounds_pediatric[group] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, group, variable)
        iqr_bounds_pediatric[group][variable] = (lower, upper)


# Set the figure size
plt.figure(figsize=(12, 8))

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting the medians for different variables using Seaborn
ax_child = sns.barplot(x='Patient population', y='value', hue='variable', data=pd.melt(df_childmod, id_vars='Patient population', var_name='variable', value_name='value'))

# Annotate bars with their median values
for p in ax_child.patches:
    ax_child.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight ="bold", color='black', xytext=(0, 5), 
                textcoords='offset points')

# Extract the groups for labeling the x-axis ticks
child_groups = ['No pediatric model', 'Pediatric model']

# Plot error bars for each variable within each group
for i, group in enumerate(child_groups):
    for j, variable in enumerate(variables):
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df_childmod[df_childmod['Patient population'] == group][median_column_name].values[0]
        lower, upper = iqr_bounds_pediatric[i][variable]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables) + 1) + (j + 1)) * bar_width -0.5
        
        # Plot error bars
        ax_child.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)


# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Customize the x-axis labels and ticks
plt.xlabel('Model specification', fontsize=14)
plt.ylabel('Median Values', fontsize=14)

# Extract counts for each health type
pediatric_counts = [f"{i} (n={pediatric_counts[i]})" for i in pediatric_types]

# Customize the x-axis ticks with health types and counts
plt.xticks(range(len(pediatric_types)), pediatric_counts, rotation=0, fontsize=14)


# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()

# Show the plot
plt.show()

#significance testing

# Usage example
columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  # Adjust as needed
group_values = [1, 0]  # Values of 'Autoimmune Disorder' column for comparison
results = mann_whitney_u_test(IV_SC, columns_to_test, 'pediatric model', group_values)

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print(f"There is a significant difference in '{col}' between pediatric and non-pediatric models.")
    else:
        print(f"There is no significant difference in '{col}' between pediatric and non-pediatric models.")
    print("--------------")

    


#----------------lin/nonlin clearance difference---------------------------------
# Initialize dictionaries to store counts for each unique NonLinCL type and for mean/median values
nonlincl_counts = {}
mean_values_nonlincl = {}
median_values_nonlincl = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    nonlincl = row['NonLinCL']
    if nonlincl in nonlincl_counts:
        nonlincl_counts[nonlincl] += 1
    else:
        nonlincl_counts[nonlincl] = 1

# changing keys of dictionary
nonlincl_counts['Only linear clearance'] = nonlincl_counts.pop(0)
nonlincl_counts['Both linear and nonlinear clearance'] = nonlincl_counts.pop(1)

# Filter the DataFrame for each NonLinCL type and calculate means or medians
nonlincl_types = ['Only linear clearance', 'Both linear and nonlinear clearance']


# For nonlin + lin clearance models
nonlincl_data = IV_SC[IV_SC['NonLinCL'] == 1]
mean_V1n_nonlincl = nonlincl_data.groupby("mAb")["V1n"].mean()
mean_V2n_nonlincl = nonlincl_data.groupby("mAb")["V2n"].mean()
mean_Qn_nonlincl = nonlincl_data.groupby("mAb")["Qn"].mean()
mean_CLn_nonlincl = nonlincl_data.groupby("mAb")["CLn"].mean()

median_V1n_nonlincl = st.median(mean_V1n_nonlincl.dropna())
median_V2n_nonlincl = st.median(mean_V2n_nonlincl.dropna())
median_Qn_nonlincl = st.median(mean_Qn_nonlincl.dropna())
median_CLn_nonlincl = st.median(mean_CLn_nonlincl.dropna())

# For only linear clearance models
lincl_data = IV_SC[IV_SC['NonLinCL'] == 0]
mean_V1n_lincl = lincl_data.groupby("mAb")["V1n"].mean()
mean_V2n_lincl = lincl_data.groupby("mAb")["V2n"].mean()
mean_Qn_lincl = lincl_data.groupby("mAb")["Qn"].mean()
mean_CLn_lincl = lincl_data.groupby("mAb")["CLn"].mean()

median_V1n_lincl = st.median(mean_V1n_lincl.dropna())
median_V2n_lincl = st.median(mean_V2n_lincl.dropna())
median_Qn_lincl = st.median(mean_Qn_lincl.dropna())
median_CLn_lincl = st.median(mean_CLn_lincl.dropna())

# Data for each type
median_values_nonlincl = {
    'NonLinCL': ['Only Linear clearance', 'Both linear and nonlinear clearance'],
    'Median V1': [median_V1n_lincl, median_V1n_nonlincl],
    'Median V2': [median_V2n_lincl, median_V2n_nonlincl],
    'Median Q': [median_Qn_lincl, median_Qn_nonlincl],
    'Median CL': [median_CLn_lincl, median_CLn_nonlincl]
}

# Create a DataFrame from the median values
df_nonlincl = pd.DataFrame(median_values_nonlincl)

# Function to calculate lower and upper bounds of IQR for a variable within a group
def calculate_iqr_bounds(data, group, variable):
    group_data = data[data['NonLinCL'] == group][variable]
    lower_bound = np.percentile(group_data.dropna(), 25)
    upper_bound = np.percentile(group_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each group (Nonlinear and Other)
groups = [1, 0]  # Assuming '1' represents 'Both linear and nonlinear clearance' and '0' represents 'Only linear clearance'
variables = ['V1n', 'V2n', 'Qn', 'CLn']

# Calculating IQR bounds
iqr_bounds = {}
for group in groups:
    iqr_bounds[group] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, group, variable)
        iqr_bounds[group][variable] = (lower, upper)


# Set the figure size
plt.figure(figsize=(12, 8))


# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2


# Plotting the medians for different variables using Seaborn
ax_nonlincl = sns.barplot(x='NonLinCL', y='value', hue='variable', data=pd.melt(df_nonlincl, id_vars='NonLinCL', var_name='variable', value_name='value'))

# Annotate bars with their median values
for p in ax_nonlincl.patches:
    ax_nonlincl.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight="bold", color='black', xytext=(0, 5), 
                textcoords='offset points')

# Extract the groups for labeling the x-axis ticks
nonlincl_groups = ['Only Linear clearance', 'Both linear and nonlinear clearance']

# Plot error bars for each variable within each group
for i, group in enumerate(nonlincl_groups):
    for j, variable in enumerate(variables):
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df_nonlincl[df_nonlincl['NonLinCL'] == group][median_column_name].values[0]
        lower, upper = iqr_bounds[i][variable]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables)+1) + (j+1)) * bar_width - 0.5
        
        # Plot error bars
        ax_nonlincl.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)


# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Customize the x-axis labels and ticks
plt.xlabel('Clearance models', fontsize=14)
plt.ylabel('Median Values', fontsize=14)

# Extract counts for each IG type
nonlincl_counts = [f"{nonlin} (n={nonlincl_counts[nonlin]})" for nonlin in nonlincl_types]

# Customize the x-axis ticks with the NonLinCL groups
plt.xticks(range(len(nonlincl_types)), nonlincl_counts, rotation=0, fontsize=14)

# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()

# Show the plot
plt.show()

#significance testing

columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  # Adjust as needed
group_values = [1, 0]  # Values of 'Autoimmune Disorder' column for comparison
results = mann_whitney_u_test(IV_SC, columns_to_test, 'NonLinCL', group_values)

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print(f"There is a significant difference in '{col}' between Nonlinear clearance included/excluded groups.")
    else:
        print(f"There is no significant difference in '{col}' between Nonlinear clearance included/excluded groups.")
    print("--------------")






#---------------------Healthy patients------------------------------------------
# Initialize dictionaries to store counts for each unique column value and for mean/median values
health_counts = {}
mean_values_health = {}
median_values_health = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    health = row['Healthy patients included']
    if health in health_counts:
        health_counts[health] += 1
    else:
        health_counts[health] = 1

# changing keys of dictionary
health_counts['No healthy patients included'] = health_counts.pop(0)
health_counts['Healthy patients included'] = health_counts.pop(1)

# Filter the DataFrame for each Health type and calculate means or medians
health_types = ['No healthy patients included', 'Healthy patients included']

# For Healthy patients included
health_data = IV_SC[IV_SC['Healthy patients included'] == 1]
mean_V1n_health = health_data.groupby("mAb")["V1n"].mean()
mean_V2n_health = health_data.groupby("mAb")["V2n"].mean()
mean_Qn_health = health_data.groupby("mAb")["Qn"].mean()
mean_CLn_health = health_data.groupby("mAb")["CLn"].mean()

median_V1n_health = st.median(mean_V1n_health.dropna())
median_V2n_health = st.median(mean_V2n_health.dropna())
median_Qn_health = st.median(mean_Qn_health.dropna())
median_CLn_health = st.median(mean_CLn_health.dropna())

# For No healthy patients included
no_health_data = IV_SC[IV_SC['Healthy patients included'] == 0]
mean_V1n_healthno = no_health_data.groupby("mAb")["V1n"].mean()
mean_V2n_healthno = no_health_data.groupby("mAb")["V2n"].mean()
mean_Qn_healthno = no_health_data.groupby("mAb")["Qn"].mean()
mean_CLn_healthno = no_health_data.groupby("mAb")["CLn"].mean()

median_V1n_healthno = st.median(mean_V1n_healthno.dropna())
median_V2n_healthno = st.median(mean_V2n_healthno.dropna())
median_Qn_healthno = st.median(mean_Qn_healthno.dropna())
median_CLn_healthno = st.median(mean_CLn_healthno.dropna())

# Data for each type
median_values_health = {
    'Patient population': ['No healthy patients included', 'Healthy patients included'],
    'Median V1': [median_V1n_healthno, median_V1n_health],
    'Median V2': [median_V2n_healthno, median_V2n_health],
    'Median Q': [median_Qn_healthno, median_Qn_health],
    'Median CL': [median_CLn_healthno, median_CLn_health]
}


# Create a DataFrame from the median values
df_healthypop = pd.DataFrame(median_values_health)

# Function to calculate lower and upper bounds of IQR for a variable within a group
def calculate_iqr_bounds(data, group, variable):
    group_data = data[data['Healthy patients included'] == group][variable]
    lower_bound = np.percentile(group_data.dropna(), 25)
    upper_bound = np.percentile(group_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each group (Healthy and No healthy patients)
groups = [1, 0]  # Assuming '0' represents 'No healthy patients included' and '1' represents 'Healthy patients included'
variables = ['V1n', 'V2n', 'Qn', 'CLn']

# Calculating IQR bounds for health data
iqr_bounds_health = {}
for group in groups:
    iqr_bounds_health[group] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, group, variable)
        iqr_bounds_health[group][variable] = (lower, upper)


# Set the figure size
plt.figure(figsize=(12, 8))

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting the medians for different variables using Seaborn
ax_health = sns.barplot(x='Patient population', y='value', hue='variable', data=pd.melt(df_healthypop, id_vars='Patient population', var_name='variable', value_name='value'))

# Annotate bars with their median values
for p in ax_health.patches:
    ax_health.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight = "bold", color='black', xytext=(0, 5), 
                textcoords='offset points')

# Extract the groups for labeling the x-axis ticks
health_groups = ['No healthy patients included', 'Healthy patients included']

# Plot error bars for each variable within each group
for i, group in enumerate(health_groups):
    for j, variable in enumerate(variables):
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df_healthypop[df_healthypop['Patient population'] == group][median_column_name].values[0]
        lower, upper = iqr_bounds_health[i][variable]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables) + 1) + (j + 1)) * bar_width -0.5
        
        # Plot error bars
        ax_health.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)


# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Customize the x-axis labels and ticks
plt.xlabel('Patient population', fontsize=14)
plt.ylabel('Median Values', fontsize=14)

# Extract counts for each health type
health_counts = [f"{i} (n={health_counts[i]})" for i in health_types]

# Customize the x-axis ticks with health types and counts
plt.xticks(range(len(health_types)), health_counts, rotation=0, fontsize=14)


# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')


plt.tight_layout()

# Show the plot
plt.show()

#significance testing

# Usage example
columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  # Adjust as needed
group_values = [1, 0]  # Values of 'Autoimmune Disorder' column for comparison
results = mann_whitney_u_test(IV_SC, columns_to_test, 'Healthy patients included', group_values)

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print(f"There is a significant difference in '{col}' between healthy study population and diseased study population groups.")
    else:
        print(f"There is no significant difference in '{col}' between healthy study population and diseased study population groups.")
    print("--------------")

#------------------------Autoimmune disease group-------------------------------
# Initialize dictionaries to store counts for each unique column value and for mean/median values
autoimmune_counts = {}
mean_values_autoimmune = {}
median_values_autoimmune = {}

# Iterate over each row in the DataFrame
for index, row in IV_SC.iterrows():
    autoimmune = row['Autoimmune Disorder']
    if autoimmune in autoimmune_counts:
        autoimmune_counts[autoimmune] += 1
    else:
        autoimmune_counts[autoimmune] = 1

# changing keys of dictionary
autoimmune_counts['No autoimmune diseased patients'] = autoimmune_counts.pop(0)
autoimmune_counts['Autoimmune diseased patients'] = autoimmune_counts.pop(1)

# Filter the DataFrame for each Autoimmune type and calculate means or medians
autoimmune_types = ['No autoimmune diseased patients', 'Autoimmune diseased patients']

# For Autoimmune diseased patients
autoimmune_data = IV_SC[IV_SC['Autoimmune Disorder'] == 1]
mean_V1n_autoimmune = autoimmune_data.groupby("mAb")["V1n"].mean()
mean_V2n_autoimmune = autoimmune_data.groupby("mAb")["V2n"].mean()
mean_Qn_autoimmune = autoimmune_data.groupby("mAb")["Qn"].mean()
mean_CLn_autoimmune = autoimmune_data.groupby("mAb")["CLn"].mean()

median_V1n_autoimmune = st.median(mean_V1n_autoimmune.dropna())
median_V2n_autoimmune = st.median(mean_V2n_autoimmune.dropna())
median_Qn_autoimmune = st.median(mean_Qn_autoimmune.dropna())
median_CLn_autoimmune = st.median(mean_CLn_autoimmune.dropna())

# For No autoimmune diseases patients
no_autoimmune_data = IV_SC[IV_SC['Autoimmune Disorder'] == 0]
mean_V1n_no_autoimmune = no_autoimmune_data.groupby("mAb")["V1n"].mean()
mean_V2n_no_autoimmune = no_autoimmune_data.groupby("mAb")["V2n"].mean()
mean_Qn_no_autoimmune = no_autoimmune_data.groupby("mAb")["Qn"].mean()
mean_CLn_no_autoimmune = no_autoimmune_data.groupby("mAb")["CLn"].mean()

median_V1n_no_autoimmune = st.median(mean_V1n_no_autoimmune.dropna())
median_V2n_no_autoimmune = st.median(mean_V2n_no_autoimmune.dropna())
median_Qn_no_autoimmune = st.median(mean_Qn_no_autoimmune.dropna())
median_CLn_no_autoimmune = st.median(mean_CLn_no_autoimmune.dropna())

# Data for each type
median_values_autoimmune = {
    'Patient population': ['No autoimmune diseased patients', 'Autoimmune diseased patients'],
    'Median V1': [median_V1n_no_autoimmune, median_V1n_autoimmune],
    'Median V2': [median_V2n_no_autoimmune, median_V2n_autoimmune],
    'Median Q': [median_Qn_no_autoimmune, median_Qn_autoimmune],
    'Median CL': [median_CLn_no_autoimmune, median_CLn_autoimmune]
}

# Create a DataFrame from the median values
df_autoimmune = pd.DataFrame(median_values_autoimmune)

# Function to calculate lower and upper bounds of IQR for a variable within a group
def calculate_iqr_bounds(data, group, variable):
    group_data = data[data['Autoimmune Disorder'] == group][variable]
    lower_bound = np.percentile(group_data.dropna(), 25)
    upper_bound = np.percentile(group_data.dropna(), 75)
    return lower_bound, upper_bound

# Calculate IQR bounds for each variable within each group (Autoimmune and No autoimmune patients)
groups = [0, 1]  # Assuming '0' represents 'No autoimmune diseases patients' and '1' represents 'Autoimmune diseased patients'
variables = ['V1n', 'V2n', 'Qn', 'CLn']

# Calculating IQR bounds for autoimmune data
iqr_bounds_autoimmune = {}
for group in groups:
    iqr_bounds_autoimmune[group] = {}
    for variable in variables:
        lower, upper = calculate_iqr_bounds(IV_SC, group, variable)
        iqr_bounds_autoimmune[group][variable] = (lower, upper)

# Set the figure size
plt.figure(figsize=(12, 8))

# Calculate the number of variables and set the bar width
num_variables = len(variables)
bar_width = 0.2

# Plotting the medians for different variables using Seaborn
ax_autoimmune = sns.barplot(x='Patient population', y='value', hue='variable', data=pd.melt(df_autoimmune, id_vars='Patient population', var_name='variable', value_name='value'))

# Annotate bars with their median values
for p in ax_autoimmune.patches:
    ax_autoimmune.annotate(f'{p.get_height():.2f}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='center', fontsize=12, fontweight='bold', color='black', xytext=(0, 5), 
                textcoords='offset points')

# Extract the groups for labeling the x-axis ticks
autoimmune_groups = ['No autoimmune diseased patients', 'Autoimmune diseased patients']

# Plot error bars for each variable within each group
for i, group in enumerate(autoimmune_groups):
    for j, variable in enumerate(variables):
        median_column_name = f"Median {variable[:-1]}"  # Construct column name without the last letter 'n'
        median = df_autoimmune[df_autoimmune['Patient population'] == group][median_column_name].values[0]
        lower, upper = iqr_bounds_autoimmune[i][variable]
        
        # Calculate x-coordinate for the error bars with an offset
        x_offset = (i * (len(variables) + 1) + (j + 1)) * bar_width - 0.5
        
        # Plot error bars
        ax_autoimmune.errorbar(x_offset, median, yerr=[[median - lower], [upper - median]], fmt='none', capsize=5, color='darkblue', zorder=5, linewidth=0.5)

# Plotting overall median values with dotted lines for each variable (thinner lines)
plt.axhline(y=median_V1n, color='blue', linestyle='--', label='Overall Median V1', linewidth=1.2)
plt.axhline(y=median_V2n, color='orange', linestyle='--', label='Overall Median V2', linewidth=1.2)
plt.axhline(y=median_Qn, color='green', linestyle='--', label='Overall Median Q', linewidth=1.2)
plt.axhline(y=median_CLn, color='red', linestyle='--', label='Overall Median CL', linewidth=1.2)

# Customize the x-axis labels and ticks
plt.xlabel('Patient population', fontsize=14)
plt.ylabel('Median Values', fontsize=14)

# Extract counts for each autoimmune type
autoimmune_counts = [f"{i} (n={autoimmune_counts[i]})" for i in autoimmune_types]

# Customize the x-axis ticks with autoimmune types and counts
plt.xticks(range(len(autoimmune_types)), autoimmune_counts, rotation=0, fontsize=14)

# Create a legend for overall median values, place it outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()


# Show the plot
plt.show()


# Significance testing

# Usage example
columns_to_test = ['V1n', 'V2n', 'Qn', 'CLn']  # Adjust as needed
group_values = [1, 0]  # Values of 'Autoimmune Disorder' column for comparison
results = mann_whitney_u_test(IV_SC, columns_to_test, 'Autoimmune Disorder', group_values)

# Display results
for col, result in results.items():
    print(f"Column: {col}")
    print(f"Mann-Whitney U Statistic: {result['Mann-Whitney U Statistic']}")
    print(f"p-value: {result['p-value']}")
    if result['p-value'] < 0.05:
        print(f"There is a significant difference in '{col}' between Autoimmune disordered or not groups.")
    else:
        print(f"There is no significant difference in '{col}' between Autoimmune disordered or not groups.")
    print("--------------")



#--------------------------Hierarchical clustering-------------------------------------


# Select the columns for clustering
columns_for_clustering = ['V1n', 'V2n', 'Qn', 'CLn']

# Extract the data for clustering
data_for_clustering = IV_SC[columns_for_clustering]

# Calculate the linkage matrix using Ward's method
Z = linkage(data_for_clustering, method='ward')

IV_SC['Autoimmune Disorder'] = IV_SC['Autoimmune Disorder'].astype(int)
IV_SC['Healthy patients included'] = IV_SC['Healthy patients included'].astype(int)
IV_SC['NonLinCL'] = IV_SC['NonLinCL'].astype(int)

# Define custom labels
custom_labels = IV_SC[['IG type', 'Illness Classification', 'Autoimmune Disorder', 'Healthy patients included', 'NonLinCL', 'mAb_Group', 'pediatric model']]
custom_labels = custom_labels.apply(lambda row: ' '.join(row.values.astype(str)), axis=1).values

# Plotting the dendrogram with custom labels and same color for leaves in the same cluster
plt.figure(figsize=(12, 8))
dendrogram(Z, labels=custom_labels, leaf_rotation=90, color_threshold=8, above_threshold_color='gray')
plt.xlabel('Cluster labels for the columns IgG type, Medical branch, Autoimmune Disorder, Healthy patients included, Nonlinear CL included, Nomenclature & Pediatric model respectively', fontsize=9)
plt.ylabel('Distance')
plt.tight_layout()
plt.show()


# Perform clustering for two clusters
n_clusters = 3
cluster_labels = fcluster(Z, n_clusters, criterion='maxclust')

# Store cluster assignments in the dataframe
IV_SC[f'Cluster_{n_clusters}'] = cluster_labels

# Create a list of categorical variables for analysis
categorical_variables = ['IG type', 'Illness Classification', 'Autoimmune Disorder', 'Healthy patients included', 'NonLinCL', 'mAb_Group', 'pediatric model']

# Store results in a dictionary
results = {}

# Compute the ratio of the presence of each class in each cluster (ignoring "other")
for variable in categorical_variables:
    # Count the occurrences of each class in each cluster
    cluster_counts = IV_SC[IV_SC[variable] != 'ther'].groupby([f'Cluster_{n_clusters}', variable]).size().unstack(fill_value=0)
    
    # Calculate the ratio of classes in each cluster
    cluster_ratios = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)
    
    # Store the results
    results[variable] = cluster_ratios

# Concatenate results into a single DataFrame
result_table = pd.concat(results.values(), keys=results.keys(), axis=1)

# Display the table
print(result_table)

#----------------------------IIV determination-------------------------------------------
df = pd.read_excel("Theo_mAbs_PK(1).xlsx", sheet_name = "Other")


selected_columns= ['Pdf', 'mAb', 'Route', 'ka', 'V1', 'CL', 'Q',
       'V2', 'V_total', 'F', 'V1_BSV %', 'V2_BSV %', 'Q_BSV %', 'CL_BSV %',
       'VM_BSV %', 'KM_BSV %', 'Ka_BSV %', 'F_BSV %', 'Prop_res %',
       'CL/F_BSV %', 'V1/F_BSV %', 'Q/F_BSV %', 'V2/F_BSV %', 'IG type',
       'Patients', 'Special patient group', 'Illness Classification',
       'Autoimmune Disorder', 'Healthy patients included','Model',
       'NonLinCL', 'V1_WT', 'CL_WT', 'Q_WT', 'V2_WT', 'WT', 'N', 'V1_Sx',
       'CL_Sx', 'V2_Sx', 'SEX', 'V1_Alb', 'CL_Alb', 'V2_Alb', 'ALB', 'CL_ADA',
       'Cov_other', 'V1n', 'CLn', 'Qn', 'V2n', 'Mcount']

data = df[selected_columns]

# selecting only 2 compartment models
data = data[data["Model"] == 2]
# selecting only lin clearance and nonlinear clearance models
data = data[data['NonLinCL'] == 1]
# Excluding models for children
data = data[data["Special patient group"]!= "Children"]

# Only including IV models
desired_routes = ["IV", "IV/SC", "SC/IV"]
IV_SC= data[data["Route"].isin(desired_routes)].reset_index(drop=True)

# Removing models with Nan values for V1n, V2n, Qn and Cln
IV_SC = IV_SC[IV_SC['V2n'].notna()]

# Correcting the Sex covariates for reference Female, which results in Male dependent covariates:
for index, value in enumerate(IV_SC["SEX"]):
    if value == "M":
        theta_old_V1Sx = IV_SC["V1_Sx"].iloc[index]
        theta_old_CLSx = IV_SC["CL_Sx"].iloc[index]
        IV_SC["V1_Sx"].iloc[index] = (1/(1+theta_old_V1Sx)) - 1
        IV_SC["CL_Sx"].iloc[index] = (1/(1+theta_old_CLSx)) - 1


# ------------------------IIV values---------------------------------------------------
# Group by "mAb" column and calculate mean for each group
mean_V1_IIV_grouped = IV_SC.groupby("mAb")['V1_BSV %'].mean()
mean_V2_IIV_grouped = IV_SC.groupby("mAb")['V2_BSV %'].mean()
mean_Q_IIV_grouped = IV_SC.groupby("mAb")['Q_BSV %'].mean()
mean_CL_IIV_grouped = IV_SC.groupby("mAb")['CL_BSV %'].mean()
mean_Vm_IIV_grouped = IV_SC.groupby("mAb")['VM_BSV %'].mean()
mean_Km_IIV_grouped = IV_SC.groupby("mAb")['KM_BSV %'].mean()

# Take the overall median of the grouped means
median_V1_IIV = st.median(mean_V1_IIV_grouped.dropna())
median_V2_IIV = st.median(mean_V2_IIV_grouped.dropna())
median_Q_IIV = st.median(mean_Q_IIV_grouped.dropna())
median_CL_IIV = st.median(mean_CL_IIV_grouped.dropna())
median_Vm_IIV = st.median(mean_Vm_IIV_grouped.dropna())
median_Km_IIV = st.median(mean_Km_IIV_grouped.dropna())


# Plotting distributions of V1 and V2 and Q and CL
plt.figure(figsize=(10, 6))
IV_SC['V1_BSV %'].plot.density(color='blue', label='V1 IIV')
IV_SC['V2_BSV %'].plot.density(color='orange', label='V2 IIV')
IV_SC['Q_BSV %'].plot.density(color='green', label='Q IIV')
IV_SC['CL_BSV %'].plot.density(color='red', label='CL IIV')
IV_SC['VM_BSV %'].plot.density(color='purple', label='VM IIV')
IV_SC['KM_BSV %'].plot.density(color='yellow', label='KM IIV')

plt.legend()
plt.show()


counts_columns = IV_SC.count()
median_values_IIV = pd.DataFrame({
    'Parameter': ['V1 IIV', 'V2 IIV', 'Q IIV', 'CL IIV', 'VM IIV', 'KM IIV'],
    'Median': [median_V1_IIV,median_V2_IIV, median_Q_IIV,median_CL_IIV, median_Vm_IIV, 
               median_Km_IIV],
    'Counts': [counts_columns.loc['V1_BSV %'], counts_columns.loc['V2_BSV %'], counts_columns.loc['Q_BSV %'],
               counts_columns.loc['CL_BSV %'], counts_columns.loc['VM_BSV %'], counts_columns.loc['KM_BSV %'],]
})

print(median_values_IIV)

median_values_IIV.to_csv('IIV.csv')



