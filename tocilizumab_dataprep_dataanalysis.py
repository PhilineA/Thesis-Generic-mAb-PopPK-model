# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:32:54 2023

@author: adolf01p
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import statistics as st
import numpy as np
import statistics as st

os.chdir("C:\\Users\\adolf01p\\OneDrive - Sanquin\\Documenten\\mAbs datasets")

data = pd.read_excel("nm_002.xlsx")
data_ref = pd.read_csv("nm_ntz_ms.csv")

columns_data= ['C', 'ID', 'TIME', 'TAD', 'DV', 'AMT', 'II', 'SS', 'MDV', 'EVID',
       'STUDYGROUP', 'SEX', 'WT', 'AGE', 'HT', 'LBM', 'BMI', 'DD', 'STARTTCZ',
       'Switcher', 'IgM', 'ACPA', 'ERO', 'Morning stiffness', 'VASpain',
       'VASgeneral', 'VASphysician', 'TJC68', 'TJC28', 'SJC66', 'SJC28',
       'DAS28_ESR', 'SDAI', 'CDAI', 'RAPID_total', 'BSE', 'BSEhl', 'CRP', 'Hb',
       'MCV', 'Leukocytes', 'Thrombocytes', 'Creatinine', 'ALT', 'MTX',
       'MTXyn']

columns_data_ref = ['NMID', 'TIME', 'TAD', 'DV', 'MDV', 'EVID', 'AMT', 'SS', 'II', 'RATE',
       'WT', 'LENGTH', 'SEX', 'AGE', 'TSFD', 'SNR', 'SMPL']

columns_del = ['C', 'LBM', 'BMI', 'DD', 'STARTTCZ',
       'Switcher', 'IgM', 'ACPA', 'ERO', 'Morning stiffness', 'VASpain',
       'VASgeneral', 'VASphysician', 'TJC68', 'TJC28', 'SJC66', 'SJC28',
       'DAS28_ESR', 'SDAI', 'CDAI', 'RAPID_total', 'BSE', 'BSEhl', 'CRP', 'Hb',
       'MCV', 'Leukocytes', 'Thrombocytes', 'Creatinine', 'ALT', 'MTX',
       'MTXyn']

data = data.rename({"ID": "NMID", "STUDYGROUP" : "SNR", "HT" : "LENGTH"})

# Deleting unknown columns for model
for i in columns_del:
    del data[i]

# transforming TIME & TAD & II from hours to days
data["TIME"] = data["TIME"]/24

for count, value in enumerate(data["TAD"]):
    if value != ".":
        data["TAD"].iloc[count] = data["TAD"].iloc[count]/24

for count, value in enumerate(data["II"]):
    if value != ".":
        data["II"].iloc[count] = data["II"].iloc[count]/24

# exporting data
data.to_csv("nm_002.csv", index= False)



#-------------------determining data characteristics-------------------------------
# Keep only the first occurrence of each unique ID
unique_data = data.drop_duplicates(subset='ID', keep='first')

patients = len(unique_data)
print(f"Number of patients: {patients}") #12

#Number of males
males = unique_data["SEX"].eq(1).sum()
print(f"Number of males out of {patients} patients: {males}") #2

#median wt
print(f"median WT: {st.median(unique_data['WT'])}") #75.5
#IQR wt
lower_IQR_wt = np.percentile(unique_data['WT'], 25)
upper_IQR_wt = np.percentile(unique_data['WT'], 75)
print(f"IQR WT: {lower_IQR_wt}-{upper_IQR_wt}") #65.75 - 82.25

#median HT
unique_data['HT'] = pd.to_numeric(data['HT'], errors='coerce')
print(f"median HT: {st.median(unique_data['HT'])}") # 167.0
#IQR wt
lower_IQR_ht = np.percentile(unique_data['HT'].dropna(), 25)
upper_IQR_ht = np.percentile(unique_data['HT'].dropna(), 75)
print(f"IQR HT: {lower_IQR_ht}-{upper_IQR_ht}") # 157.0-176.0

#median age
print(f"median AGE: {st.median(unique_data['AGE'])}") #63
#IQR wt
lower_IQR_age = np.percentile(unique_data['AGE'].dropna(), 25)
upper_IQR_age = np.percentile(unique_data['AGE'].dropna(), 75)
print(f"IQR AGE: {lower_IQR_age}-{upper_IQR_age}")  # 56.25-71.5



