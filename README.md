This repository contains the (non-confidental) datasets and python scripts used for my Master thesis titled: "Meta-Analysis Insights into Monoclonal Antibody Population Pharmacokinetics:
Development of a Generic Model"

The supplementary material is stored in "Supplementary.docx"

The collected PopPK models for the meta-analysis and their characteristics are stored in "Theo_mAbs_PK(1).xlsx"

The meta-analysis and clustering on the colllected models in "Theo_mAbs_PK(1).xlsx" is perfromed using "Meta_analysis_and_clustering.py"

"dOFV_plots.py" is the script used to plot the dOFV in the corresponding dOFV excel sheets.

"MAPE_plots.py" plots the MAPE from the corresponding MAPE Excel files

The output tables from the NONMEM model runs using either Pirana or Python are extracted in the "Extracting_output_tables.py" script

The "Residual_plots.py" plots several residual measures using the model outcomes and extracted output tables.

"tocilizumab_dataprep_dataanalysis.py" is used to pre-proces the dataset for tocilizumab to make it suitable for the NONMEM models.
