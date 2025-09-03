README for Predicting-Chemicals Project

This project accompanies the MSc Data Science dissertation:
"From Structures to Risk: Predicting p53-Mediated Carcinogenicity in Atmospheric Chemicals Using Machine Learning."

Author: Jashwanthikaa Duraikannu Rajendiran
Course: MSc Data Science
Institution: University of York

--------------------------------------------------
Project Overview
--------------------------------------------------
This project develops and validates machine learning models (XGBoost) to predict the 
p53-mediated carcinogenicity of atmospheric chemicals. It integrates:

*EPA CompTox ATG_p53_CIS bioassay data (toxicology labels).
*Cheminformatics fingerprints (MACCS, PubChem, Klekota-Roth, SMARTS, ring system).
*MassBank (LCSB contributor) MS/MS spectra for external validation.
*Real limonene SOA chamber data to demonstrate environmental application.

--------------------------------------------------
Project Structure
--------------------------------------------------
The project is organised into the following main directories:

**Code/** -> All R scripts for preprocessing, fingerprint generation, model building, validation, and plotting.
**Data/** -> Raw, intermediate, and processed datasets (CompTox assay, MassBank, CAS-matched sets).
**Results/** -> Model outputs, trained models, confusion matrices, and prediction files.
**Figs/** -> Visual outputs (ROC curves, SHAP plots, feature importance, threshold optimisation).
**Real Sample/** -> Validation on limonene SOA chamber data (spectra, SIRIUS predictions, plots).
**validation compound/** -> Compound-level external validation folders, one per compound.
**Docs/** -> Supporting documentation and dissertation text.
**Predicting-Chemicals.Rproj** -> RStudio project file for reproducibility.
**FinalReport.docx** -> Dissertation document draft.
**README.txt** -> This file.

--------------------------------------------------
Code (Scripts)
--------------------------------------------------
The Code/ folder contains R scripts used for different stages of the pipeline:

**Preprocessing.R/** -> Initial cleaning of raw assay data, filtering compounds, removing missing values.

**Filter compounds.R/** -> Applies CHNOS-only filtering and prepares final subset of compounds.

**chemical compound separation.R/** -> Separates compounds based on structural and elemental criteria.

**CAS Finder.R/** -> Extracts CAS numbers from datasets and matches them with CompTox/MassBank.

**Compounds_Details.R/** -> Generates metadata for compounds, including molecular descriptors.

**RCDK Fingerprint Calculation.R/** -> Generates multiple molecular fingerprints (MACCS, PubChem, Klekota-Roth, SMARTS, ring system).

**Removing Correlated Fingerprints.R / Removing Correlated Fingerprints1.R/** -> Identifies and removes highly correlated fingerprint features to reduce redundancy.

**model building 2.R / model building 3.R / model_building(1-3).R / yes_no model building.R/** -> Multiple iterations of model development using XGBoost with different data splits and thresholds.

**best tuning in 100.R / Thershold(0.65)100.R / Thershold(0.65)200.R / Thershold(0.65)300.R/** -> Hyperparameter tuning experiments with different rounds and probability thresholds.

**Final_Model.R / Carcinogenicity Model 2 (1).R/** -> Final optimised XGBoost models trained on the balanced dataset.

**finalmodelwithpaneldistr.R /** -> Panel distribution plots for visualising active/inactive separation.

**model_comparisons.R/** -> Comparison of multiple models (100, 200, 300 rounds; different thresholds).

**Shap Analysis.R/** -> Computes SHAP values for model interpretability.

**XGBoostfinalplots.R / finalplots.R / final_Modelwsplots.R / figexplain.R/** -> Generates final visualisations including ROC curves, calibration plots, and feature importance.

**Heatmap.R/** -> Heatmap of feature correlations or importance.

**counts.R/** -> Counts active/inactive distributions in datasets.

**Finding Validation(set).R / LCSB code.R / Multiclass MB code.R/** -> Validation code linking model predictions with external LC-MS data (MassBank, SIRIUS).

**SIRIUS_Predicted_Fingerprints.R/** ?-> Integrates fingerprint predictions obtained from SIRIUS spectral analysis.

**RealSampleLimoneneNegative.R / Limonenebasedplots.R/** -> Real-sample testing on limonene SOA chamber data and related visualisations.

**Validating Actual and Predicted.R/** -> Compares predicted vs. actual activity for external validation compounds.

**Chemical analysis.R / data bank Finder.R/** -> Utility scripts for chemical composition analysis and data matching across databases.

--------------------------------------------------
Data
--------------------------------------------------
This folder contains the datasets used for model training, preprocessing, and validation.

**Key Files

*This project uses publicly available data from the **EPA CompTox Chemicals Dashboard**:  
-> Main dashboard: [https://comptox.epa.gov/dashboard/](https://comptox.epa.gov/dashboard/)  
-> Assay endpoint: [ATG_p53_CIS](https://comptox.epa.gov/dashboard/assay-endpoints/ATG_p53_CIS)  
  (p53 activation assay, used for labelling compounds as Active/Inactive).

The initial download included all assay compounds, but a subset was removed during preprocessing 
(e.g., non-CHNOS compounds, missing values, and highly correlated features).  

*FinalData.xlsx  
  -> This is the **main processed dataset** used for model training and evaluation.  
  -> Created after all preprocessing steps, including:  
     * Removal of missing values  
     * Filtering of CHNOS-only compounds  
     * Calculation of cheminformatics fingerprints  
     * Removal of highly correlated and near-zero variance features  
  ->This dataset represents the final balanced version used for XGBoost model building.  

*FinalData(1).xlsx  
  -> This is the same as **FinalData**, but with **column names adjusted according to SIRIUS output**.  
  -> Used for external validation, where SIRIUS-generated fingerprints and identifiers were matched against the processed dataset.  
  -> Ensures compatibility between model input and SIRIUS validation features.  

**Supporting Files

* Dataset.xlsx / DataF.xlsx / D1.xlsx -> Intermediate datasets generated during filtering and feature preparation.  
* classification.csv -> Full ATG_p53_CIS assay activity classification (raw binary labels).  
* CCD-Batch-Search_*.csv -> CAS number batch search results used to cross-match compounds with MassBank.  

*All_319/ -> Subset of 319 compounds cross-matched with MassBank (CAS-based).  
* CC1/ -> Contains SIRIUS-extracted candidate structures for validation.  
* **MassBank-data-2025.05.1.zip** ??? Full MassBank archive (public release).  
  -> Source: [MassBank-data GitHub release (2025.05.1)](https://github.com/MassBank/MassBank-data/releases/tag/2025.05.1)  
  -> This archive contains data from multiple contributors.  
  -> For this project, only the **Luxembourg Centre for Systems Biomedicine (LCSB)** contributor subset was extracted and used for external validation.

--------------------------------------------------
Results
--------------------------------------------------

The results/ directory contains model outputs, prediction tables, and saved models:

** Model training outputs

Best Tuning on training data.xlsx -> Hyperparameter tuning summary.

model building 3 results.xlsx, model training 2 results.xlsx -> Intermediate model runs.

Final training (train/test).xlsx -> Performance results for final split training/testing.

-> Validation outputs

Validation results.xlsx -> External validation results on independent data.

Sirius_predictions.xlsx / Sirius_predictions_names.xlsx -> Predictions on SIRIUS-processed spectral data.

Prediction_Comparison.xlsx / Prediction_Comparison_Unique.xlsx -> Comparison across different models.

Formula_Summary_Unique.xlsx -> Summary of predicted active/inactive chemical formulas.

-> Final trained XGBoost models (.rds)

Saved trained models: xgb_final_model_trained_on_all_data.rds, xgb_final_model_trained_on_all_datas.rds,
xgb_final_model_trained_on_all_dataF1.rds, xgb_final_model_trained_on_all_dataFinal.rds,
xgb_final_model_trained_on_all_datatrail.rds.

These allow reproducibility without retraining.

-> Prediction tables

Multiple files (xgb_model_predictions_lessdata*, xgb_predictions_trained_on_all_data*) -> Model predictions under different datasets and thresholds.

xgboost_predictions_with_metadata.xlsx -> Predictions merged with compound metadata.

xgb_multiclass_predictions.xlsx -> Exploratory multiclass run.

-> Evaluation metrics

confusion_matrix_final_model1.xlsx -> Confusion matrix of final tuned model.

--------------------------------------------------
Figures
--------------------------------------------------
The figs/ directory contains all plots and visualisations generated during the project. These figures were used for model evaluation, interpretation, and final reporting.

** Subfolders

*figs_final/ -> Final versions of plots selected for inclusion in the dissertation.

*figs_with_SMILES/ -> Plots showing predicted activities alongside molecular structures (SMILES visualisations).

*Plots(Models)/ -> Figures from model training and testing runs (different nrounds, thresholds).

*Plots(Modelsb)/ -> Alternative model comparison plots from exploratory builds.

*plots_crisp/ -> Cleaned and simplified plots for presentation.

*SHAP/ ->  SHAP summary plots showing feature contributions across all compounds.

*SHAP Analysis/ -> Detailed SHAP breakdown plots (per-feature contribution visualisations).

-> Standalone figures

*Actual vs Predicted.png / Actual vs Predicted (file) -> Shows comparison between predicted and observed activity classes.

*ROC_curve_final_model.png -> ROC curve for the best-performing final tuned model.

*roc_xgboost1.pdf -> Summary ROC curves from XGBoost model runs.

*thershold optimization plot.png -> Threshold optimisation curve comparing trade-offs in sensitivity/specificity.

*intro img1.png -> Supporting figure used in the Introduction section of the dissertation.

-> These figures illustrate:

*Model discrimination ability (ROC curves).

*The effect of threshold adjustments (sensitivity vs. specificity trade-off).

*Comparison of multiple model builds (100, 200, 300 rounds).

*Model interpretability (SHAP feature contribution plots).

*Validation results (predicted vs. actual activity).

*Structural comparisons (SMILES-annotated plots).

--------------------------------------------------
Real Sample Validation
--------------------------------------------------
The real_sample/ directory contains the data, scripts, and results from the application of the final trained model to limonene secondary organic aerosol (SOA) chamber experiments.

*LimoneneNegative.mgf -> Raw MS/MS spectra from the limonene SOA experiment.

*LimoneneNegative/ -> Folder containing curated negative mode MS/MS spectra.

*Structures/ -> Structural candidates identified via SIRIUS.

*Literature_vs_SIRIUS_summary.csv -> Comparison table between SIRIUS-predicted formulas and literature reports (e.g., Witkowski & Gierczak, 2017).

*plot_active_inactive_by_formula.png -> Predicted active/inactive compounds grouped by elemental formula.

*plot_counts_by_formula.png -> Counts of predicted actives/inactives across all candidate formulas.

*RealSampleLimoneneNegative.R -> Script for processing limonene chamber SOA data.

*Realsamplepredictions.R -> Script generating probability predictions for limonene SOA candidate structures.

*Sirius_predictionLimoneneNegative.csv -> Fingerprint predictions extracted from SIRIUS.

*xgb_final_model_trained_on_all_data.rds -> Final trained model (full dataset) applied to real sample SOA validation.

--------------------------------------------------
Validation Compounds
--------------------------------------------------
The validation_compounds/ directory contains compound-level predictions and SIRIUS outputs used for external validation of the XGBoost model.

**Each subfolder (e.g., 1_12_Hydroxyoctadecanoic_acid, 2_17_Methyltestosterone, 3_1_1_2_Trimethyl_1H_benzo_e_indole) corresponds to a single compound tested during external validation.

**Inside each compound folder:

SIRIUS output files (predicted molecular formulas and candidate structures).

Predicted fingerprints derived from MS/MS spectra.

XGBoost model predictions of carcinogenicity (probability scores and binary active/inactive classification).

**These compounds were selected by matching CAS numbers between the EPA CompTox dataset and the MassBank LCSB contributor dataset.

**The purpose of this folder is to demonstrate compound-level interpretability and reproducibility of the validation workflow.

--------------------------------------------------
Software Environment
--------------------------------------------------
*R version: 4.3.1
*Key packages: rcdk, CHNOSZ, caret, xgboost, Matrix, pROC, ggplot2, here, renv
*Reproducibility: Environment locked with `renv`.

To restore environment:
```R
install.packages("renv")
renv::restore()