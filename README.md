# Target Enrichment Analysis

## 00_prep
Install VICTOR R pacage (self-made handling outputs of a Victor plate reader), and setup path to the directories `src_dir`, `data_dir`,`rda_dir`, which contains source files, data (txt, excel) files, R object files.

## 01_MS_info
Compile compound annotation (nominal targets) and compound plate map (96 well -> 384 well) of the compound library purchased form the MicroSource Discovery Systems.

## 02_screening_info
First, compile screening data, i.e., assemble AlamarBlue cell viability, compound annotation (ID, structure), annotation (experimental conditions) and normalize viability in [compile_screening_data.r](compile_screening_data.r). Next, visualize the results to QC [checking_raw_intensity.r](checking_raw_intensity.r). Then, compute difference in AUC as a measure of effects of screening compounds in [compute_auc.r](compute_auc.r).

## 03_SEA_ChEMBL
First, compile the ligand-target relationships for each screening compounds, predicted by SEA [2013-04-10-trim-dictionary-MS.r](2013-04-10-trim-dictionary-MS.r). Next, filter the predictions by different E-value thresholds.

## 04_gsea
First, compound annotation (`dicts`) and screening data (`eff`) are combined into one object in [gsea_prep.r](gsea_prep.r). Next, run GSEA (using R package `fgsea`) in [gsea.r](gsea.r) and [master_gsea.sh](master_gsea.sh). GSEA was run in HPC using Slurm scheduler. Next, the data are assembled and summarized in plots in [assemble_gsea.r](assemble_gsea.r). Finally, overlap index is calculated to find target proteins bound by a similar set of compounds [summary_eff_scatter.r](summary_eff_scatter.r).

## 05_summary_compounds
Identify hit compounds that inhibit ferroptosis and/or necrosis in [exemplars.r](exemplars.r).
