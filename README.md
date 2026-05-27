Fedrug
A computational method for predicting ferroptosis-inducing compounds based on LINCS 2020 data.

Data Download
The raw data used in this study are from the LINCS L1000 project. Due to large file sizes and access restrictions, the data must be downloaded directly from the official source.
Source
LINCS 2020 Level 3 Data
 Website: https://clue.io/releases/data-dashboard
 Required Files:
   `level3_beta_trt_cp_n1805898x12328.gctx` (Treatment profiles)
   `level3_beta_ctl_n188708x12328.gctx` (Control profiles)

Download Instructions
1.  Register/Login: Visit the [LINCS Data Dashboard](https://clue.io/releases/data-dashboard) and register for a free account (required for data access).
2.  Locate Files: After logging in, navigate to the "L1000" or "Level 3" data section. Use the search/filter function to find the two `.gctx` files listed above.
3.  Download: Download both files. These files are very large (multiple GBs) and may require a stable internet connection.
4.  Organize: After downloading, place both `.gctx` files into the `data/` folder within your working directory.


Gene Lists
The five ferroptosis core gene sets (G1–G5) are provided in `Genelist_use.Rdata`.

load("Genelist_use.Rdata")


Analysis Scripts
Two R scripts are provided to process the data and calculate Fedrug scores.

Script 1: 1.linc2020_GSEA.R
Purpose: Process LINCS 2020 data and perform GSEA for each compound.

Input:
data/level3_beta_trt_cp_n1805898x12328.gctx (Treatment data)
data/level3_beta_ctl_n188708x12328.gctx (Control data)
Genelist_use.Rdata (Five core gene sets)

Output:
GSEA results for each compound across G1–G5 (NES and P values)

Script 2: 2.score.R
Purpose: Calculate Fedrug scores for compounds using the ridge regression model.

Input:
GSEA results from Script 1

Output:
Fedrug scores for each compound
Fedrug-score formula:

Fedrug score = -0.1718 + 0.3841×G1 + 0.3314×G2 + 0.3418×G3 + 0.3776×G4 + 0.3540×G5

where G1–G5 are NES values from GSEA.

For questions or collaboration: dianyating@outlook.com
