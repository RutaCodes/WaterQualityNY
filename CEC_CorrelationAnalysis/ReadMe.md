*CEC_corr_analysis.R* - code that analyzes correlation patterns among samples as well as detected compounds.

*CEC_use_info_dataset.csv* - information about contaminants that were detected in samples.

Dataset with detected concentration values will be shared after this project is published.

# Correlation Analysis

The spearman rank correlation coefficient (rho) was used to analyze monotonicity of observed changes and quantify similarity in CEC occurrence patterns. Spearman rank 
correlation was chosen since it does not require the relationship between variables to be linear, as temporal CEC magnitude changes are often non-linear. First, we 
calculated spearman rank correlation coefficients among the most commonly detected compounds (detection frequency > 45%) to determine co-occurrence patterns. For this 
analysis, we combined all samples taken during two rounds of sampling. Strong positive coefficients are interpreted to indicate similarities in CEC origin as changing 
hydroclimatic conditions affected those compounds in a comparable way resulting in strong positive monotonic pattern.

We also calculated spearman rank correlation coefficients between each sample (across multiple rounds of sampling in each watershed). In this comparison, we used 
concentrations of all detected compounds shared between two samples to determine similarity between samples, and note that number of detected CEC changes from sample 
to sample. Although we expected to see clear differences among samples from different watersheds due to varying source dynamics, low correlation among samples from the 
same watershed would indicate that temporal changes in hydroclimatic conditions may have strong impact on CEC detection.  

Across both analyses, non-detect concentration values were identified using pairwise detection and removed prior to correlation coefficient calculations. Spearman 
rank correlation coefficient values were calculated using rcorr function in R (version 1.1.463), which allows the correlation values to reflect patterns among detected 
concentration without being affected by non-detects using pairwise detection. As we implemented Spearman correlation analysis for only the most commonly detected 
compounds, we further minimized the effects of non-detects on our correlation results. 
