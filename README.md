# Surprise EFPM
Zhang, Z., & Rosenberg, M. D. (2023). Brain network dynamics predict moments of surprise across contexts (p. 2023.12.01.569271). bioRxiv. https://doi.org/10.1101/2023.12.01.569271

## Part I: Surprise edge-fluctuation based predictive models (EFPM) masks
The mask for high and low surprise network can be found under EFPM/output/masks/surprise_EFPM <br> You can apply them to your own brain data and calculate the strength of this model with the function under EFPM/func/get_EFPM_score.m

## Part II: Code for running the EFPM pipeline
Main script containing the code for building EFPM and obtaining EFPM strength in external data is /EFPM/EFPM.m. <br>
Helpful functions are in /EFPM/func <br>
This script performs the following steps:
1. Load example edge time series and behavioral data.
2. Generate correlation matrix between behavioral variable of interest and edge time series.
3. Perform edge selection based on correlation matrix using cross validation.
   - Significance of selected edges can be assessed by comparing the Rho value on held-out fold with permutation-generated null.
4. Obtain selected edges for the EFPM as positive and negative masks.
5. Apply EFPM masks to external dataset and calculate EFPM strength time series.
   - To test external generalizability of the EFPM, you can then select a statistical method to assess the relationship between the EFPM strength time series and the time series of the behavioral variable of interest for your own research question.
