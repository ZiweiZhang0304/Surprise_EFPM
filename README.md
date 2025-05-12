# Surprise EFPM
Zhang, Z., & Rosenberg, M. D. (2024). Brain network dynamics predict moments of surprise across contexts. Nature Human Behaviour, 1–15. https://doi.org/10.1038/s41562-024-02017-0

## Part I: Surprise edge-fluctuation based predictive models (EFPM) masks
The mask for high and low surprise network can be found under EFPM/output/masks/surprise_EFPM <br> These masks can be applied to your own brain data to obtain the strength of this model. The function EFPM/func/get_EFPM_score.m (See Part II 5.) achieves this.

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
  

All scripts were written, run, and tested in matlab 2021b. No non-standard hardware is required. Typical install time will take less than a minute. Expected run time for demo takes about 10 mins.
