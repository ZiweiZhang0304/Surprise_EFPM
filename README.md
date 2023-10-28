# Surprise_CFPM

# Code for running the Surprise CFPM Pipeline

Code associated with building co-fluctuation based predictive models. Helpful functions are in /CFPM/func

Main script containing the code for building CFPM and obtaining CFPM strength in external data is /CFPM/CFPM.m. This script performs the following steps:

1. Load example edge time series and behavioral data.
2. Generate correlation matrix between behavioral variable of interest and edge time series.
3. Perform edge selection based on correlation matrix using cross validation.
   - Significance of selected edges can be assessed by comparing the Rho value on held-out fold with permutation-generated null.
4. Obtain selected edges for the CFPM as positive and negative masks.
5. Apply CFPM masks to external dataset and calculate CFPM strength time series.
   - To test external generalizability of the CFPM, you can then select a statistical method to assess the relationship between the CFPM strength time series and the time series of the behavioral variable of interest for your own research question.
