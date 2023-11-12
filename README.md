# Surprise EFPM

# Code for running the EFPM Pipeline

Zhang, Z., Rosenberg, M. D. (in prep). Brain network dynamics predict moments of surprise across contexts.

Code associated with building the edge-fluctuation based predictive models (EFPM). Helpful functions are in /EFPM/func

Main script containing the code for building EFPM and obtaining EFPM strength in external data is /EFPM/EFPM.m. This script performs the following steps:
1. Load example edge time series and behavioral data.
2. Generate correlation matrix between behavioral variable of interest and edge time series.
3. Perform edge selection based on correlation matrix using cross validation.
   - Significance of selected edges can be assessed by comparing the Rho value on held-out fold with permutation-generated null.
4. Obtain selected edges for the EFPM as positive and negative masks.
5. Apply EFPM masks to external dataset and calculate EFPM strength time series.
   - To test external generalizability of the EFPM, you can then select a statistical method to assess the relationship between the EFPM strength time series and the time series of the behavioral variable of interest for your own research question.
