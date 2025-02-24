# SQR-with-time-series-data

## Smoothing Quantile Regressions (SQR) with Time Series Data

This repository contains an implementation of smoothing quantile regression for time series data. The main objectives are:

### Simulation Study and Monte Carlo Analysis
The main script ([MCS-lambdas-general.R](MCS-lambdas-general.R)) runs a Monte Carlo simulation with 5000 replications using a custom time series model. The simulated paths are generated via quantile functions and then used to estimate regression parameters using both standard quantile regression ([quantreg](https://cran.r-project.org/web/packages/quantreg/)) and a smoothing approach via the [conquer](https://cran.r-project.org/web/packages/conquer/) package.

### Lambda Function for Smoothing
A set of predefined lambda functions are available to control the smoothing. These are defined in [choose_lambda.R](choose_lambda.R), where a choice (e.g., `"hacovercos"`) determines the transformation applied to the quantile index \(\tau\).

### Estimation and Bandwidth Tuning
The simulation computes coefficient estimates, bias, and mean squared error (MSE) for different bandwidth levels. Several arrays store the estimates, and additional diagnostic plots are generated using packages such as `ggplot2`, `cowplot`, and `latex2exp`.

### Output Files
The script saves CSV checkpoints and RDS files for the simulated coefficient arrays and produces summary plots that include the mean, bias, and variance per quantile across replications.

## Repository Main Structure

- **[choose_lambda.R](choose_lambda.R)**  
  Contains the function [`lambda(tau, choice)`](choose_lambda.R) that selects the transformation for smoothing quantile regressions.

- **[MCS-lambdas-general.R](MCS-lambdas-general.R)**  
  The main simulation script that sets up the working directories, defines the time series model parameters, simulates sample paths, performs quantile regressions, computes diagnostic metrics (bias, variance, MSE), and exports the results.

- **README.md**  
  This file.

## How to Use

1. Set your working directories in [MCS-lambdas-general.R](MCS-lambdas-general.R) to ensure that checkpoint outputs and plots are correctly saved.
2. Choose the appropriate lambda function by setting the `choice` parameter (available options are defined in [choose_lambda.R](choose_lambda.R)).
3. Run the simulation:
   
   ```r
   source("MCS-lambdas-general.R")
   ```

4. Check the output folder for CSV files, RDS objects, and generated plots summarizing the performance of the smoothing quantile regression.

This repository is intended for researchers and practitioners interested in smoothing techniques in quantile regression, particularly in the context of time series data.

Feel free to adjust any details (e.g., chosen parameter values or file paths) to match your exact experimental setup.
