# M-quantile Area-level Models for Robust Small Area Estimation without Reliance on Design-based Variances
### María Bugallo (mbugallo@umh.es)
Center of Operations Research - CIO, Miguel Hernández University of Elche. Address: Edificio Torretamarit. Avenida de la Universidad, s/n 03202 Elche (Alicante), Spain.
### María José Lombardía (maria.jose.lombardia@udc.es) and Alexandro Aneiros-Batista (alexandro.aneiros.batista@udc.es)
Research Center on Information and Communication Technologies, University of A Coruña. Address: Campus de Elviña s/n, 15071 A Coruña (A Coruña), Spain.

# Short description
This repository contains R scripts to perform Monte Carlo simulations comparing small area estimation methods, including Hajek, FH, robust FH, MQ, FHMQ and AMQ models, with evaluation via RRMSE and RBIAS.

# Scripts
- **`AreaLevelMQ.R`** – Main script for running the simulations. Generates population data, performs sampling, computes all small area predictors, and calculates analytical and bootstrap MSE estimators for the AMQ models.

- **`readMSEresults.R`** – Script for reading previously saved MSE results in the fold `MSEresultsB`, summarizing performance metrics (RRMSE, RBIAS), and producing boxplots for visual comparison across models and predictors.

- **`AuxFunctions.R`** – Auxiliary functions used in the main script, including prediction, MSE computations, covariance estimation, grid fitting, and MAD calculation.

# Reproducibility
- Key implementation details are documented in the scripts.
- Simulation outputs and MSE results are saved as CSV files to allow full replication of the analyses.

# Usage
1. Run `AreaLevelMQ.R` to perform the simulation experiments.
2. Use `readMSEresults.R` to read the outputs and generate plots.
3. `AuxFunctions.R` is required and automatically sourced by the main script.
