# Using horseshoe prior for incorporating multiple historical control data in randomized controlled trials

This repository contains the code used to completion the case study in a manuscript (being under peer review) entitled "Using horseshoe prior for incorporating multiple historical control data in randomized controlled trials" by Tomohiro Ohigashi, Kazushi Maruo, Takashi Sozu and Masahiko Gosho.
The following files are included:

- `BIN_case_study.R`: Code for completing a clinical trial example with binary endpoint reported in manuscript.

- `Surv_case_study.R`: Code for completing a clinical trial example with time-to-event endpoint. Note – since the raw data for the example reported in manuscript is not publicly available, this code generates a simulated dataset and analyzes it.

- `Surv_cov_case_study.R`: Code for completing a clinical trial example with time-to-event endpoint with or without cthe ovariate adjustment. Note – since the raw data for the example reported in manuscript is not publicly available, this code generates a simulated dataset and analyzes it.

- `Stan_BIN`, `Stan_Surv`, `Stan_Surv_Cov`: Directory containing stan codes called in each sample code.

## Software Dependencies
All code and analyses generated in R version 4.1.1 and used Stan via `cmdstnr` package (version 0.4.0).

