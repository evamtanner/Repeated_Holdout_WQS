# Code to Implement Repeated Holdout Validation for Weighted Quantile Sum Regression
* Weighted Quantile Sum (WQS) regression is a method used to examine chemical mixtures in relation to a health outcome of interest. 
* Data is typically randomly partitioned into a single training set for estimating weights, and a single test set for assessing the association between between the weighted index and the outcome. However, this may produce unstable estimates in finite samples.
* Repeated holdout validation solves this problem by partitioning data m times (>=100) to simulate a distribution of validated results.

## Files 
* Repeated Holdout Validation for WQS.R - Implement & compile repeated holdout WQS results and compare to no validation R code
* Weight Uncertainty Plot.sas - Visualize uncertainty in identifying chemicals of concern SAS code
* WQS with IPWs & Repeated Holdout.R - Edited gwqs function to incorporate inverse probability weights R code
* ISEE 2019 Preconference Workshop Tutorial
* Simulated SELMA data for the ISEE 2019 Preconference Workshop

## Software Requirements
### Repeated Holdout WQS:
* R version 3.5.1 or greater
* Packages: gWQS, purrr, broom, matrixStats, plyr, xlsx
### Weight Uncertainty Plot
* SAS 9.4 (14.3)

## Authors
Eva M Tanner & Chris Gennings, Icahn School of Medicine at Mount Sinai, New York, USA

[![DOI](https://zenodo.org/badge/183286526.svg)](https://zenodo.org/badge/latestdoi/183286526)

### Collaborators
* Carl-Gustaf Bornehag, Karlstad University, Karlstad, Sweden

### Acknowledgement of Genius R Skills
* Jonathan Heiss & Anu Joshi, Icahn School of Medicine at Mount Sinai, New York, USA

### Article & Study Links

[Early prenatal exposure to suspected endocrine disruptor mixtures is associated with lower IQ at age seven](https://doi.org/10.1016/j.envint.2019.105185)

[The SELMA Study](http://selmastudien.se/)

Licensed under the [MIT License](LICENSE).
