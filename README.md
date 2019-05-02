#Code to Implement Repeated Holdout Validation for Weighted Quantile Sum Regression
*Weighted Quantile Sum (WQS) regression is a method used to examine chemical mixtures in relation to a health outcome of interest. 
*Data is typically randomly partitioned into a single training set for estimating weights, and a single test set for assessing the association between between the weighted index and the outcome. However, this may produce unstable estimates in finite samples.
*Repeated holdout validation solves this problem by partitioning data m times (>=100) to simulate a distribution of validated results.

##Software Requirements
*R version 3.5.1 or greater
*Packages: gWQS, purrr, broom, matrixStats, plyr, xlsx

