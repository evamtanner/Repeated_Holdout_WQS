######################################################################################################
# ISEE 2019, Utrecht, Netherlands
# Preconference Workshop - Mixures Analysis with Weighted Quantile Sum regression & its Extensions
# Tutorial on WQS with Repeated Holdout Validation

# Authors: Eva Tanner & Chris Gennings
#          Icahn School of Medicine at Mount Sinai
# Contact: eva.tanner@mssm.edu

# This tutorial made using R v3.6.1 & gWQS v1.1.1
######################################################################################################

install.packages("gWQS")
install.packages("purrr")
install.packages("broom")
install.packages("matrixStats")
install.packages("plyr")
library("gWQS")
library("purrr")
library("broom")
library("matrixStats")
library("plyr")

# Provided Dataset:
SELMA <-read.csv(url("https://raw.githubusercontent.com/evamtanner/Repeated_Holdout_WQS/master/ISEE%202019%20Tutorial%20-%20Simulated%20SELMA%20Data.csv"), 
                 header=T)
## Simulated from SELMA cohort (Bornehag et al. 2012, doi:10.1111/j.1365-3016.2012.01314.x) for privacy
## Exposure: 10 suspected Endocrine Disrupting Chemicals (EDCs) measured in prenatal 
### urine (creatinine-corrected), serum, and plasma
## Outcome: IQ scores assessed at age 7
## Confounders: maternal age, IQ, and weight
## Simulations based on variable distributions and correlation structure (SAS PROC SIMNORNMAL)

EDCmix = c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10') #assign mixture variables
set.seed(213684854) 

######################################################################################################
# Run WQS with single partition per Carrico et al. (2015, doi:10.1007/s13253-014-0180-3)
## We will consider chemicals potentially concerning if weights > 0.1 (100%/10 chemicals)
split.1 = gwqs(IQ ~ mom_age + mom_IQ + mom_weight, 
               mix_name = EDCmix, 
               data = SELMA, 
               q = 5, # Quintiled exposure
               validation = 0.6, # 40% trainging data, 60% testing data 
               b = 10, # This bootstraps the training set for weight estimation
                       # 100 bootstraps recommended, but gWQS runs faster with fewer
               b1_pos = F, b1_constr = T, # We hypothesize higher prenatal EDCs related to lower IQ
               family = "gaussian", 
               plots = T, tables = T)

# Rerun WQS with a single partition. 
## Do you get similar WQS index estimates and potential chemicals of concern (weight > 0.1)?
### If results are stable, you can report results for a single partition
### If results differ (e.g. lead to different conclusions), you may need alternative methods
### Unstable results may be due to low power, unrepresentative partitions, and/or ungeneralizeable weights
######################################################################################################

######################################################################################################
# One solution to stabilize results is to train and test on the same data.
split.0 = gwqs(IQ ~ mom_age + mom_IQ + mom_weight, 
               mix_name = EDCmix, 
               data = SELMA, 
               q = 5, 
               validation = 0, # train & test on same data
               b = 10, # 100 bootstraps recommended, but gWQS runs faster with fewer
               b1_pos = F, b1_constr = T, 
               family = "gaussian", 
               plots = T, tables = T)
## Do you get similar WQS index estimates and potential chemicals of concern (weight > 0.1)?
### Unfortunately, training and testing on the same data may produce overly optimistic, 
### ungeneralizeable results (Yarkoni & Westfall 2017, doi: 10.1177/17456916117693393)
######################################################################################################

######################################################################################################
# Repeated Holdout Validation is another solution
## This bootstraps the single-partition WQS 100 times to create a distribution of validated results
rep.hold = lapply(1:100,function(i){ # 100 repetitions
                  gwqs(IQ ~ mom_age + mom_IQ + mom_weight, 
                       mix_name = EDCmix, 
                       data = SELMA, 
                       q = 5, 
                       validation = 0.6, # 40% trainging data, 60% testing data
                       b = 10, # 100 bootstraps recommended, but gWQS runs faster with fewer
                       b1_pos = F, b1_constr = T, 
                       family = "gaussian", 
                       plots = F, tables = F)})

# Function to compile results from 100 repeated holdouts
f_compile <- function(gWQSresult){ # gWQSresult = gWQS results name for repeated holdout (a list)
  
  # WQS index beta estimates
  tmp1 = map(gWQSresult,"fit") %>% map(tidy) # convert fit results to a dataframe
  tmp1 = map2(tmp1,1:length(tmp1),function(x,y){x$rep=y; return(x)}) # add repetition index number
  tmp1 = do.call("rbind",tmp1) # bind all repetitions together
  tmp1 = subset(tmp1, term=="wqs") # subset results to just the wqs-index estimate
  
  # Chemical weights
  tmp3 = map(gWQSresult,"final_weights") # convert weight results to a dataframe
  tmp3 = map2(tmp3,1:length(tmp3),function(x,y){x$rep.rep=y; return(x)}) # add repetition index number
  tmp3 = do.call("rbind",tmp3) # bind all repetitions together
  tmp3 = reshape(tmp3, timevar="mix_name", idvar="rep.rep", direction="wide") # put weights in wide format
  names(tmp3) = sub(".*\\.", "", names(tmp3)) # remove 'mean_weight.' 
  
  # Merge estimates & weights
  tmp4 = merge(tmp1,tmp3,by="rep") 
  tmp4 <- subset(tmp4, select=-c(term,statistic,std.error,p.value)) # drop unecessary variables
}

# Apply f_compile function to put repeated holdout results into dataframe
rep.hold.df = f_compile(rep.hold)

# Check the WQS index estimate distribution - if not approximately ~N then more repetitions needed
hist(rep.hold.df$estimate)

# Calculate the nonparametric WQS index estimate & 95% CI from repeated holdout sampling distributon
median(rep.hold.df$estimate) 
quantile(rep.hold.df$estimate,0.025) 
quantile(rep.hold.df$estimate,0.975)
## How does this estimate compare to the single and no partition methods?

# Calculate mean weights
colMeans(rep.hold.df[,3:12])
######################################################################################################
