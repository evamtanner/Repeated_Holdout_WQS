######################################################################################################
### Example of WQS with Repeated Holdout Validation for 26 EDCs & IQ in SELMA ###
######################################################################################################

#load required packages
library("gWQS")
library("purrr")
library("broom")
library("matrixStats")
library("plyr")
library("xlsx")
library("boot")

#set up basic gWQS 
setwd("J:\\PM\\SELMA\\WQS\\Analysis for Manuscript - March 2019")
data = read.csv("wqsdata_April_2019.csv", header=TRUE)

X26 = c('MEP','MBP','MBzP','DEHP','DINP','MOiNCH','MHiDP','MCiNP',
         'DPP','TCP','BPS','BPA','BPF','Triclosan','PBA','OHPH',
         'PCB','HCB','Nonachlor','DDT_DDE','PFNA','PFDA','PFUnDA','PFHxS','PFOA','PFOS')

set.seed(213684854)

######################################################################################################
#Repeated Holdout
###use lapply function to repeat gWQS desired #times
######################################################################################################
IQ_rep = lapply(1:100,function(i){
  gwqs(IQ ~ mom_edu + mom_IQ + mom_smoke_3 + mom_weight + birth_week + Male,
       mix_name = X26, data = data, 
       q = 10, validation = 0.6, b = 100, b1_pos = F, b1_constr = T, 
       family = "gaussian", wqs2 = F, plots = F, tables = F)
})

######################################################################################################
#Organize Repeated holdout Validation Results
###Create a compile function to compile wqs-index estimate and weight results from N repetitions  
######################################################################################################

f_compile <- function(gWQSresult){ #gWQSresult = gWQS results name for repeated holdout (a list)
  
  #WQS index beta estimates
  tmp1 = map(gWQSresult,"fit") %>% map(tidy) #convert fit results to a dataframe
  tmp1 = map2(tmp1,1:length(tmp1),function(x,y){x$rep=y; return(x)}) #add repetition index number
  tmp1 = do.call("rbind",tmp1) #bind all repetitions together
  tmp1 = subset(tmp1, term=="wqs") #subset results to just the wqs-index estimate
  
  #WQS index IQR
  ##Note that N repetitions must be changed manually for this part to work
  ###this shows example for 100 repetitions
  tmp2 = map(gWQSresult,"wqs") 
  tmp2 = matrix(unlist(tmp2), nrow=100, byrow=T) #convert wqs index values (for M observations) to a matrix
  tmp = data.frame(cbind(matrix(1:100, nrow=100, ncol=1),rowIQRs(tmp2))) #convert matrix to dataframe & calculate IQR
  tmp2 = rename(tmp, c('X1'='rep','X2'='wqs_iqr'))
  
  #Chemical weights
  tmp3 = map(gWQSresult,"final_weights") #convert weight results to a dataframe
  tmp3 = map2(tmp3,1:length(tmp3),function(x,y){x$rep.rep=y; return(x)}) #add repetition index number
  tmp3 = do.call("rbind",tmp3) #bind all repetitions together
  tmp3 = reshape(tmp3, timevar="mix_name", idvar="rep.rep", direction="wide") #put weights in wide format
  names(tmp3) = sub(".*\\.", "", names(tmp3)) #remove 'mean_weight.' 
  
  #Merge estimates, IQR, weights
  tmp4 = merge(tmp1,tmp2,by="rep")
  tmp4$est_iqr = tmp4$estimate*tmp4$wqs_iqr #calculate interpreted estimate (change in y per IQR change in WQS index)
  tmp5 = merge(tmp4,tmp3,by="rep") 
  tmp5 <- subset(tmp5, select=-c(term,statistic,std.error,p.value)) #drop unecessary variables
}

######################################################################################################

#apply the f_compile function to compile repeated holdout results
IQ_100 = f_compile(IQ_rep)

######################################################################################################
#Compile Summary Statistics from Full Dataset and Repeated Holdout Results
###Output Results from Repeated Holdout and All Summary Stats
######################################################################################################

f_summary <- function(gWQSresult,FullResult,ValResult,xlsxFileName,xlsxSheetNameEst,xlsxSheetNameWeight,xlsxSheetNameValWeights){ 
  #gWQSresult = gWQS results name for repeated holdout (a list)
  #FullResult = gWQS results name for full sample (a list)
  #ValResult = compiled repeated holdout results name (a dataframe)
  #xlsxFileName = excel file name
  #xlsxSheetNameEst = excel file sheet to summary of estimate results
  #xlsxSheetNameWeight = excel file sheet to summary of weight results
  #xlsxSheetNameValWeights = excel file sheet name for compiled repeated holdout weight results
  
  #Full Dataset
  tmp1 = tidy(FullResult$fit) #convert full fit into dataframe
  tmp1 = rename(tmp1, c('estimate'='est_full')) 
  tmp1$LL_full = tmp1$est_full - (1.96*tmp1$std.error) #calculate 95% CI
  tmp1$UL_full = tmp1$est_full + (1.96*tmp1$std.error)
  
  tmp1$est_full_iqr = tmp1$est_full*IQR(FullResult$wqs) #calculate interpretted estimate
  tmp1$LL_full_iqr = tmp1$LL_full*IQR(FullResult$wqs) #calculate 95% CI for interpretted estimate
  tmp1$UL_full_iqr = tmp1$UL_full*IQR(FullResult$wqs)
  
  tmp1 = subset(tmp1, term=="wqs") #subset results to just the wqs-index estimate
  tmp1 = subset(tmp1, select=-c(term,std.error,statistic,p.value)) #drop unecessary variables
  
  tmp2 = as.data.frame(FullResult$final_weights) #convert weights to dataframe
  tmp2 = rename(tmp2,c("mean_weight"="weight_full"))
  
  #Validated Results
  tmp1$est_val_mean = mean(ValResult$estimate) #estimate based on mean of sampling distribution
  tmp1$LL_val_mean = tmp1$est_val_mean - (1.96*(sd(ValResult$estimate))) #symmetric 95% CI
  tmp1$UL_val_mean = tmp1$est_val_mean + (1.96*(sd(ValResult$estimate)))
  
  tmp1$est_val_med = median(ValResult$estimate) #estimate based on median of sampling distribution
  tmp1$est_val_LL = quantile(ValResult$estimate,0.025) #assymmetric 95% CI
  tmp1$est_val_UL = quantile(ValResult$estimate,0.975)
  
  tmp1$est_val_mean_iqr = mean(ValResult$est_iqr) #mean interpretted estimate of sampling distribution
  tmp1$LL_val_mean_iqr = tmp1$est_val_mean_iqr - 1.96*(sd(ValResult$est_iqr))
  tmp1$UL_val_mean_iqr = tmp1$est_val_mean_iqr + 1.96*(sd(ValResult$est_iqr))
  
  tmp1$est_val_med_iqr = median(ValResult$est_iqr) #median interpretted estimate of sampling distribution
  tmp1$est_val_LL_iqr = quantile(ValResult$est_iqr,0.025)
  tmp1$est_val_UL_iqr = quantile(ValResult$est_iqr,0.975)
  
  tmp3 = data.frame(weight_rep=colMeans(ValResult[,5:30])) #this is repetitive, but needs weights in long format for plots, must change based on #X
  labels = rownames(tmp3)
  as.data.frame(labels)
  tmp3 = data.frame(labels,tmp3)
  tmp4 = data.frame(times_bad=colSums(ValResult[,5:30] > (1/length(X26))))
  labels = rownames(tmp4)
  as.data.frame(labels)
  tmp4 = data.frame(labels,tmp4)
  tmp3 = merge(tmp3,tmp4,by="labels") #merge stats from val
  tmp3 = rename(tmp3,c("labels"="mix_name"))
  
  #merge weight stats
  tmp2 = merge(tmp2,tmp3,by="mix_name")
  
  #weights in long format for plots
  tmp5 = map(gWQSresult,"final_weights") #convert weight results to a dataframe
  tmp5 = map2(tmp5,1:length(tmp5),function(x,y){x$rep.rep=y; return(x)}) #add repetition index number
  tmp5 = do.call("rbind",tmp5) #bind all repetitions together
  tmp6 = rename(tmp5, c('rep.rep'='rep'))
  
  #Add to Excel results file
  write.xlsx(tmp1, file=xlsxFileName, sheetName=xlsxSheetNameEst, row.names=F, col.names=T, append=T)
  write.xlsx(tmp2, file=xlsxFileName, sheetName=xlsxSheetNameWeight, row.names=F, col.names=T, append=T)
  write.xlsx(tmp6, file=xlsxFileName, sheetName=xlsxSheetNameValWeights, row.names=F, col.names=T, append=T)
}

######################################################################################################

#apply the f_summary function to summarize results from the full sample & repeated holdout 
##& output all results to to Excel
####f_summary(FullResult,ValResult,xlsxSheetNameEst,xlsxSheetNameWeight,xlsxSheetNameVal)
setwd("J:\\PM\\SELMA\\WQS\\Analysis for Manuscript - March 2019\\Results")
    
f_summary(IQ_rep,IQ_full,IQ_100,"WQS-IQ X26 Smoke3 4-16-19.xlsx","IQest","IQweight","IQ100")






