#Purpose: This program batch corrects raw metabolomics feature tables, creates a data table indicating whether #a feature value was originally non-detect and imputed as min value for that feature/2, and illustrative plots and summary statistics.

#Author: Alan J. Fossa (alan_fossa@brown.edu)

#Loading required packages----
library(tidyverse)
library(devtools)
library(readxl)
library(conflicted)
#devtools::install_github("dengkuistat/WaveICA_2.0",host="https://api.github.com")
library(WaveICA2.0)
library(gridExtra)

#Set working dir----
setwd("C:\\Users\\afossa\\Brown Dropbox\\Alan Fossa\\HOME_MIREC_metabolomics\\")

#Load in raw feature tables and create unique row identifier----
##Negative ion mode peak area----
neg<-read.csv("data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.csv") %>% 
  mutate(
    feature_id=paste(mz,RT,sep="__")
  ) %>% 
  select(-X)#Delete first column of rownames

#Get detection frequencies for participant samples(%)----
neg$detect_prop<-neg %>% 
  dplyr::select(contains("SAMPLE")) %>% 
  apply(MARGIN = 1, FUN = function(x){sum(x>0)/ncol(.)})

#Plot number of retained features across detection thresholds----
thresholds<-seq(from=0.10,to=0.75,by=0.05)

features_retained<-vector(mode="integer",length=length(thresholds))

for (i in seq_along(thresholds)){
  features_retained[i]<-sum(neg$detect_prop>=thresholds[i])
}

plot(thresholds,features_retained,xlab="Detection cutoff",ylab="Features retained")

#Create non-detect flag---
neg_w_flg<-neg %>% 
  mutate(across(contains("HOME_"),.fns=~if_else(.x>0,"Detect","Non-detect"),.names = "{.col}_detect"))

#Batch correct using waveICA2.0----