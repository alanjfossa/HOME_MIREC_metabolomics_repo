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
#neg_wide<-read.csv("data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.csv")

#saveRDS(neg_wide,file="data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.rds")

neg_long<-readRDS("data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.rds")

row.names(neg_long)<-paste(neg_long$mz,neg_long$RT,sep="__")#Set peak identifier as row names.

length(unique(row.names(neg_long)))==nrow(neg_long)#Are all peak identifiers unique? Yes!

#Drop needless columns including "PRIME" samples----
neg_thin<-neg_long %>% 
  select(!c(X,mz,RT,freqPeakArea,medianPeakArea,Flag,contains("PRIME")))

rm(neg_long)

gc()

#Filter out peaks with all zeros as sample values and set non-detects (0) to missing----

neg_thin$area_sum<-apply(neg_thin,1,sum)#Column of peak area sums

neg_thin<-neg_thin %>% 
  dplyr::filter(area_sum!=0) %>% 
  select(-area_sum) %>% 
  mutate(across(everything(),~if_else(.x==0,NA,.x)))

gc()

#Transpose dataset,replace non-detect with min/sqrt(2), and log transform----
neg_wide<-t(neg_thin) %>% 
  as.data.frame() %>% 
  mutate(across(everything(),~if_else(is.na(.x),log(min(.x,na.rm=T)/sqrt(2)),log(.x)),.names = "ln_{.col}"))

#Import meta data----
meta<-read_excel("data\\Home Study BatchCorrectionFinal..xlsx",sheet="C18NEG MS") %>% 
  dplyr::filter(!is.na(SAMPLES)) %>% 
  mutate(sample_mzml=str_replace(SAMPLES,".cdf",".mzML")) %>% 
  arrange(SAMPLES)

table(meta$CLASS)

# #Get detection frequencies for participant samples(%)----
# detect_prop<-neg %>% 
#   dplyr::select(contains("SAMPLE")) %>% 
#   apply(MARGIN = 1, FUN = function(x){mean(x>0)})
# 
# #Plot number of retained features across detection thresholds----
# thresholds<-seq(from=0.10,to=0.75,by=0.05)
# 
# features_retained<-vector(mode="integer",length=length(thresholds))
# 
# for (i in seq_along(thresholds)){
#   features_retained[i]<-sum(neg$detect_prop>=thresholds[i])
# }
# 
# png("HOME_MIREC_metabolomics_repo\\output\\features_retained.png",type="cairo")
# 
# plot(thresholds,features_retained,xlab="Detection cutoff",ylab="Features retained")
# 
# dev.off()

#Create non-detect flag---
neg_w_flg<-neg %>% 
  mutate(across(contains("HOME_"),.fns=~if_else(.x>0,"Detect","Non-detect"),.names = "{.col}_detect"))

#Batch correct using waveICA2.0----

