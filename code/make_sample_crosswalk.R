#Purpose: This R script takes raw sequence and batch files from mass spec analysis and creates crosswalk files that can be used to link sample labels in metabolomics feature tables to HOME and MIREC participant study data. 

#Author: Alan J. Fossa (alan_fossa@brown.edu)

#Loading required packages----
library(tidyverse)
library(readxl)
library(conflicted)

#Set working dir----
setwd("C:\\Users\\afossa\\Brown Dropbox\\Alan Fossa\\HOME_MIREC_metabolomics")

#HOME cohort negative mode----

##List of sequence table files----
home_neg_seq_files<-paste(getwd(),"data\\ID -Batch Files",list.files("data/ID -Batch Files",pattern="C18NEG"),sep="/")

##Import and concatenate tables and make useful variable names----
home_neg_seq_table<-map_dfr(neg_seq_files,read_csv) %>% 
  rename(
    sample_label=Filename,
    barcode=`Sample name`
  ) %>% 
  select(sample_label,barcode)

##Load in crosswalk between barcode and participant ID from batch files----
home_neg_batch_table<-read_excel("data\\HOME Study Batches.20230524.xlsx") %>% 
  rename_all(tolower) %>% 
  rename(
    timepoint=sample_type
  )

##Merge sequence and batch data and export to flat file----  
home_neg_crosswalk<-inner_join(
  home_neg_seq_table,
  home_neg_batch_table,
  by="barcode"
)

write_csv(home_neg_crosswalk,"data\\home_neg_crosswalk.csv")
