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
setwd("C:/Users/afossa/Brown Dropbox/Alan Fossa/HOME_MIREC_metabolomics")

#Load in raw feature tables and create unique row identifier----
#neg_wide<-read.csv("data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.csv")

#saveRDS(neg_wide,file="data\\IDSL.IPA_NEG\\peak_alignment\\peak_area.rds")

neg_long<-readRDS("data\\peak_area.rds")

row.names(neg_long)<-paste(neg_long$mz,neg_long$RT,sep="__")#Set peak identifier as row names.

length(unique(row.names(neg_long)))==nrow(neg_long)#Are all peak identifiers unique? Yes!

#Drop needless columns including "PRIME" samples----
neg_thin<-neg_long %>% 
  select(!c(X,mz,RT,freqPeakArea,medianPeakArea,Flag,contains("PRIME")))

rm(neg_long)

#Filter out peaks with all zeros as sample values and set non-detects (0) to missing----

neg_thin$area_sum<-apply(neg_thin,1,sum)#Column of peak area sums

neg_thin<-neg_thin %>% 
  dplyr::filter(area_sum!=0) %>% 
  select(-area_sum) %>% 
  mutate(across(everything(),~if_else(.x==0,NA,.x)))

#Create a table of non-detect flags----
# neg_flg<-neg_thin %>% 
#   mutate(across(everything(),.fns=~if_else(is.na(.x),"Non-detect","Detect")))
# 
# saveRDS(neg_flg,"data\\neg_flg.rds")
# 
# rm(neg_flg)

#Replace non-detects with min peak area/sqrt(2) and log transform----
neg_thin$area_min<-apply(neg_thin,1,min,na.rm=T)#Column of peak area mins

neg_thin_ln<-neg_thin %>% 
  mutate(across(!area_min,~if_else(is.na(.x),log(area_min/sqrt(2)),log(.x)))) %>% 
  select(-area_min)

#Transpose dataset----
neg_wide_ln<-t(neg_thin_ln) %>% as.data.frame() 

neg_wide_ln$sample_mzml<-row.names(neg_wide_ln)#Set row names as column in dataset----

#Import meta data----
meta<-read_excel("data\\HOME_negative_batch_and_run_order.xlsx",sheet="C18NEG MS") %>% 
  dplyr::filter(!is.na(SAMPLES)) %>% 
  mutate(
    sample_mzml=str_replace(SAMPLES,".cdf",".mzML"),
    batch_num=as.numeric(str_remove(BATCH,"B"))
    ) %>% 
  arrange(batch_num,ORDER) %>% 
  mutate(order_cont=row.names(.)) %>% #Create a continuous run order (currently grouped within batch)
  select(-c(BATCH,ORDER,SAMPLES))

#Merge meta data into feature table----
neg_wide_ln_w_meta<-left_join(
  neg_wide_ln,
  meta,
  by="sample_mzml"
) %>% 
  arrange(order_cont)

#Use WaveICA2.0 to correct for batch effects----
Order<-neg_wide_ln_w_meta$order_cont

correction<-WaveICA_2.0(
  select(neg_wide_ln_w_meta,!all_of(names(meta))),
  Injection_Order=Order,
  alpha=0,
  Cutoff=0.1,
  K=10
)

corrected<-correction$data_wave %>% as.data.frame()

##Merge meta data back in----
corrected_w_meta<-bind_cols(corrected,arrange(meta,order_cont))

#Asses batch correction in peaks detected in >=20% of participant samples----

###Get names of participant samples----
part_samples<-meta[meta$CLASS=="UC",]$sample_mzml

###Load in detection table and select only participant samples----
part_samples_flg<-readRDS("data\\neg_flg.rds") %>% 
  select(all_of(part_samples))

part_samples_flg$pct_detect<-apply(part_samples_flg,1,FUN=function(x){mean(x=="Detect")})

###Get list of peaks detected in >=10% of participant samples----
peaks_pct<-row.names(part_samples_flg[part_samples_flg$pct_detect>=0.20,])

###Fit pc models----
####Raw feature table----

raw_20<-neg_wide_ln_w_meta[,peaks_pct]

pc_model_raw<-prcomp(
  raw_20,
  rank.=3,
  scale.=T
)

neg_wide_ln_w_meta[,c("PC1","PC2","PC3")]<-predict(pc_model_raw,newdata = neg_wide_ln_w_meta)

####Corrected feature table----

corrected_20<-corrected_w_meta[,peaks_pct]

pc_model_corrected<-prcomp(
  corrected_20,
  rank.=3,
  scale.=T
)

corrected_w_meta[,c("PC1","PC2","PC3")]<-predict(pc_model_corrected,newdata = corrected_w_meta)

###Vizualize batch effects----

####Raw feature table----
raw_plot_a<-neg_wide_ln_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC1,y=PC2,color=as.factor(batch_num),alpha=0.5))+
  guides(color="none")+
  theme_classic()

raw_plot_b<-neg_wide_ln_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC2,y=PC3,color=as.factor(batch_num)),alpha=0.5)+
  guides(color="none")+
  theme_classic()

raw_plot_c<-neg_wide_ln_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC1,y=PC3,color=as.factor(batch_num)),alpha=0.5)+
  guides(color="none")+
  theme_classic()

grid.arrange(raw_plot_a,raw_plot_b,raw_plot_c,nrow=2,top="Uncorrected negative mode output")

####Corrected feature table----
corrected_plot_a<-corrected_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC1,y=PC2,color=as.factor(batch_num)),alpha=0.5)+
  guides(color="none")+
  theme_classic()

corrected_plot_b<-corrected_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC1,y=PC3,color=as.factor(batch_num)),alpha=0.5)+
  guides(color="none")+
  theme_classic()

corrected_plot_c<-corrected_w_meta %>% 
  ggplot()+
  geom_point(aes(x=PC2,y=PC3,color=as.factor(batch_num)),alpha=0.5)+
  guides(color="none")+
  theme_classic()

grid.arrange(corrected_plot_a,corrected_plot_b,corrected_plot_c,nrow=2,top="Corrected negative mode output")

#Transpose corrected data back to long form, back transform log transformation, and export----
row.names(corrected_w_meta)<-corrected_w_meta$sample_mzml
  
corrected_long<-corrected_w_meta %>% 
  select(contains("__")) %>% #Select only peak columns. 
  t() %>% 
  as.data.frame() %>% 
  mutate(
    across(everything(),exp),
    mz__rt=row.names(.)
    ) 

write_rds(corrected_long,"data\\peak_area_waveICA2.0_corrected.rds")
