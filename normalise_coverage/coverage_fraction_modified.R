
library(tidyverse)
library(cowplot)
library(magick)
library(ggpubr)

#######
#Set the dir to the strain of interest.
#######
setwd("~/Documents/natural_isolates_methylation/")
########

unfiltered_rds_filename <- "results/nd_merged/E1_1_37_difference.RDS" #commandArgs(trailingOnly=T)[1]

motif_sites_filename<- "data/motif_locations/DAM_E1_positions.tab" #commandArgs(trailingOnly = T)[2]

unfiltered_rds <- readRDS(unfiltered_rds_filename)
motif_sites <- read.table(motif_sites_filename)

filenames<-as_vector(unfiltered_rds_filename)%>%
  str_remove("_difference.RDS")%>%
  str_remove("results/nd_merged/")

motif<-motif_sites_filename%>%
  str_remove("data/motif_locations/")%>%
  str_remove("_\\w{2}_positions.tab")

strain<-as_vector(unfiltered_rds_filename)%>%
  str_remove("_\\w{1}_\\w{2}_difference.RDS")%>%
  str_remove("results/nd_merged/")


motif_sites<-motif_sites%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%  
  filter(Contig == "contig_1_pilon")%>%
  mutate(Site=as.numeric(Site))%>%
  drop_na(Site)%>%
  add_column(Pval=0)%>%
  add_column(WGA_cov=0)%>%
  add_column(Nat_cov=0)

#Filter the rds files to the first 1million reads 
clean_rds<-unfiltered_rds%>%
  mutate(u_test_pval_log=log10(u_test_pval))%>%
  filter(contig == "contig_1_pilon")%>%
  filter(!is.na(u_test_pval))

# and a window size to look at to see if there are modified sites
# this window is *exclusive* so a window of 2 will look 1bp up and 
# 1 bp downstream
window <- 2


#This was mainly used on the pre filtered RDS sites, hence it was important to check the site location, In this case all positions sites will be present. 
#This actually means, for some sites a slightly different sites p values may have been used, say if the exact site did not pass the filter, but its neighbor did. 

for (i in 1:length(motif_sites$Site)) {
  
  #i wanted to make sure this was done in the same way as p values were selected, so that coverage points would correspond with the p values that were used. 
  #Make a mini dataframe of the RDS sites which match the position of interest
  if (dim(subset(clean_rds, position < (motif_sites[i,2]+ window) & position > (motif_sites[i,2]-window)))[1]>0){
   
    df<-subset(clean_rds, position < (motif_sites[i,2]+ window) & position > (motif_sites[i,2]-window))
    #Filter out the other contigs to only use contig one
    df<- df%>%
      filter(contig=="contig_1_pilon")
    
    #As with the P values, find the absolute distance from the exact location for each site in the window
    pos.abs<-abs(df$position-motif_sites[i,2])
    
    #Find which sites are the closest
    pos.min<- which(pos.abs==min(pos.abs))
    
    #Get the P values for the sites which are the closest
    p.values<-df[pos.min,"u_test_pval_log"]
    
    #Take the lowest P value from among the closest sites
    p.value.min<-which(p.values==min(p.values))
    
    #Get the index for the pvalue which was used. 
    p.value.used<-pos.min[p.value.min]
    motif_sites[i,3] <- df[pos.min[p.value.min],"u_test_pval_log"]
    #Use this index to take the coverage values from the df
    wga_cov<-df[p.value.used,"N_wga"]
    nat_cov<-df[p.value.used,"N_nat"]
    #Save the coverage values to the motif sites DF. 
    motif_sites[i,4]<-wga_cov
    motif_sites[i,5]<-nat_cov  
  } 

}

#Save a copy of the motif sites and Pvals as generating this takes a long time.
write_delim(motif_sites, paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",filenames,"_",motif,"_coverage.tab"), delim = "\t")




