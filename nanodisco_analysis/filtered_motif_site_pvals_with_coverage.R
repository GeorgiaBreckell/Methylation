
library(tidyverse)
library(cowplot)
library(magick)

# data files
setwd("~/Documents/natural_isolates_methylation/")

filtered_rds_filename <-  commandArgs(trailingOnly=T)[1] #"results/nanodisco/filtered_rds/F6/F6_2_37_DAM_filtered_rds.tab" #

motif_sites_filename<-  commandArgs(trailingOnly = T)[2] #"data/motif_locations/DAM_F6_positions.tab"  #

filtered_rds <- read_delim(filtered_rds_filename,"\t", escape_double = FALSE, trim_ws = TRUE)
motif_sites <- read.table(motif_sites_filename)

filename<-as_vector(filtered_rds_filename)%>%
  str_remove("results/nanodisco/filtered_rds/\\w{2}/")%>%
  str_remove("_filtered_rds.tab")

strain<-as_vector(filtered_rds_filename)%>%
  str_remove("/\\w{2}_\\w{1}_\\w{2}_\\w*_filtered_rds.tab")%>%
  str_remove("results/nanodisco/filtered_rds/")

#Filter motif sites of NA values, rename columns, add extra columns for output later  
motifs_sites_clean<-motif_sites%>%
  rename(Contig=V1)%>%
  rename(Site=V2)%>%
  mutate(Contig=gsub("contig_1_pilon","contig_1",Contig))%>%
  mutate(Contig=gsub("contig_2_pilon","contig_2",Contig))%>%
  mutate(Contig=gsub("contig_3_pilon","contig_3",Contig))%>%
  mutate(Contig=gsub("contig_6_pilon","contig_6",Contig))%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)%>%
  add_column(wga_cov=0)%>%
  add_column(nat_cov=0)

# add a window size to look at to see if there are modified sites
# this window is *exclusive* so a window of 2 will look 1bp up and 
# 1 bp downstream
window <- 2

#Add a vector to store the fraction values for each strain 
motif_fraction_test<-vector()

for (i in 1:length(motifs_sites_clean$Site)) {
  
  if(dim(subset(filtered_rds, position < (motifs_sites_clean[i,2]+ window) & position > (motifs_sites_clean[i,2]-window)))[1]>0)  {
    # if it does, we switch the bit to indicate it's modified
    motifs_sites_clean[i,3] <- 1
    
    
    #Then find the P.value of matching sites 
    #record the matching sites and there meta data
    pos.mod<-subset(filtered_rds, position < (motifs_sites_clean[i,2]+window) & position > (motifs_sites_clean[i,2]-window))
    #find the absolute difference from a dam site to a modified site
    pos.abs<-abs(pos.mod$position-motifs_sites_clean[i,2])
    #find the closest site within the window
    pos.min<-which(pos.abs==min(pos.abs))
    #get the P.values of the closest site/s
    p.values<-pos.mod[pos.min,"u_test_pval_log"]
    #Take the lowest p.value from the closest sites 
    p.value.min<-which(p.values==min(p.values))
    #Save the P.value with the site. 
    motifs_sites_clean[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
    no_zero<-subset(motifs_sites_clean,Pval!=0)
    #get the coverage values for the site
    wga_cov<-pos.mod[pos.min,"N_wga"]
    nat_cov<-pos.mod[pos.min,"N_nat"]
    motifs_sites_clean[i,5]<-wga_cov
    motifs_sites_clean[i,6]<-nat_cov
  }
}


write_delim(motifs_sites_clean,paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",filename,"_site_pvals.tab"),delim = "\t")

#motif_fraction_test <- sum(motifs_sites_clean[,3])/length(motifs_sites_clean[,3])

#write(motif_fraction_test,paste0("results/nanodisco/fraction_modified/",strain,"/",filename,"_fraction_modified.tab"))


