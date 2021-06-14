
library(tidyverse)
library(cowplot)
library(magick)

#######
#Set the dir to the strain of interest.
#######
setwd("~/Documents/natural_isolates_methylation/results/nanodisco/p_value_cutoff_subsample/H8/")
########
###Change these to match the strain of interest.
#######
dam <- read.table(file="~/Documents/natural_isolates_methylation/data/motif_locations/DAM_SC11A_H8_positions.tab")
dcm <- read.table(file="~/Documents/natural_isolates_methylation/data/motif_locations/DCM_SC11A_H8_positions.tab")


# data files
nd_results<-list.files(pattern = "*.RDS" )
nd_list<- lapply(nd_results,readRDS)

filenames<-as_vector(nd_results)%>%
  str_remove("_difference.RDS")


#filter the rds files to the first 1million reads in a new list
firstmil_rds<- list() 

for (i in 1:length(nd_list)){
  firstmil_rds[[i]]<-nd_list[[i]]%>%
    mutate(u_test_pval_log=log10(u_test_pval))%>%
    subset(position < 1000000)%>%
    filter(contig != "contig_2_pilon")%>%
    filter(!is.na(u_test_pval))
}

#To test only Dam sites in first 1mil bases, filter same as damc but to 1mil
dam1m<-dam%>%
  #separate(V1,c("Contig","Site"),sep="_pilon_")%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%  
  filter(Contig != "contig_2")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)

dcm1m<-dcm%>%
  #separate(V1,c("Contig","Site"),sep="_pilon_")%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%
  filter(Contig != "contig_2_pilon")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)



#Filter dam sites to only include chromosome and add site column, and remove NA values 
#damc<-dam%>%
#  separate(V1,c("Contig","Site"),sep="_pilon_")%>%
#  filter(Contig != "contig_2")%>%
#  filter(!is.na(Site))%>%
#  add_column(Site_mod=0)%>%
#  add_column(Pval=0)




#Create a list of random sites to create a null for dam and dcm 
damr<-dam%>%
  #separate(V1,c("Contig","Site"),sep="_pilon_")%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%  
  filter(Contig != "contig_2")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  mutate(Site=sample(1:1e6,length(Site)))%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)

dcmr<-dcm%>%
  #separate(V1,c("Contig","Site"),sep="_pilon_")%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%
  filter(Contig != "contig_2")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  mutate(Site=sample(1:1e6,length(Site)))%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)

#Run for both dam and dcm 

for (i in 1:length(dcmr$Site)) {
  
  if(dim(subset(firstmil_rds[[1]], position < (dcmr[i,2]+window) & position > (dcmr[i,2]-window)))[1]>0)  {
    # if it does, we switch the bit to indicate it's modified
    dcm1m[i,3] <- 1
    
    
    #Then find the P.value of matching sites 
    #record the matching sites and there meta data
    pos.mod<-subset(firstmil_rds[[1]], position < (dcmr[i,2]+window) & position > (dcmr[i,2]-window))
    #find the absolute difference from a dam site to a modified site
    pos.abs<-abs(pos.mod$position-dcmr[i,2])
    #find the closest site within the window
    pos.min<-which(pos.abs==min(pos.abs))
    #get the P.values of the closest site/s
    p.values<-pos.mod[pos.min,"u_test_pval_log"]
    #Take the lowest p.value from the closest sites 
    p.value.min<-which(p.values==min(p.values))
    #Save the P.value with the site. 
    dcmr[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
    no_zero<-subset(dcmr,Pval!=0)
  }
}


for (i in 1:length(damr$Site)) {
  
  if(dim(subset(firstmil_rds[[1]], position < (damr[i,2]+window) & position > (damr[i,2]-window)))[1]>0)  {
    # if it does, we switch the bit to indicate it's modified
    dam1m[i,3] <- 1
    
    
    #Then find the P.value of matching sites 
    #record the matching sites and there meta data
    pos.mod<-subset(firstmil_rds[[1]], position < (damr[i,2]+window) & position > (damr[i,2]-window))
    #find the absolute difference from a dam site to a modified site
    pos.abs<-abs(pos.mod$position-damr[i,2])
    #find the closest site within the window
    pos.min<-which(pos.abs==min(pos.abs))
    #get the P.values of the closest site/s
    p.values<-pos.mod[pos.min,"u_test_pval_log"]
    #Take the lowest p.value from the closest sites 
    p.value.min<-which(p.values==min(p.values))
    #Save the P.value with the site. 
    damr[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
    no_zero<-subset(damr,Pval!=0)
  }
}



#Sort these by pvalue after creating for plotting later 

damr_test<-damr%>%
  arrange(Pval)

dcmr_test<-dcmr%>%
  arrange(Pval)


# and a window size to look at to see if there are modified sites
# this window is *exclusive* so a window of 2 will look 1bp up and 
# 1 bp downstream
window <- 1

# scoot through the whole set of sites

#filter_p_values_0_0001 <- function(nanodisco_rds){
#  filtered_list <-vector('list',length(nanodisco_rds))

#  for(i in seq_along(nanodisco_rds)){
#    filtered_list[[i]]<-nanodisco_rds[[i]]%>%
#      mutate(u_test_pval_log=log10(u_test_pval))%>%
#      filter(u_test_pval_log< -4.0)
#  }
#  filtered_list
#}

#Function to check if sites match between the nanodisco output and known motif sites, and store the associated p value 
#rds will be a list of rds files for each sample, motif_sites is the DAM or DCM sites
find_motif_pval_dam <- function(rds,motif_sites) {
  
  for (i in 1:length(motif_sites$Site)) {
    if(dim(subset(rds, position < (motif_sites[i,2]+window) & position > (motif_sites[i,2]-window)))[1]>0){
      motif_sites[i,3] <- 1
      
      
      #Then find the P.value of matching sites 
      #record the matching sites and there meta data
      pos.mod<-subset(rds, position < (motif_sites[i,2]+window) & position > (motif_sites[i,2]-window))
      #find the absolute difference from a dam site to a modified site
      pos.abs<-abs(pos.mod$position-motif_sites[i,2])
      #find the closest site within the window
      pos.min<-which(pos.abs==min(pos.abs))
      #get the P.values of the closest site/s
      p.values<-pos.mod[pos.min,"u_test_pval_log"]
      #Take the lowest p.value from the closest sites 
      p.value.min<-which(p.values==min(p.values))
      #Save the P.value with the site. 
      motif_sites[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
      no_zero<-subset(motif_sites,Pval!=0)
      dam_motif_pvals<-motif_sites
    }
    
  }
  dam_motif_pvals 
}



#same funtion but for DCM motifs
find_motif_pval_dcm <- function(rds,motif_sites) {
  
  for (i in 1:length(motif_sites$Site)) {
    if(dim(subset(rds, position < (motif_sites[i,2]+window) & position > (motif_sites[i,2]-window)))[1]>0){
      motif_sites[i,3] <- 1
      
      
      #Then find the P.value of matching sites 
      #record the matching sites and there meta data
      pos.mod<-subset(rds, position < (motif_sites[i,2]+window) & position > (motif_sites[i,2]-window))
      #find the absolute difference from a dam site to a modified site
      pos.abs<-abs(pos.mod$position-motif_sites[i,2])
      #find the closest site within the window
      pos.min<-which(pos.abs==min(pos.abs))
      #get the P.values of the closest site/s
      p.values<-pos.mod[pos.min,"u_test_pval_log"]
      #Take the lowest p.value from the closest sites 
      p.value.min<-which(p.values==min(p.values))
      #Save the P.value with the site. 
      motif_sites[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
      no_zero<-subset(motif_sites,Pval!=0)
      dcm_motif_pvals<-motif_sites
    }
    
  }
  dcm_motif_pvals 
}

#Create the lists to save everything to. 
dam_p_value_distribution_1m<-list()
dcm_p_value_distribution_1m<-list()


#function_test<-find_motif_pval(test_rds[[i]],dam1m)

#loop through the function for every file and save each output to a list element for dam and dcm 
for (i in 1:length(firstmil_rds)){
  print(filenames[[i]])
  dam_p_value_distribution_1m[[i]]<-find_motif_pval_dam(firstmil_rds[[i]],dam1m)
  
}

for (i in 1:length(firstmil_rds)){
  print(filenames[[i]])
  dcm_p_value_distribution_1m[[i]]<-find_motif_pval_dcm(firstmil_rds[[i]],dcm1m)
  
}

# Create the list to save the plots to
dam_plots <- list()
dcm_plots<- list()


for (i in 1:length(dam_p_value_distribution_1m)){
  df<-dam_p_value_distribution_1m[[i]]%>%
    arrange(Pval)
  dam_plots[[i]]<-ggplot(df,aes(x=Pval,y=(1:length(Pval)/length(Pval)) ))+
    geom_point()+
    geom_vline(xintercept = -4.60, colour= "red", size= 1.5)+
    geom_point(data = damr_test,colour="red")+
    xlim(-60,0)+
    labs(title = paste(filenames[i],"DAM first million sites, unfiltered, window of 1"), x="T test P value", y= "")
}

for (i in 1:length(backup_dcm_p_value_distribution_1m)){
  df<-backup_dcm_p_value_distribution_1m[[i]]%>%
    arrange(Pval)
  dcm_plots[[i]]<-ggplot(df,aes(x=Pval,y=(1:length(Pval)/length(Pval)) ))+
    geom_point()+
    geom_vline(xintercept = -4.78, colour= "red", size= 1.5)+
    geom_point(data = dcmr_test,colour="red")+
    xlim(-60,0)+
    labs(title = paste(filenames[i],"DCM first million sites, unfiltered, window of 1 "), x="T test P value", y= "")
}


all_dam_plots<-plot_grid(plotlist = dam_plots,ncol=1)

all_dam_plots

all_dcm_plots<-plot_grid(plotlist = dcm_plots,ncol=1)

all_dcm_plots

all_plots<-plot_grid(all_dam_plots,all_dcm_plots)
all_plots

ggsave(plot=all_plots,filename = "pvalue_cutoff_tests_H8_6_window1.pdf", width = 18, height = 20, units = "in", path="./"  )

