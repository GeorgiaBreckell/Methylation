
library(tidyverse)
library(cowplot)
library(magick)

#######
#Set the dir to the strain of interest.
#######
setwd("~/Documents/natural_isolates_methylation/")
########

unfiltered_rds_filename <-commandArgs(trailingOnly=T)[1]

motif_sites_filename<- commandArgs(trailingOnly = T)[2]

filter_perc <- as.numeric(commandArgs(trailingOnly = T)[3])



unfiltered_rds <- readRDS(unfiltered_rds_filename)
motif_sites <- read.table(motif_sites_filename)


filenames<-as_vector(unfiltered_rds_filename)%>%
  str_remove("_difference.RDS")%>%
  str_remove("results/nd_merged/\\w{2}/")

motif<-motif_sites_filename%>%
  str_remove("data/motif_locations/")%>%
  str_remove("_\\w{2}_positions.tab")

strain<-as_vector(unfiltered_rds_filename)%>%
  str_remove("/\\w{2}_\\w{1}_\\w{2}_difference.RDS")%>%
  str_remove("results/nd_merged/")

# and a window size to look at to see if there are modified sites
# this window is *exclusive* so a window of 2 will look 1bp up and 
# 1 bp downstream
window <- 2


#Filter the rds files to the first 1million reads 
firstmil_rds<-unfiltered_rds%>%
  mutate(u_test_pval_log=log10(u_test_pval))%>%
  subset(position < 1000000)%>%
  filter(contig != "contig_2_pilon")%>%
  filter(!is.na(u_test_pval))

#To test only motif sites in first 1mil bases
motif_sites_1m<-motif_sites%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%  
  filter(Contig != "contig_2")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)

#Create a list of random sites to create a null for dam and dcm 
motif_random<-motif_sites%>%
  rename("Contig"=V1)%>%
  rename("Site"=V2)%>%  
  filter(Contig != "contig_2")%>%
  mutate(Site=as.numeric(Site))%>%
  filter(Site < 1000000)%>%
  mutate(Site=sample(1:1e6,length(Site)))%>%
  drop_na(Site)%>%
  add_column(Site_mod=0)%>%
  add_column(Pval=0)


#Find p values of random sites 

for (i in 1:length(motif_random$Site)) {
  
  if(dim(subset(firstmil_rds, position < (motif_random[i,2]+ window) & position > (motif_random[i,2]-window)))[1]>0)  {
    # if it does, we switch the bit to indicate it's modified
    motif_random[i,3] <- 1
    
    
    #Then find the P.value of matching sites 
    #record the matching sites and there meta data
    pos.mod<-subset(firstmil_rds, position < (motif_random[i,2]+window) & position > (motif_random[i,2]-window))
    #find the absolute difference from a dam site to a modified site
    pos.abs<-abs(pos.mod$position-motif_random[i,2])
    #find the closest site within the window
    pos.min<-which(pos.abs==min(pos.abs))
    #get the P.values of the closest site/s
    p.values<-pos.mod[pos.min,"u_test_pval_log"]
    #Take the lowest p.value from the closest sites 
    p.value.min<-which(p.values==min(p.values))
    #Save the P.value with the site. 
    motif_random[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
    no_zero<-subset(motif_random,Pval!=0)
  }
}

#Sort these by pvalue after creating for plotting later 

motif_random<-motif_random%>%
  arrange(Pval)


for (i in 1:length(motif_sites_1m$Site)) {
  
  if(dim(subset(firstmil_rds, position < (motif_sites_1m[i,2]+ window) & position > (motif_sites_1m[i,2]-window)))[1]>0)  {
    # if it does, we switch the bit to indicate it's modified
    motif_sites_1m[i,3] <- 1
    
    
    #Then find the P.value of matching sites 
    #record the matching sites and there meta data
    pos.mod<-subset(firstmil_rds, position < (motif_sites_1m[i,2]+window) & position > (motif_sites_1m[i,2]-window))
    #find the absolute difference from a dam site to a modified site
    pos.abs<-abs(pos.mod$position-motif_sites_1m[i,2])
    #find the closest site within the window
    pos.min<-which(pos.abs==min(pos.abs))
    #get the P.values of the closest site/s
    p.values<-pos.mod[pos.min,"u_test_pval_log"]
    #Take the lowest p.value from the closest sites 
    p.value.min<-which(p.values==min(p.values))
    #Save the P.value with the site. 
    motif_sites_1m[i,4] <- pos.mod[pos.min[p.value.min],"u_test_pval_log"]
    no_zero<-subset(motif_sites_1m,Pval!=0)
  }
}


df<-motif_sites_1m%>%
    arrange(Pval)

pval_cutoff<-motif_random[as.integer(length(motif_random[,4])*filter_perc),4]

print(length(motif_random[,4]))
print(filter_perc)
print(pval_cutoff)

motif_plot<-ggplot(df,aes(x=Pval,y=(1:length(Pval)/length(Pval)) ))+
  geom_point()+
  geom_vline(xintercept = pval_cutoff , colour= "red", size= 1.5)+
  geom_point(data = motif_random,colour="red")+
  geom_text(x=-40, y=0.8, label=paste("Pval cutoff= ",pval_cutoff),family="Times",size=8)+
  xlim(-60,0)+
  labs(title = paste(filenames,"DAM first million sites, unfiltered"), x="T test P value", y= "")

motif_plot

#I need to save the plot 

ggsave(plot=motif_plot,filename = paste0("pvalue_cutoff_",filenames,"_",motif,".pdf"), width = 9, height = 10, units = "in", path="./results/nanodisco/cutoff_plots"  )


#filter the rds base on the cutoff 

filtered_rds<-unfiltered_rds%>%
  mutate(u_test_pval_log=log10(u_test_pval))%>%
  filter(u_test_pval_log<pval_cutoff)

write_delim(filtered_rds,paste0("results/nanodisco/filtered_rds/",strain,"/",filenames,"_",motif,"_filtered_rds.tab"),delim = "\t")




