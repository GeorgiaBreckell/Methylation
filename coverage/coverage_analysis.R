library(tidyverse)
library(ggpubr)



####Coverage analysis and plots 

#Read in the coverage and p value data for each sample

pvals_coverage_filename<-"results/nanodisco/motif_pvals_with_coverage/H8/H8_1_M9_DAM_coverage.tab\t"

pvals_cov<-H8_1_M9_DAM_coverage <- read_table2(pvals_coverage_filename)

MS_filename <-"results/nanodisco/genome_wide_methylation/H8_DAM_all_10000/10000_window_DAM_MS_H8_1_M9_contig_1.csv"
  
MS<- read_csv(MS_filename,col_names = FALSE)

TS_filename<-"results/nanodisco/genome_wide_methylation/H8_DAM_all_10000/10000_window_DAM_TS_H8_1_M9_contig_1.csv"

TS<- read_csv(TS_filename,col_names = FALSE)

### Naming Variables ####

sample<-MS_filename%>%
  str_remove("results/nanodisco/genome_wide_methylation/\\w{16}/[0-9]{4,5}_window_\\w{3,5}_MS_")%>%
  str_remove("_contig_\\w{1}.csv")

window<-TS_filename%>%
  str_remove("results/nanodisco/genome_wide_methylation/\\w{16}/")%>%
  str_remove("_\\w{21}_contig_\\w{1}.csv")

strain<-sample%>%
  str_remove("_\\d_\\w{2}")

motif<-MS_filename%>%
  str_remove("results/nanodisco/genome_wide_methylation/\\w{16}/10000_window_")%>%
  str_remove("_MS_\\w{7}_contig_\\w{1}.csv")

fraction_mod<-TS%>%
  cbind(MS$X2)%>%
  rename("Modified_sites" = "MS$X2")%>%
  rename("Window" = "X1")%>%
  rename("No.Sites" = "X2")%>%
  mutate(fraction_modified=Modified_sites/No.Sites)%>%
  mutate(rolling_mean=zoo::rollmedian(fraction_modified,k=5,fill=NA))%>%
  drop_na()



##### Tidy the data ######

#make sure only contig one is present
pvals_cov<-pvals_cov%>% filter(Contig=="contig_1_pilon")%>%
  add_column(windows=0)

#Make the windows which match with the per window genome wide analysis
end_ind<-length(pvals_cov$Site)
end_val<-as_vector(pvals_cov[end_ind,2])
end_buffer<-end_val+10000
windows<-seq(0,end_buffer,by=10000)

#Add the window sizes to the df
for (i in 1:length(windows)){
  w<-which(pvals_cov$Site>windows[i]&pvals_cov$Site<windows[i+1])
  for (j in w){
    pvals_cov[j,6]<-windows[i]
  }
}

#make two new dataframes for WGA and nat coverage (possibly nat coverage is redundant at this point) with the median coverage value for each window
#Wga coverage is all smoothed with a rolling mean in the same way the fraction modified data is.
wga_cov_medians<-pvals_cov%>%
  group_by(windows)%>%
  summarize(wga_cov_mean=median(WGA_cov))%>%
  mutate(rolling_mean_wga=zoo::rollmedian(wga_cov_mean,k=5,fill=NA))%>%
  drop_na()

nat_cov_medians<-pvals_cov%>%
  group_by(windows)%>%
  summarize(nat_cov_mean=mean(Nat_cov))%>%
  mutate(rolling_mean_nat=zoo::rollmedian(nat_cov_mean,k=5,fill=NA))%>%
  drop_na()

cov_frac_windows<-fraction_mod%>%
  cbind(nat_cov_medians[2:3])%>%
  cbind(wga_cov_medians[2:3])

### Plots and metrics ###

#Correlation of Fraction of modified sites vs WGA coverage for each window 

fract_cov_scatter<-ggplot(cov_frac_windows,aes(x=wga_cov_mean,y=rolling_mean))+
  geom_point()+
  labs(x="WGA coverage", y="Fraction modified",title=paste0(sample," ",motif," fraction modified vs median WGA coverage per 10Kbp window"))+ 
  stat_cor(method = "spearman", label.x = .1, label.y = .77)

ggsave(plot=fract_cov_scatter,filename=paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",sample,"_",motif,"_fraction_mod_vs_WGA.pdf"))

#Histogram of WGA coverage for modified sites 

mod_sites_cov<-pvals_cov%>%
  filter(Pval!=0)

WGA_mod_hist<-ggplot(mod_sites_cov,aes(x=WGA_cov))+geom_histogram(binwidth = 10)+
  xlim(-1,1500)+
  labs(title=paste0("Modified sites WGA coverage ",sample))

ggsave(plot=WGA_mod_hist,filename=paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",sample,"_",motif,"_WGA_mod_hist.pdf"))


#Histogram of WGA coverage for unmodified sites 

unmod_sites_cov<-pvals_cov%>%
  filter(Pval==0)

WGA_unmod_hist<-ggplot(unmod_sites_cov,aes(x=WGA_cov))+geom_histogram(binwidth = 10)+
  xlim(-1,1500)+
  labs(title=paste0("Un-Modified sites WGA coverage ",sample))

ggsave(plot=WGA_unmod_hist,filename=paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",sample,"_",motif,"_WGA_unmod_hist.pdf"))


#WGA autocorrelation N1 and N2 scatter plots 


WGA_N0<-wga_cov_medians$wga_cov_mean
WGA_N0<-WGA_N0[1:(length(wga_cov_medians$wga_cov_mean)-2)]

WGA_N<-wga_cov_medians$wga_cov_mean
WGA_N<-WGA_N[1:length(wga_cov_medians$wga_cov_mean)-1]

WGA_N1<-wga_cov_medians$wga_cov_mean
WGA_N1<-WGA_N1[2:length(wga_cov_medians$wga_cov_mean)]

WGA_N2<-wga_cov_medians$wga_cov_mean
WGA_N2<-WGA_N2[3:length(wga_cov_medians$wga_cov_mean)]

WGA1_scatter<-cbind(WGA_N,WGA_N1)
WGA1_scatter<-as_tibble(WGA1_scatter)

WGA2_scatter<-cbind(WGA_N0,WGA_N2)
WGA2_scatter<-as_tibble(WGA2_scatter)

WGA1_scatter_plot<-ggplot(WGA1_scatter,aes(x=WGA_N,y=WGA_N1))+geom_point()+
  stat_cor(method = "spearman", label.x = 10, label.y = 120)+
  labs(title=paste0(sample," ",motif," WGA Coverage Autocorrelation N vs N+1"))

WGA2_scatter_plot<-ggplot(WGA2_scatter,aes(x=WGA_N0,y=WGA_N2))+geom_point()+
  stat_cor(method = "spearman", label.x = 10, label.y = 120)+
  labs(title=paste0(sample," ",motif," WGA Coverage Autocorrelation N vs N+2"))

WGA_auto_cov<-cowplot::plot_grid(WGA1_scatter_plot,WGA2_scatter_plot,ncol=2)

ggsave(plot=WGA_auto_cov,filename=paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",sample,"_",motif,"_WGA_cov_autocorrelation.pdf"))



#Native Autocorrelation N1 and N2 scatter plots 

NAT_N0<-nat_cov_medians$nat_cov_mean
NAT_N0<-NAT_N0[1:(length(nat_cov_medians$nat_cov_mean)-2)]

NAT_N<-nat_cov_medians$nat_cov_mean
NAT_N<-NAT_N[1:length(nat_cov_medians$nat_cov_mean)-1]

NAT_N1<-nat_cov_medians$nat_cov_mean
NAT_N1<-NAT_N1[2:length(nat_cov_medians$nat_cov_mean)]

NAT_N2<-nat_cov_medians$nat_cov_mean
NAT_N2<-NAT_N2[3:length(nat_cov_medians$nat_cov_mean)]

N1_scatter<-cbind(NAT_N,NAT_N1)
N1_scatter<-as_tibble(N1_scatter)

N2_scatter<-cbind(NAT_N0,NAT_N2)
N2_scatter<-as_tibble(N2_scatter)

nat1_scatter_plot<-ggplot(N1_scatter,aes(x=NAT_N,y=NAT_N1))+geom_point()+
  stat_cor(method = "spearman", label.x = 10, label.y = 120)+
  labs(title=paste0(sample," ",motif," Native Coverage Autocorrelation N vs N+1"))

nat2_scatter_plot<-ggplot(N2_scatter,aes(x=NAT_N0,y=NAT_N2))+geom_point()+
  stat_cor(method = "spearman", label.x = 10, label.y = 120)+
  labs(title=paste0(sample," ",motif," Native Coverage Autocorrelation N vs N+2"))

Nat_auto_cov<-cowplot::plot_grid(nat1_scatter_plot,nat2_scatter_plot)

ggsave(plot=Nat_auto_cov,filename=paste0("results/nanodisco/motif_pvals_with_coverage/",strain,"/",sample,"_",motif,"_native_cov_autocorrelation.pdf"))





