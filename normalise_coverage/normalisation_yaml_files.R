library(seqinr)
library(yaml)

#Script to make the yaml files for the normalisation snakemake. 
#generates a yaml file with the genome divided into set sized windows, one file per contig. 


filename<- commandArgs(trailingOnly=T)[1] #"data/nanodisco_genomes/F6_polished_genome.fasta" 
genome<-read.fasta(file=filename)

strain<-filename%>%
  str_remove("data/nanodisco_genomes/")%>%
  str_remove("_polished_genome.fasta")

for (i in 1:length(genome)){
  contig<-genome[[i]]
  len<-length(contig)
  wins_start<-seq(1,len,by=5000)
  wins_end<-seq(5000,len,by=5000)
  wins_end<-append(wins_end,len)
  windows<-cbind(wins_start,wins_end)
  windows<-as_tibble(windows)%>%
    unite(windows,sep="-")
  write_yaml(windows, paste0("results/",strain,"_windows_contig_",i,".yaml"))
}


#Every column header becomes a key, and the values become the values for the Yaml, 
#Ideally I would have an extra column with a single value being strain:F6 etc. 
#But if I have a column with one values the others get filled as NA and write_yaml writes NA for each one
#So I've been adding the strain key:value to the top of each file before using it. 

#No issues with the Indenting, and the use of scientific notation when running the snakemake

