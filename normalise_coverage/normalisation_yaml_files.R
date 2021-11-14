library(seqinr)
library(yaml)

#Script to make the yaml files for the normalisation snakemake. 
#generates a yaml file with the genome divided into set sized windows, one file per contig. 


F6_genome<-read.fasta(file="data/nanodisco_genomes/F6_polished_genome.fasta")

F6_contig_1<-F6_genome[[1]]

F6_len_1<-length(F6_contig_1)

F6_wins_start_1<-seq(1,F6_len_1,by=5000)

F6_wins_end_1<-seq(5000,F6_len_1,by=5000)

F6_wins_end_1<-append(F6_wins_end_1,F6_len_1)

windows_contig_1<-cbind(F6_wins_start_1,F6_wins_end_1)

windows_contig_1<-as_tibble(windows_contig_1)%>%
  unite(windows_contig_1,sep="-")

F6_contig_4<-F6_genome[[4]]

F6_len_4<-length(F6_contig_4)

F6_wins_start_4<-seq(1,F6_len_4,by=5000)

F6_wins_end_4<-seq(5000,F6_len_4,by=5000)

F6_wins_end_4<-append(F6_wins_end_4,F6_len_4)

windows_contig_4<-cbind(F6_wins_start_4,F6_wins_end_4)

windows_contig_4<-as_tibble(windows_contig_4)%>%
  unite(windows_contig_4,sep="-")


#Every column header becomes a key, and the values become the values for the Yaml, 
#Ideally I would have an extra column with a single value being strain:F6 etc. 
#But if I have a column with one values the others get filled as NA and write_yaml writes NA for each one
#So I've been adding the strain key:value to the top of each file before using it. 
#windows[1,2]<-"F6"

write_yaml(windows_contig_1, "results/F6_windows_contig_1.yaml")
write_yaml(windows_contig_4, "results/F6_windows_contig_4.yaml")


#No issues with the Indenting, and the use of scientific notation when running the snakemake

