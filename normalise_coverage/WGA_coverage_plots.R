library(readr)
SC12A_E1_coverage <- read_table2("Documents/natural_isolates_methylation/results /QC/SC12A_E1_coverage.txt", col_names = FALSE)
SC12A_B2_coverage <- read_table2("Documents/natural_isolates_methylation/results /QC/SC12A_B2_coverage.txt", col_names = FALSE)
SC11A_H8_coverage <- read_table2("Documents/natural_isolates_methylation/results /QC/SC11A_H8_coverage.txt", col_names = FALSE)

SC11A_H8_chromosome <- SC11A_H8_coverage[which(SC11A_H8_coverage[,1]=="contig_1_pilon"),]
SC11A_H8_plasmid <- SC11A_H8_coverage[which(SC11A_H8_coverage[,1]=="contig_2_pilon"),]

SC11A_H8_filtered_chr_cov <- filter(SC11A_H8_chromosome, X2 %% 1000 == 1)
SC11A_H8_filtered_plas_cov <- filter(SC12A_B2_plasmid, X2 %% 100 == 1)

SC12A_E1_chromosome <- SC12A_E1_coverage[which(SC12A_E1_coverage[,1]=="contig_1_pilon"),]
SC12A_E1_plasmid <- SC12A_E1_coverage[which(SC12A_E1_coverage[,1]=="contig_2_pilon"),]

SC12A_E1_filtered_chr_cov <- filter(SC12A_E1_chromosome, X2 %% 1000 == 1)
SC12A_E1_filtered_plas_cov <- filter(SC12A_E1_plasmid, X2 %% 100 == 1)

SC12A_B2_chromosome <- SC12A_B2_coverage[which(SC12A_B2_coverage[,1]=="contig_1_pilon"),]
SC12A_B2_plasmid <- SC12A_B2_coverage[which(SC12A_B2_coverage[,1]=="contig_2_pilon"),]

SC12A_B2_filtered_chr_cov <- filter(SC12A_B2_chromosome, X2 %% 1000 == 1)
SC12A_B2_filtered_plas_cov <- filter(SC12A_B2_plasmid, X2 %% 100 == 1)

SC12A_B2_chromosome_plot<- ggplot(SC12A_B2_filtered_chr_cov, aes(x=X2,y=X3)) + geom_point() + geom_smooth() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw() +labs(title="SC12A_B2 Chromosome Covergae, WGA  ")
SC12A_B2_plasmid_plot<- ggplot(SC12A_B2_filtered_plas_cov, aes(x=X2,y=X3)) + geom_point() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw()+ labs(title="SC12A_B2 Plasmid Covergae, WGA reads ")

SC12A_E1_chromosome_plot<- ggplot(SC12A_E1_filtered_chr_cov, aes(x=X2,y=X3)) + geom_point() + geom_smooth() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw() +labs(title="SC12A_E1 Chromosome Covergae, WGA reads ")
SC12A_E1_plasmid_plot<- ggplot(SC12A_E1_filtered_plas_cov, aes(x=X2,y=X3)) + geom_point() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw()+ labs(title="SC12A_E1 Plasmid Covergae, WGA reads ")

SC11A_H8_chromosome_plot<- ggplot(SC11A_H8_filtered_chr_cov, aes(x=X2,y=X3)) + geom_point() + geom_smooth() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw() +labs(title="SC11A_H8 Chromosome Covergae, WGA reads ")
SC11A_H8_plasmid_plot<- ggplot(SC11A_H8_filtered_plas_cov, aes(x=X2,y=X3)) + geom_point() + ylab(label="Coverage") + xlab(label="Genome position") + theme_bw()+ labs(title="SC11A_H8 Plasmid Covergae, WGA reads")

