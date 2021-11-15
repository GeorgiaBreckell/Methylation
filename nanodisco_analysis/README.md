# Processing Nanodisco output 

Nanodisco outputs a .RDS file for each sample containing reported statistics for both the forward and reverse strand of each site in the genome. 

The snakefile **Nanodisco_analysis.Snakefile** is used to parse this. 

This snakefile requires the .RDS file for each sample, as well as the genomic positions of motifs to be investigated. 

Genomic positions of motifs of interest are identified with the script **Parsing_motifs_in_genomes.py** and described in the [Identify MTase genes section](). 