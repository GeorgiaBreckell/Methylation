Normalising reads from sequencing experiments for use in Nanodisco 

1. For each strain, a reference genome is required. 

2. Reads are normalised by mapping reads to windows of the genome, a yaml file describing these windows is required.  
This can be produced with the R script "normalisation_yaml_files.R", which will produce a yaml file for each contig in the genome.  
These are then merged into a single file with the following format:
    
    ```
    strain:
        - B2

    sample:
        - B2_2_96       # This will be changed for each sample run. 

    contig_1_windows:
        - 1-5000
        - 5001-10000
    etc

    contig_2_windows:
        - 1-5000
        - 5001-10000
    etc
    ```