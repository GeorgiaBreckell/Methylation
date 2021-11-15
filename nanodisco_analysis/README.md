# Processing Nanodisco output 

Nanodisco outputs a .RDS file for each sample reporting measured and calculated statisitics for both the forward and reverse stand of each site in the genome. As below:
|contig | position  |      dir   |  strand | N_wga |  N_nat  | mean_diff       |t_test_pval   |  u_test_pval   |  u_test_pval_log|
|---|---|---|---|---|---|---|---|---|---|

We use the u_test_pval to determine whether or not a site is modified by first generating a null model of expected pvals for each sample and MTase combination. From this we establish a unique p-value cut off for each sample, representing a 5% false discovery rate which allows us to classify sites as modified or unmodified in a binary way.

The snakefile **Nanodisco_analysis.Snakefile** is used for this analysis. This snakefile requires the .RDS file for each sample, as well as the genomic positions of motifs to be investigated. 

Genomic positions of motifs of interest are identified with the script **Parsing_motifs_in_genomes.py** and described in the [Identify MTase genes section](https://github.com/GeorgiaBreckell/Methylation/tree/master/identify_MTase_genes). 
 
Output from this pipeline is:
- {sample}_{MTase}_filtered_rds.tab   
    in the same format as the original .RDS but with an extra column indicating if a site passed the filter or not.  
- {sample}_{MTase}_site_pvals.tab  
    a .tab file reporting on the pval, WGA and native coverage, and whether a site is modified, for each genomic location of a given MTase's target motifs. 
-  a .pdf plot for visual inspection of the imposed p-value cutoff, showing the null model, observed p-value distributions for each sample/MTase combination and the imposed p-value cutoff. 
