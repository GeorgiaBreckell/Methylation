# Creating a snakefile that will link together all of the files/steps that I currently run on the .RDS files 
# after creating them in Nanodisco

#1. Take the raw Nanodisco file and filter it to the appropriate p value cut off. 
# 	Generate a plot to show the cutoff but ideally all filtering will be automatic. 

# Current file for this: 
# fist_mil_filter_pval.R

#These scripts have currently been written to take a Dir of files rather than the single files that snakemake works with 
#Therefore will need to duplicate and tweak the scripts. 

#Input for first mil filter
#	Motif position files 
#	RDS files
# 	Currently makes a random df within the R script but i think the list of random sites would be better as an input. 

#Output 
#	plot of real vs null p values for each motif site
#	Filtered tab file for each RDS file, and motif. 

#Challenges 
#Automate the 5% discovery for the cutoff
# Include the 5% line and value on the plots for record keeping. 
#Change the scripts to run on single files at a time. 

configfile: 
	"methylation_samples.yml"

rule all:
	input:
		expand("results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_filtered_rds.tab", sample=config["sample"],motif=config["motif"], strain=config["strain"])

rule nanodisco_filter_rds
	input:
		sites="data/motif_sites/{strain}/{motif}_positions.tab",
		random_sites = "data/motif_sites/random_sites/{strain}_{motif}.tab",
		unfiltered_rds="results/nanodisco/sample_RDS_files/{sample}_difference.RDS"
	params:
		filter_level= "5"
	output:
		filtered_rds="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_filtered_rds.tab"
	shell:
		"Rscript bin/filter_rds.R {input.unfiltered_rds} {input.sites} {input.random_sites} {params.filter_level} "

