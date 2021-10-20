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

#2. 
#After filtering the reads I output a tab file that contains only the genome sites that meet the p.value cut off 

# I then use that file to determine the fraction of modified sites for each methyltransferase 

# All of this is currently only done for DAM and DCM 

#3. I then find the fraction of modified sites in each genomic region and plot this. 

#4. Un-complete step is to overlay the fraction of modified sites with an anotated genome. 

#Ideally run on agnes, therefore, dir structure is results/nd_merged/{strain}_{sample}_difference.rds

configfile: 
	"bin/nanodisco_analysis/methylation_analysis.yml"

rule all:
	input:
		#expand("results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_fwd_site_pvals.tab", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		#expand("results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_rev_site_pvals.tab", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		expand("results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_site_pvals.tab", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		#expand("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_fwd.done", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		#expand("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_rev.done", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		expand("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}.done", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		#expand("results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}_fwd.pdf", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		#expand("results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}_rev.pdf", sample=config["sample"],motif=config["motif"], strain=config["strain"])
		#expand("results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}.pdf", sample=config["sample"],motif=config["motif"], strain=config["strain"])

rule nanodisco_filter_rds:
	input:
		sites="data/motif_locations/{motif}_{strain}_positions.tab",
		unfiltered_rds="results/nd_merged/{strain}/{sample}_difference.RDS"
	params:
		filter_level= "0.05"
	output:
		filtered_rds_fwd="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_fwd_filtered_rds.tab",
		filtered_rds_rev="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_rev_filtered_rds.tab"
	shell:
		"Rscript bin/nanodisco_analysis/filter_rds.R {input.unfiltered_rds} {input.sites} {params.filter_level} "

rule nanodisco_filter_rds_global:
	input:
		sites="data/motif_locations/{motif}_{strain}_positions.tab",
		unfiltered_rds="results/nd_merged/{strain}/{sample}_difference.RDS"
	params:
		filter_level= "0.05"
	output:
		filtered_rds="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_filtered_rds.tab",
		motif_pvals_with_coverage="results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_site_pvals.tab"
	shell:
		"Rscript bin/nanodisco_analysis/filter_rds_all.R {input.unfiltered_rds} {input.sites} {params.filter_level} "

rule all_pvals_forward: 
	input:
		sites="data/motif_locations/{motif}_{strain}_positions.tab",
		filtered_rds="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_fwd_filtered_rds.tab"
	output:
		"results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_fwd_site_pvals.tab"
	shell:
		"Rscript bin/nanodisco_analysis/filtered_motif_site_pvals.R  {input.filtered_rds} {input.sites}"

rule all_pvals_reverse: 
	input:
		sites="data/motif_locations/{motif}_{strain}_positions.tab",
		filtered_rds="results/nanodisco/filtered_rds/{strain}/{sample}_{motif}_rev_filtered_rds.tab"
	output:
		"results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_rev_site_pvals.tab"
	shell:
		"Rscript bin/nanodisco_analysis/filtered_motif_site_pvals.R  {input.filtered_rds} {input.sites}"

rule genome_wide_mods:
	input:
		mod_sites="results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_site_pvals.tab",
		genome="data/nanodisco_genomes/{strain}_polished_genome.fasta"
	output:
		touch("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}.done")
	shell:
		"python bin/nanodisco_analysis/Genome_wide_modifications_all.py {input.genome} {input.mod_sites}"

rule genome_wide_mods_forward:
	input:
		mod_sites="results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_fwd_site_pvals.tab",
		genome="data/nanodisco_genomes/{strain}_polished_genome.fasta"
	output:
		touch("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_fwd.done")
	shell:
		"python bin/nanodisco_analysis/Genome_wide_modifications.py {input.genome} {input.mod_sites}"

rule genome_wide_mods_reverse:
	input:
		mod_sites="results/nanodisco/motif_pvals/{strain}/{sample}_{motif}_rev_site_pvals.tab",
		genome="data/nanodisco_genomes/{strain}_polished_genome.fasta"
	output:
		touch("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_rev.done")
	shell:
		"python bin/nanodisco_analysis/Genome_wide_modifications.py {input.genome} {input.mod_sites}"

rule genome_wide_plots:
	input:
		"results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}.done"
	output:
		main="results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}.pdf",
		ACF="results/nanodisco/methylation_ACF_plots/{strain}_10000/{sample}_{motif}.pdf"
	params:
		mod_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_MS_{sample}_contig_1.csv",
		all_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_TS_{sample}_contig_1.csv"
	shell:
		"Rscript bin/nanodisco_analysis/methylation_pattern_plots_all.R {params.mod_sites} {params.all_sites}"


rule genome_wide_plots_forward:
	input:
		"results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_fwd.done"
	output:
		main="results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}_fwd.pdf",
		ACF="results/nanodisco/methylation_ACF_plots/{strain}_10000/{sample}_{motif}_fwd.pdf"
	params:
		mod_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_fwd_MS_{sample}_contig_1.csv",
		all_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_fwd_TS_{sample}_contig_1.csv"
	shell:
		"Rscript bin/nanodisco_analysis/methylation_pattern_plots.R {params.mod_sites} {params.all_sites}"

rule genome_wide_plots_reverse:
	input:
		"results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}_rev.done"
	output:
		main="results/nanodisco/methylation_pattern_plots/{strain}_10000/{sample}_{motif}_rev.pdf",
		ACF="results/nanodisco/methylation_ACF_plots/{strain}_10000/{sample}_{motif}_rev.pdf"
	params:
		mod_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_rev_MS_{sample}_contig_1.csv",
		all_sites="results/nanodisco/genome_wide_methylation/10000_window_{motif}_rev_TS_{sample}_contig_1.csv"
	shell:
		"Rscript bin/nanodisco_analysis/methylation_pattern_plots.R {params.mod_sites} {params.all_sites}"

