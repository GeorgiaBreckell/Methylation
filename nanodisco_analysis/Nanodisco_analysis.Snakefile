
configfile: 
	"Bin/nanodisco_analysis/methylation_analysis.yml"

rule all:
	input:
		expand("results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_site_pvals.tab", sample=config["sample"],motif=config["motif"], strain=config["strain"]),
		expand("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}.done", sample=config["sample"],motif=config["motif"], strain=config["strain"]),

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
		"Rscript Bin/nanodisco_analysis/filter_rds_all.R {input.unfiltered_rds} {input.sites} {params.filter_level} "

rule genome_wide_mods:
	input:
		mod_sites="results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_site_pvals.tab",
		genome="data/nanodisco_genomes/{strain}_polished_genome.fasta"
	output:
		touch("results/nanodisco/genome_wide_methylation/markers/{strain}/{sample}_{motif}.done")
	shell:
		"python Bin/nanodisco_analysis/Genome_wide_modifications_all.py {input.genome} {input.mod_sites}"

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
