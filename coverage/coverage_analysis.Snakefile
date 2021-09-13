#Snakemake for the coverage analysis 
configfile: "configs/nano_config.yml"

rule all: 
	input: 
		expand("results/nanodisco/coverage_analysis/{strain}/{sample}_{motif}_native_cov_autocorrelation.pdf", strain=config["strain"],sample=config["sample"],motif=config["motif"])

rule Nanodisco_pvals_and_coverage:
	input:
		rds="results/nd_merged/{strain}/{sample}_difference.RDS",
		motif_sites="data/motif_locations/{motif}_{strain}_positions.tab"

	output:
		"results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_coverage.tab"

	shell:
		"Rscript bin/coverage_fraction_modified.R {input.rds} {input.motif_sites}"


rule Coverage_analysis:
	input:
		pvals_coverage_filename = "results/nanodisco/motif_pvals_with_coverage/{strain}/{sample}_{motif}_coverage.tab",
		MS_filename = "results/nanodisco/genome_wide_methylation/{strain}_{motif}_all_10000/10000_window_{motif}_MS_{sample}_contig_1.csv",
		TS_filename = "results/nanodisco/genome_wide_methylation/{strain}_{motif}_all_10000/10000_window_{motif}_TS_{sample}_contig_1.csv"
  	
	output:
		"results/nanodisco/coverage_analysis/{strain}/{sample}_{motif}_native_cov_autocorrelation.pdf"
	shell:
		"Rscript bin coverage_analysis.R {input.pvals_coverage_filename} {input.MS_filename} {input.TS_filename}"

