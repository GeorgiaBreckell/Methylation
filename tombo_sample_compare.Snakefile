configfile:
    "configs/methyl_config.yml"
    
rule all:
	input: expand("results/tombo_compare/{strain}/modified_motifs_meme/meme.html", strain=config["strain"]), #expand("{strain}/fast5/fast5_reads_grouped.txt", strain=FAST5),

rule resquiggle_WGA:
	input:
		genome=("data/genomes/{strain}_polished_genome.fasta"),
		fast5=("data/WGA_basecalled_fast5/{strain}/") #Ideally this should read ("{strain}/fast5/reads")
	output:
		touch("results/tombo_compare/{strain}/WGA_resquiggle_complete.txt")
	benchmark:
		"results/tombo_compare/{strain}/WGA_resquiggle_benchmark.txt"
	log:
		"results/tombo_compare/{strain}/WGA_resquiggle.log"
	shell:
		"tombo resquiggle {input.fast5} {input.genome} --processes 4 --num-most-common-errors 5 --overwrite"
		" 2> {log}"
		
rule resquiggle_native:
	input:
		genome=("data/genomes/{strain}_polished_genome.fasta"),
		fast5=("data/native_basecalled_fast5/{strain}/") #Ideally this should read ("{strain}/fast5/reads")
	output:
		touch("results/tombo_compare/{strain}/native_resquiggle_complete.txt")
	benchmark:
		"results/tombo_compare/{strain}/benchmarks/native_resquiggle_benchmark.txt"
	log:
		"results/tombo_compare/{strain}/benchmarks/native_resquiggle.log"
	shell:
		"tombo resquiggle {input.fast5} {input.genome} --processes 4 --num-most-common-errors 5 --overwrite"
		" 2> {log}"

rule tombo_detect_methylation:
	input:
		native=("data/native_basecalled_fast5/{strain}/"),
		native_resquiggled=("results/tombo_compare/{strain}/native_resquiggle_complete.txt"),
		WGA=("data/WGA_basecalled_fast5/{strain}/"),
		WGA_resquiggled=("results/tombo_compare/{strain}/WGA_resquiggle_complete.txt")
	output:
		("results/tombo_compare/{strain}/sample_compare.tombo.stats")
	benchmark:
		("results/tombo_compare/{strain}/benchmarks/tombo_detect_methylation_benchmark.txt")
	log:
		("results/tombo_compare/{strain}/benchmarks/tombo_detect_methylation.log")
	params:
		statistics_basename=("results/tombo_compare/{strain}/sample_compare")
	shell:
		"tombo detect_modifications model_sample_compare --fast5-basedirs {input.native} --control-fast5-basedirs {input.WGA} --statistics-file-basename {params.statistics_basename} --processes 4"
		

rule significantly_modifed_sites:
	input: 
		fast5=("data/native_basecalled_fast5/{strain}/"),
		stats_filename=("results/tombo_compare/{strain}/sample_compare.tombo.stats")
	output:
		sites=("results/tombo_compare/{strain}/sample_compare.fasta")
	benchmark:
		("results/tombo_compare/{strain}/benchmarks/tombo_significantly_modified_benchmark.txt")
	log:
		("results/tombo_compare/{strain}/benchmarks/significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 10000 --num-bases 50"
		

rule motif_finder:
	input:
		sites=("results/tombo_compare/{strain}/sample_compare.fasta")
	output:
		file=("results/tombo_compare/{strain}/modified_motifs_meme/meme.html")
	params:
		motifs=directory("results/tombo_compare/{strain}/modified_motifs_meme")
	benchmark:
		("results/tombo_compare/{strain}/benchmarks/motif_finder_benchmark.txt")
	log:
		("results/tombo_compare/{strain}/benchmarks/modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 -maxsize 1500000 "
		##max size set very high due to increased number of sites and increased bp output for each site. Maxsize is measured in characters. 
