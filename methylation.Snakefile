configfile:
    "methyl_config.yml"

#FAST5, =glob_wildcards("{strain}/fast5/")
#####ISSUE, for some reason the wildcard {strain} is being saved as SC12A_*/fast5, therefore 
##### rather than defining the reads input as {strain}/fast5/reads, if i just define it as 
##### {strain}/reads it reads this as the full path of "{strain}/fast5/reads" I dont know why
##### and this shouldnt be happening. 

rule all:
	input: expand("{strain}/modified_motifs_meme/meme.html", strain=config["strain"]), #expand("{strain}/fast5/fast5_reads_grouped.txt", strain=FAST5),

rule resquiggle:
	input:
		genome=("reference_genomes/{strain}_polished_genome.fasta"),
		fast5=("{strain}/fast5/reads/") #Ideally this should read ("{strain}/fast5/reads")
	output:
		touch("{strain}/fast5/resquiggle_complete.txt")
	benchmark:
		"{strain}/resquiggle_benchmark.txt"
	log:
		"{strain}/resquiggle.log"
	shell:
		"tombo resquiggle {input.fast5} {input.genome} --processes 4 --num-most-common-errors 5 --overwrite"
		" 2> {log}"
		

rule tombo_detect_methylation:
	input:
		fast5=("{strain}/fast5/reads/"),
		resquiggled=("{strain}/fast5/resquiggle_complete.txt")
	output:
		("{strain}/tombo/denovo_detection.tombo.stats")
	benchmark:
		("{strain}/benchmarks/tombo_detect_methylation_benchmark.txt")
	log:
		("{strain}/benchmarks/tombo_detect_methylation.log")
	params:
		statistics_basename=("{strain}/tombo/denovo_detection")
	shell:
		"tombo detect_modifications de_novo --fast5-basedirs {input.fast5} --statistics-file-basename {params.statistics_basename} --processes 4"
		

rule significantly_modifed_sites:
	input: 
		fast5=("{strain}/fast5/reads/"),
		stats_filename=("{strain}/tombo/denovo_detection.tombo.stats")
	output:
		sites=("{strain}/tombo/de_novo_testing.fasta")
	benchmark:
		("{strain}/benchmarks/tombo_significantly_modified_benchmark.txt")
	log:
		("{strain}/benchmarks/significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000"
		

rule motif_finder:
	input:
		sites=("{strain}/tombo/de_novo_testing.fasta")
	output:
		file=("{strain}/modified_motifs_meme/meme.html")
	params:
		motifs=directory("{strain}/modified_motifs_meme")
	benchmark:
		("{strain}/benchmarks/motif_finder_benchmark.txt")
	log:
		("{strain}/benchmarks/modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 "
		




