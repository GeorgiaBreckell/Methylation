configfile:
    "methyl_config.yml"

#FAST5, =glob_wildcards("{strain}/fast5/")
#####ISSUE, for some reason the wildcard {strain} is being saved as SC12A_*/fast5, therefore 
##### rather than defining the reads input as {strain}/fast5/reads, if i just define it as 
##### {strain}/reads it reads this as the full path of "{strain}/fast5/reads" I dont know why
##### and this shouldnt be happening. 

rule all:
	input: 
		expand("{strain}/dcm_modified_motifs_meme/meme.html", strain=config["strain"]),
		expand("{strain}/dam_modified_motifs_meme/meme.html", strain=config["strain"]),
		expand("{strain}/6ma_modified_motifs_meme/meme.html", strain=config["strain"]),
		expand("{strain}/5mc_modified_motifs_meme/meme.html", strain=config["strain"]),
	 #expand("{strain}/fast5/fast5_reads_grouped.txt", strain=FAST5),

rule resquiggle:
	input:
		genome=("reference_genomes/{strain}_polished_genome.fasta"),
		fast5=("{strain}/fast5/{strain}_fast5_reads/") #Ideally this should read ("{strain}/fast5/reads")
	output:
		touch("{strain}/fast5/resquiggle_complete.txt")
	benchmark:
		"{strain}/resquiggle_benchmark.txt"
	log:
		"{strain}/resquiggle.log"
	shell:
		"tombo resquiggle {input.fast5} {input.genome} --processes 4 --num-most-common-errors 5"
		" 2> {log}"
		

rule tombo_detect_methylation_5mc:
	input:
		fast5=("{strain}/fast5/{strain}_fast5_reads"),
		resquiggled=("{strain}/fast5/resquiggle_complete.txt")
	output:
		("{strain}/tombo/5mc_detection.5mC.tombo.stats")
	benchmark:
		("{strain}/benchmarks/tombo_5mc_detect_methylation_benchmark.txt")
	log:
		("{strain}/benchmarks/tombo_5mc_detect_methylation.log")
	params:
		statistics_basename=("{strain}/tombo/5mc_detection"),
		per_read_statistics_basename=("{strain}/tombo/5mc_per_read_detection")
	shell:
		"tombo detect_modifications alternative_model --alternate-bases 5mC --fast5-basedirs {input.fast5} --statistics-file-basename {params.statistics_basename} --per-read-statistics-basename {params.per_read_statistics_basename} --processes 4"
		
rule tombo_detect_methylation_dam:
	input:
		fast5=("{strain}/fast5/{strain}_fast5_reads"),
		resquiggled=("{strain}/fast5/resquiggle_complete.txt")
	output:
		("{strain}/tombo/dam_detection.dam.tombo.stats")
	benchmark:
		("{strain}/benchmarks/tombo_dam_detect_methylation_benchmark.txt")
	log:
		("{strain}/benchmarks/tombo_dam_detect_methylation.log")
	params:
		statistics_basename=("{strain}/tombo/dam_detection"),
		per_read_statistics_basename=("{strain}/tombo/dam_per_read_detection")
	shell:
		"tombo detect_modifications alternative_model --alternate-bases dam --fast5-basedirs {input.fast5} --statistics-file-basename {params.statistics_basename} --per-read-statistics-basename {params.per_read_statistics_basename} --processes 4"
		
rule tombo_detect_methylation_dcm:
	input:
		fast5=("{strain}/fast5/{strain}_fast5_reads"),
		resquiggled=("{strain}/fast5/resquiggle_complete.txt")
	output:
		("{strain}/tombo/dcm_detection.dcm.tombo.stats")
	benchmark:
		("{strain}/benchmarks/tombo_dcm_detect_methylation_benchmark.txt")
	log:
		("{strain}/benchmarks/tombo_dcm_detect_methylation.log")
	params:
		statistics_basename=("{strain}/tombo/dcm_detection"),
		per_read_statistics_basename=("{strain}/tombo/dcm_per_read_detection")
	shell:
		"tombo detect_modifications alternative_model --alternate-bases dcm --fast5-basedirs {input.fast5} --statistics-file-basename {params.statistics_basename} --per-read-statistics-basename {params.per_read_statistics_basename} --processes 4"
		
rule tombo_detect_methylation_6ma:
	input:
		fast5=("{strain}/fast5/{strain}_fast5_reads"),
		resquiggled=("{strain}/fast5/resquiggle_complete.txt")
	output:
		("{strain}/tombo/6ma_detection.6mA.tombo.stats")
	benchmark:
		("{strain}/benchmarks/tombo_6ma_detect_methylation_benchmark.txt")
	log:
		("{strain}/benchmarks/tombo_6ma_detect_methylation.log")
	params:
		statistics_basename=("{strain}/tombo/6ma_detection"),
		per_read_statistics_basename=("{strain}/tombo/6ma_per_read_detection")
	shell:
		"tombo detect_modifications alternative_model --alternate-bases 6mA --fast5-basedirs {input.fast5} --statistics-file-basename {params.statistics_basename} --per-read-statistics-basename {params.per_read_statistics_basename} --processes 4"


rule significantly_modifed_sites_5mc:
	input: 
		fast5=("{strain}/fast5/{strain}_fast5_reads/"),
		stats_filename=("{strain}/tombo/5mc_detection.5mC.tombo.stats")
	output:
		sites=("{strain}/tombo/5mc_testing.fasta")
	benchmark:
		("{strain}/benchmarks/tombo_5mc_significantly_modified_benchmark.txt")
	log:
		("{strain}/benchmarks/5mc_significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000"
		

rule motif_finder_5mc:
	input:
		sites=("{strain}/tombo/5mc_testing.fasta")
	output:
		file=("{strain}/5mc_modified_motifs_meme/meme.html")
	params:
		motifs=directory("{strain}/5mc_modified_motifs_meme")
	benchmark:
		("{strain}/benchmarks/5mc_motif_finder_benchmark.txt")
	log:
		("{strain}/benchmarks/5mc_modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 "

rule significantly_modifed_sites_6ma:
	input: 
		fast5=("{strain}/fast5/{strain}_fast5_reads/"),
		stats_filename=("{strain}/tombo/6ma_detection.6mA.tombo.stats")
	output:
		sites=("{strain}/tombo/6ma_testing.fasta")
	benchmark:
		("{strain}/benchmarks/tombo_6ma_significantly_modified_benchmark.txt")
	log:
		("{strain}/benchmarks/6ma_significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000"
		

rule motif_finder_6ma:
	input:
		sites=("{strain}/tombo/6ma_testing.fasta")
	output:
		file=("{strain}/6ma_modified_motifs_meme/meme.html")
	params:
		motifs=directory("{strain}/6ma_modified_motifs_meme")
	benchmark:
		("{strain}/benchmarks/6ma_motif_finder_benchmark.txt")
	log:
		("{strain}/benchmarks/6ma_modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 "

rule significantly_modifed_sites_dam:
	input: 
		fast5=("{strain}/fast5/{strain}_fast5_reads/"),
		stats_filename=("{strain}/tombo/dam_detection.dam.tombo.stats")
	output:
		sites=("{strain}/tombo/dam_testing.fasta")
	benchmark:
		("{strain}/benchmarks/tombo_dam_significantly_modified_benchmark.txt")
	log:
		("{strain}/benchmarks/dam_significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000"
		

rule motif_finder_dam:
	input:
		sites=("{strain}/tombo/dam_testing.fasta")
	output:
		file=("{strain}/dam_modified_motifs_meme/meme.html")
	params:
		motifs=directory("{strain}/dam_modified_motifs_meme")
	benchmark:
		("{strain}/benchmarks/dam_motif_finder_benchmark.txt")
	log:
		("{strain}/benchmarks/dam_modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 "

rule significantly_modifed_sites_dcm:
	input: 
		fast5=("{strain}/fast5/{strain}_fast5_reads/"),
		stats_filename=("{strain}/tombo/dcm_detection.dcm.tombo.stats")
	output:
		sites=("{strain}/tombo/dcm_testing.fasta")
	benchmark:
		("{strain}/benchmarks/tombo_dcm_significantly_modified_benchmark.txt")
	log:
		("{strain}/benchmarks/dcm_significantly_modified.log")
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.fast5} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000"
		

rule motif_finder_dcm:
	input:
		sites=("{strain}/tombo/dcm_testing.fasta")
	output:
		file=("{strain}/dcm_modified_motifs_meme/meme.html")
	params:
		motifs=directory("{strain}/dcm_modified_motifs_meme")
	benchmark:
		("{strain}/benchmarks/dcm_motif_finder_benchmark.txt")
	log:
		("{strain}/benchmarks/dcm_modified_motifs.log")
	shell:
		"meme -oc {params.motifs} -dna -mod zoops {input.sites} -nmotifs 5 "
		




