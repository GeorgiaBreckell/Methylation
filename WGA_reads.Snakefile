#I need something that will take as input each assembly, compare it to the Rebase GOLD DB and then transform the output into a "new"
#database to feed into Seqfinder and Abricate. 
##################################################################################
###		QC for WGA reads, 
###		Align reads to genome (done)
###		Filter reads that map (done)
###     Plot coverage of the read mapping (done) 
###  	Run Tombo "de novo" to check for methylation in reads (done)
###  	Run Tombo WGA model with raw reads and Nanodisco
###     
###		Georgia Breckell  8.08.2020
###		!!!!GET COVERAGE PLOTTING SCRIPT FROM ASSEMBLIES BIN!!!!!
###
###
###################################################################################

configfile:
    "strains.yml"

rule all:
    input:
        #expand("data/native_reads/{strain}/fast5_moved.marker", strain=config["strain"]),
        #expand("data/WGA_reads/{strain}/fast5_moved.marker", strain=config["strain"]),
        expand("data/native_raw_fast5/{strain}_basecalled", strain=config["strain"]),
        expand("data/WGA_raw_fast5/{strain}_basecalled", strain=config["strain"]),
        #expand("results/WGA/QC/{strain}_mapped_read_stats.txt", strain=config["strain"]),
        #expand("results/WGA/QC/{strain}_unmapped_read_stats.txt", strain=config["strain"]),
        #expand("results/WGA/QC/{strain}_coverage.pdf",strain=config["strain"]),
        #expand("results/WGA/Tombo/{strain}_de_novo_testing.fasta", strain=config["strain"]),
        #expand("results/WGA/Tombo/{strain}_native_compare_detection.tombo.stats", strain=config["strain"])

#Guppy basecalling by strain
rule guppy_basecaller_WGA:
	input:
		"data/WGA_raw_fast5/{strain}/"
	output:
		touch("data/WGA_raw_fast5/{strain}_basecalled")
	params:
		"data/WGA_reads/{strain}",
		workspace="data/WGA_reads/{strain}/workspace/"
	shell:
		"guppy_basecaller -i {input}  -s {params}  --flowcell FLO-MIN106 --kit SQK-RBK004 --fast5_out -x \"cuda:0\" "

rule guppy_basecaller_native:
	input:
		"data/native_raw_fast5/{strain}/reads/"
	output:
		touch("data/native_raw_fast5/{strain}_basecalled")
	params:
		"data/native_reads/{strain}",
		workspace="data/native_reads/{strain}/workspace/"
	shell:
		"guppy_basecaller -i {input}  -s {params}  --flowcell FLO-MIN106 --kit SQK-RBK004 --fast5_out -x \"cuda:0\" "

rule move_basecalled_fast5:
	input:
		WGA="data/WGA_reads/{strain}/workspace/",
		native="data/native_reads/{strain}/workspace/",
	output:
		WGA=touch("data/WGA_reads/{strain}/fast5_moved.marker"),
		native=touch("data/native_reads/{strain}/fast5_moved.marker")
	params:
		WGA="data/WGA_basecalled_fast5/{strain}",
		native="data/native_basecalled_fast5/{strain}"
	run:
		shell("mv {input.WGA} {params.WGA}"),
		shell("mv {input.native} {params.native}")

#Minimap align Rebase GOLD DB to genome. 
rule Align_WGA_to_genome:
    input:
        genome="data/genomes/{strain}_polished_genome.fasta",
        reads="data/WGA_reads/{strain}_wga.fastq.gz"
    output:
        "results/WGA/QC/{strain}_wga_aligned.sam"
    shell:
        "minimap2 -a --secondary=no {input.genome} {input.reads} -o {output}"

#Samtools select genes that aligned to genome.  
rule sort_samfile:
    input:
        "results/WGA/QC/{strain}_wga_aligned.sam"
    output:
        "results/WGA/QC/{strain}_sorted.bam" 
    shell:
        "samtools sort {input} > {output}" 

rule isolate_mapped_WGA_reads:
	input: 
		"results/WGA/QC/{strain}_sorted.bam"
	output: 
		"results/WGA/QC/{strain}_mapped.bam"
	shell:
		"samtools view -F 4 -b -o {output} {input}"

rule isolate_unmapped_WGA_reads:
	input: 
		"results/WGA/QC/{strain}_sorted.bam"
	output: 
		"results/WGA/QC/{strain}_unmapped.bam"
	shell:
		"samtools view -f 4 -b -o {output} {input}"

rule mapped_reads_fasta:
	input: 
		"results/WGA/QC/{strain}_mapped.bam"
	output: 
		"results/WGA/{strain}_mapped_reads.fasta"
	shell: 
		"samtools fasta {input} > {output}"

rule unmapped_reads_fasta:
	input: 
		"results/WGA/QC/{strain}_unmapped.bam"
	output: 
		"results/WGA/{strain}_unmapped_reads.fasta"
	shell: 
		"samtools fasta {input} > {output}"

rule mapped_reads_stats: 
	input:
		"results/WGA/{strain}_unmapped_reads.fasta"
	output: 
		"results/WGA/QC/{strain}_unmapped_read_stats.txt"
	shell:
		"seqkit stats -aT {input} > {output}"

rule unmapped_reads_stats: 
	input:
		"results/WGA/{strain}_mapped_reads.fasta"
	output: 
		"results/WGA/QC/{strain}_mapped_read_stats.txt"
	shell:
		"seqkit stats -aT {input} > {output}"

rule aligned_read_coverage:
	input: 
		"results/WGA/QC/{strain}_sorted.bam"
	output:
		"results/WGA/QC/{strain}_coverage.txt"
	shell:
		"samtools depth {input} > {output}"

rule Plotting_coverage: 
	input: 
		"results/WGA/QC/{strain}_coverage.txt"
	output: 
		"results/WGA/QC/{strain}_coverage.pdf"
	shell:
		"Rscript ./scripts/coverage_plots_g.R {input}"		

### Using tombo de novo model to check the WGA reads do not have methylation, ie the experiement processed correctly. 

rule resquiggle:
	input:
		genome="data/genomes/{strain}_polished_genome.fasta",
		wga="data/WGA_fast5/{strain}/" 
	output:
		touch("results/WGA/Tombo/{strain}_resquiggle_complete.marker")
	benchmark:
		"benchmark/{strain}/resquiggle_benchmark.txt"
	log:
		"logs/{strain}/resquiggle.log"
	shell:
		"tombo resquiggle {input.wga} {input.genome} --processes 4 --num-most-common-errors 5"
		" 2> {log}"

rule tombo_detect_methylation_denovoQC:
	input:
		wga="data/WGA_fast5/{strain}/",
		resquiggled="results/WGA/Tombo/{strain}_resquiggle_complete.marker"
	output:
		"results/WGA/Tombo/{strain}_denovo_detection.tombo.stats"
	benchmark:
		"benchmark/{strain}/tombo_detect_methylation_benchmark.txt"
	log:
		"logs/{strain}/tombo_detect_methylation.log"
	params:
		statistics_basename="results/WGA/Tombo/{strain}_denovo_detection"
	shell:
		"tombo detect_modifications de_novo --fast5-basedirs {input.wga} --statistics-file-basename {params.statistics_basename} --processes 4"

rule significantly_modifed_sites_denovoQC:
	input: 
		wga="data/WGA_fast5/{strain}/",
		stats_filename="results/WGA/Tombo/{strain}_denovo_detection.tombo.stats"
	output:
		sites="results/WGA/Tombo/{strain}_de_novo_testing.fasta"
	benchmark:
		"benchmark/{strain}/tombo_significantly_modified_benchmark.txt"
	log:
		"logs/{strain}/significantly_modified.log"
	shell:
		"tombo text_output signif_sequence_context --fast5-basedirs {input.wga} --statistics-filename {input.stats_filename} --sequences-filename {output.sites} --num-regions 5000" 

rule tombo_compare_model:
	input:
		wga="data/WGA_fast5/{strain}/",
		resquiggled="results/WGA/Tombo/{strain}_resquiggle_complete.txt",
		native="data/native_fast5/{strain}/"
	output:
		"results/WGA/Tombo/{strain}_native_compare_detection.tombo.stats"
	benchmark:
		"benchmark/{strain}/tombo_compare_methylation_benchmark.txt"
	log:
		"logs/{strain}/tombo_compare_methylation.log"
	params:
		statistics_basename="results/WGA/Tombo/{strain}_native_compare"
	shell:
		"tombo detect_modifications model_sample_compare --fast5-basedirs {input.wga} --minimum-test-reads 30 --control-fast5-basedirs {input.native} --statistics-file-basename {params.statistics_basename} --processes 4"







###Nanodisco 

rule Nanodisco_preprocess_WGA:
	input:
		wga="data/fast5/{strain}/reads/",
		resquiggled="results/WGA/Tombo/{strain}_resquiggle_complete.txt"
	output:
		"results/WGA/Tombo/{strain}_denovo_detection.tombo.stats"
	benchmark:
		"benchmark/{strain}/tombo_detect_methylation_benchmark.txt"
	log:
		"logs/{strain}/Nanodisco_preprocess_WGA.log"
	params:
		ref_genome="results/WGA/Tombo/{strain}_denovo_detection"
	shell:
		"nanodisco preprocess -p 4 -f {input.wga} -s {strain}_WGA -o {output.dir} -r {params.ref_genome}"

