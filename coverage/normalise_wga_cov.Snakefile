
configfile:
    "configs/F6_windows.yaml"

rule all:
    input:
        expand("even_WGA_coverage/{strain}/{sample}_contig_1_fastq_markers/{window}_filtered.marker", strain=config["strain"], window=config["contig_1_windows"],sample=config["sample"])
        expand("even_WGA_coverage/{strain}/{sample}_contig_2_fastq_markers/{window}_filtered.marker", strain=config["strain"], window=config["contig_2_windows"],sample=config["sample"])
        expand("even_WGA_coverage/{strain}/{sample}_contig_3_fastq_markers/{window}_filtered.marker", strain=config["strain"], window=config["contig_3_windows"],sample=config["sample"])
        expand("even_WGA_coverage/{strain}/{sample}_contig_4_fastq_markers/{window}_filtered.marker", strain=config["strain"], window=config["contig_4_windows"],sample=config["sample"])


#Randomly subsample 1000 Nanopore long reads. Generates a new fastq containing the subsampled reads.
rule extract_reads_per_window_1:
    input:
        mapped="even_WGA_coverage/all_reads_mapping/{strain}/{sample}_all_reads_sorted_mapping.bam"
    params:
        "contig_1_pilon:{window}"
    output:
        "even_WGA_coverage/{strain}/{sample}_contig_1_fastqs/{window}_all_reads.fastq"
    shell:
        "bash Bin/filtering_by_window.sh {input.mapped} {params} {output} "

rule filter_fastq_1:
    input:
        "even_WGA_coverage/{strain}/{sample}_contig_1_fastqs/{window}_all_reads.fastq"
    output:
        touch("even_WGA_coverage/{strain}/{sample}_contig_1_fastq_markers/{window}_filtered.marker")
    params:
        "even_WGA_coverage/{strain}/{sample}_contig_1_even_coverage.fastq"
    shell:
        "seqtk sample {input} 78 >> {params}"

#####
rule extract_reads_per_window_2:
    input:
        mapped="even_WGA_coverage/all_reads_mapping/{strain}/{sample}_all_reads_sorted_mapping.bam"
    params:
        "contig_2_pilon:{window}"
    output:
        "even_WGA_coverage/{strain}/{sample}_contig_2_fastqs/{window}_all_reads.fastq"
    shell:
        "bash Bin/filtering_by_window.sh {input.mapped} {params} {output} "

rule filter_fastq_2:
    input:
        "even_WGA_coverage/{strain}/{sample}_contig_2_fastqs/{window}_all_reads.fastq"
    output:
        touch("even_WGA_coverage/{strain}/{sample}_contig_2_fastq_markers/{window}_filtered.marker")
    params:
        "even_WGA_coverage/{strain}/{sample}_contig_2_even_coverage.fastq"
    shell:
        "seqtk sample {input} 89 >> {params}"
#####
rule extract_reads_per_window_3:
    input:
        mapped="even_WGA_coverage/all_reads_mapping/{strain}/{sample}_all_reads_sorted_mapping.bam"
    params:
        "contig_3_pilon:{window}"
    output:
        "even_WGA_coverage/{strain}/{sample}_contig_3_fastqs/{window}_all_reads.fastq"
    shell:
        "bash Bin/filtering_by_window.sh {input.mapped} {params} {output} "


rule filter_fastq_3:
    input:
        "even_WGA_coverage/{strain}/{sample}_contig_3_fastqs/{window}_all_reads.fastq"
    output:
        touch("even_WGA_coverage/{strain}/{sample}_contig_3_fastq_markers/{window}_filtered.marker")
    params:
        "even_WGA_coverage/{strain}/{sample}_contig_3_even_coverage.fastq"
    shell:
        "seqtk sample {input} 91 >> {params}"
#####
rule extract_reads_per_window_4:
    input:
        mapped="even_WGA_coverage/all_reads_mapping/{strain}/{sample}_all_reads_sorted_mapping.bam"
    params:
        "contig_4_pilon:{window}"
    output:
        "even_WGA_coverage/{strain}/{sample}_contig_4_fastqs/{window}_all_reads.fastq"
    shell:
        "bash Bin/filtering_by_window.sh {input.mapped} {params} {output} "

rule filter_fastq_4:
    input:
        "even_WGA_coverage/{strain}/{sample}_contig_4_fastqs/{window}_all_reads.fastq"
    output:
        touch("even_WGA_coverage/{strain}/{sample}_contig_4_fastq_markers/{window}_filtered.marker")
    params:
        "even_WGA_coverage/{strain}/{sample}_contig_4_even_coverage.fastq"
    shell:
        "seqtk sample {input} 160 >> {params}"

