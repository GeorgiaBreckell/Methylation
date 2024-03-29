###NanoDisco Snakemake 

configfile:
    "nano_config.yml"

container: "shub://fanglab/nanodisco"

rule all:
    input:
        expand("/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{sample}_native.sorted.bam", strain=config["strain"],sample=config["sample"]),
        expand("/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{strain}_WGA.sorted.bam",strain=config["strain"]),
        expand("/scratch/georgia/methylation/results/nd_difference/{strain}/{sample}_complete.marker",strain=config["strain"],sample=config["sample"]),
        expand("/scratch/georgia/methylation/results/nd_merged/{strain}/{sample}_difference.RDS",strain=config["strain"],sample=config["sample"])

rule NanoDisco_Preprocess_native:
    input:
        fast5="/scratch/georgia/methylation/data/native_basecalled_fast5/{strain}_{sample}_native_fast5/",
        reference="/scratch/georgia/methylation/data/genomes/{strain}_polished_genome.fasta"
    output:
        "/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{sample}_native.sorted.bam"
    params:
        strain_ID="{sample}_native",
        out_dir="/scratch/georgia/methylation/results/nd_preprocessed/{strain}"
    shell:
        "nanodisco preprocess -p 16 -f {input.fast5} -s {params.strain_ID} -o {params.out_dir} -r {input.reference}"

rule NanoDisco_Preprocess_WGA:
    input:
        fast5="/scratch/georgia/methylation/data/WGA_basecalled_fast5/{strain}",
        reference="/scratch/georgia/methylation/data/genomes/{strain}_polished_genome.fasta"
    output:
        "/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{strain}_WGA.sorted.bam"
    params:
        strain_ID="{strain}_WGA",
        out_dir="/scratch/georgia/methylation/results/nd_preprocessed/{strain}"
    shell:
        "nanodisco preprocess -p 16 -f {input.fast5} -s {params.strain_ID} -o {params.out_dir} -r {input.reference}"

rule NanoDisco_difference: 
    input:
        nd_preprocessed_native="/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{sample}_native.sorted.bam", 
        nd_preprocessed_WGA="/scratch/georgia/methylation/results/nd_preprocessed/{strain}/{strain}_WGA.sorted.bam",   
        reference="/scratch/georgia/methylation/data/genomes/{strain}_polished_genome.fasta"   
    output:
        touch("/scratch/georgia/methylation/results/nd_difference/{strain}/{sample}_complete.marker")
    params:
        nd_preprocessed="/scratch/georgia/methylation/results/nd_preprocessed/{strain}",
        WGA="{strain}_WGA",
        native="{sample}_native",
        out_dir="/scratch/georgia/methylation/results/nd_difference/{strain}/{sample}"
    shell:
        "nanodisco difference -nj 15 -nc 1 -p 16 -i {params.nd_preprocessed} -o {params.out_dir} -w {params.WGA} -n {params.native} -r {input.reference}"

rule NanoDisco_merge:
    input:
        "/scratch/georgia/methylation/results/nd_difference/{strain}/{sample}_complete.marker"
    output:
        "/scratch/georgia/methylation/results/nd_merged/{strain}/{sample}_difference.RDS"
    params:
        sample="{sample}",
        out_dir="/scratch/georgia/methylation/results/nd_merged/{strain}",
        in_dir="/scratch/georgia/methylation/results/nd_difference/{strain}/{sample}"
    shell:
        "nanodisco merge  -d {params.in_dir} -o {params.out_dir} -b {params.sample}"


      


