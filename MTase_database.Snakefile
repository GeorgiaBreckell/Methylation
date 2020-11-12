#I need something that will take as input each assembly, compare it to the Rebase GOLD DB and then transform the output into a "new"
#database to feed into Seqfinder and Abricate. 
##################################################################################
###		Creating Custom MTase database
###         
###     Georgia Breckell  27.07.2020
###
###
###
###################################################################################

configfile:
    "configs/strains.yml"

rule all:
    input:
        expand("results/rebase/{strain}_rebase.paf", strain=config["strain"]),
        expand("results/rebase/{strain}_MTase.pooled", strain=config["strain"]),
	expand("results/abricate/{strain}_abricate.tsv", strain=config["strain"])


#Minimap align Rebase GOLD DB to genome. 
rule Rebase_Minimap:
    input:
        genome="data/genomes/{strain}_polished_genome.fasta"
    output:
        "results/rebase/{strain}_rebase.paf"
    shell:
        "minimap2 {input} rebenz_fixed.fasta -o {output}"

#Samtools select genes that aligned to genome.  
rule pool_MTases:
    input:
        "results/rebase/{strain}_rebase.paf"
    output:
        touch("results/rebase/{strain}_MTase.pooled")
    shell:
        "python MTase_filter.py {input}" 

rule Abricate:
    input: 
        "data/genomes/{strain}_polished_genome.fasta"
    output: 
        touch("results/abricate/{strain}_abricate.tsv")
    params:
        "results/abricate/abricate.tab"
    shell: 
        "abricate --db V_genes {input} >> {params}" 
