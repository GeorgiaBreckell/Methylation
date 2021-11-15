## Normalising reads from sequencing experiments for use in Nanodisco 

1. For each strain, a reference genome is required. 

2. Reads are normalised by mapping reads to windows of the genome, a yaml file describing these windows is required.  
This can be produced with the R script *normalisation_yaml_files.R*, which will produce a yaml file for each contig in the genome.  
```
Rscript normalisation_yaml_files.R {reference_genome.fasta}
```
Before use update the Rscript with the reference genome file naming structure, the output path, and the taret window size. 

3.  Merged each contigs yaml file into a single file with the following format:
    
    ```
    strain:
        - B2

    sample:
        - B2_2_96       # This will be changed for each sample run. 

    contig_1_windows:
        - 1-5000
        - 5001-10000
    etc

    contig_2_windows:
        - 1-5000
        - 5001-10000
    etc
    ```
4. Map all reads onto the genome with no secondary mapping, and extract mapped reads into a fastq file
```
minimap2 -x map-ont --secondary=no -a reference_genome.fasta reads.fastq > all_reads_mapping.sam 
samtools fastq all_reads_mapping.sam > all_reads_mapped.fastq
```
5. For each contig, generate a fastq of mapped reads and read coverage details
``` 
samtools sort all_reads_mapping.sam -o all_reads_sorted_mapping.bam
samtools index all_reads_sorted_mapping.bam
samtools view -h all_reads_sorted_mapping.bam {contig_name} | samtools fastq - > {contig_name}_mapped.fastq
samtools view -h all_reads_sorted_mapping.bam {contig_name} | samtools depth - > {contig_name}_mapped_coverage.txt
```
Repeat for each contig in the genome

6. Calculate the 5th percentile coverage for each contig, I used R. 
```
###
In a R consol
###
samples <-c(1,2,3,4) #number of contigs
for (i in samples){
  print(i)
  cov<-read_delim(paste0("contig_",i,"_mapped_coverage.txt"),col_names=FALSE)
  print(quantile(cov$X3,probs=c(0.05,0.1)))
  }
```
Use this to decide on a coverage target for each sample. I used one coverage target for all samples within a strain.

7. Calculate the mean read lengths for each contigs reads 

```
seqkit stats {contig_name}_mapped.fastq 
```

8. Based on the read lengths and coverage target, estimate the number of reads to select per genomic window.  
Formula:  
```
(Window size/mean read length) * target coverage = estimated number of reads per window
```

9. Update *normalise_WGA_cov.Snakefile* with the target number of reads for each contig and run 
```
snakemake -s normalise_WGA_cov.Snakefile --cores 1 
```
** Only use 1 core ** despite each window being processed individually, the subsampling rule merges all fastq into a single output and if these are run in parallel the resulting fastq file will be corrupt.  

Output files will be a dir of marker files for each contig indicating each windows reads were subsampled, and a single fastq for each contig containing the subsampled reads. 

10. Confirm target coverage was achieved by mapping the subsampled reads onto each contig. 
```


