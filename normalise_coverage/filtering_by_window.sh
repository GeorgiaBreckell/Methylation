#!/bin/bash


# $1 = all reads mapped bam file
# $2 = window
# $3 = output file


#Manipulating the  window variable to allow use in the awk command

IFS='-'
read -a splitwin <<< $2

IFS=':'
read -a splitcont <<< ${splitwin[0]}

echo "${splitcont[1]}"
echo "${splitwin[0]}"
echo "${splitwin[1]}"

#Filteriing reads 
echo "samtools view -h $1 $2 | awk -v  FS='\t' '$1 ~"^@" || ($4 >= '"${splitcont[1]}"' && $4 <= '"${splitwin[1]}"' )' | samtools fastq  - >  $3"

samtools view -h $1  $2 | awk -v  FS='\t' '$1 ~"^@" || ($4 >= '"${splitcont[1]}"' && $4 <= '"${splitwin[1]}"' )' | samtools fastq  - >  $3
