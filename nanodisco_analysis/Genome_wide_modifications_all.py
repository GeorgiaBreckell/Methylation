#####
#Script which takes the filtered CSV files outputted from R and returns the number of methylated sites
#in a genome across a series of windows. 

#import re to search for phrases in the text. 
import csv
import sys
import re
import numpy as np
import pandas as pd 

#Input Files
genome = sys.argv[1] #  for snakemake
modified_sites = sys.argv[2] #

#Global variables
sample = re.sub("_[A-Z0-9]{3,4}_site_pvals.tab"," ", modified_sites)
sample = re.sub("results/nanodisco/motif_pvals_with_coverage/\w{2}/"," ",sample) #Trimmed sample name
sample = sample.strip()

print(sample)

MTase = re.sub("_site_pvals.tab"," ", modified_sites)
MTase = re.sub("results/nanodisco/motif_pvals_with_coverage/\w{2}/\w{2}_\w{1}_\w{2}_"," ",MTase)
MTase = MTase.strip()

print(MTase)

def contig_lengths(genome):
        with open(genome) as g:
                contigs=0
                contig_len_d={}
                for line in g:
                        if not line.startswith(">"):
                                contigs += 1
                                bp= line.strip("\n")
                                contig_size = len(bp)
                                contig_len_d[contigs]=contig_size
                return(contig_len_d)

#Creates a dictionary of contig lengths, keys are numbered contig counts, values are lengths in bp of each contig. 

def mod_site_detection(sites,modified_sites,window_count,contig):
    for line in range(len(sites)):
        #Check if the contig is the same as the contig being run, so it does this for EVERY line in the the file
        if (sites.at[line,"Site"] < win_end):
            if (sites.at[line,"Contig"].split("_")[1] == str(contig) and sites.at[line,"Site"] >= win_start and sites.at[line,"Site"] < win_end ):
                #print("contig is " + sites.at[line,"Contig"].split("_")[1]+ " position is " + str(sites.at[line,"Site"]))
                window_count += 1
                #print("Sites in window is " + str(window_count))
                    
                if (sites.at[line,"Site"] >= win_start and sites.at[line,"Site"] < win_end and sites.at[line,"Site_mod"]!= 0):
                    modified_sites += 1
                    #print("modified sites in window is " + str(modified_sites))
                    
    return(modified_sites,window_count)
    
                #print(window_count)
                #return(modified_sites)    


window=4000
jump=4000

contigs = contig_lengths(genome)
print("The contigs are " + str(contigs))
# Creates a dictionary variable to use later when looping through the genomes by contig to ensure the right number nad length of contigs are evaluated

#Take the modified sites CSV output from R and stores in a dictionary a list of all modified sites. 
#with open(modified_sites, newline='') as modified_csv:
    #csv_reader = csv.reader(modified_csv,delimiter="\t")
sites= pd.read_csv(modified_sites,sep="\t")
len_sites=len(sites)
print(sites.head(3))
line_count=0
modified_positions = []
#window_total_sites = {}
#window_modified_sites = {}


#I need to loop through each window? and run the checking of sites as a function that gets run for each window, 
#so when it exits it restarts from the first site in that window each time? 


for contig in contigs:
    window_total_sites = {}
    window_modified_sites = {}
    print("Running loop for contig " + str(contig))
    contig_size = contigs[contig] #the contig size is equal to the value in the contigs dictionary associated with that key.
    #print(contig_size) 
    win_start = 0  
    while win_start < contig_size:
        if win_start + window < contig_size:
            win_end = win_start + window
        else: 
            win_end = contig_size
        #print(win_end)
        modified_sites = 0 
        window_count = 0 
 
        mod_sites, total_sites =mod_site_detection(sites,modified_sites,window_count,contig)
        #print(mod_sites)
        #print(total_sites)       
                #if (sites.at[line,"Site"] >= win_start and sites.at[line,"Site"] < win_end):
                    #print("position is " + str(sites.at[line,"Site"]))
                    #window_count += 1
                    #print("Sites in window is " + str(window_count))
                
                #if (sites.at[line,"Site"] >= win_start and sites.at[line,"Site"] < win_end and sites.at[line,"Site_mod"]!= 0):
                    #print("1 means modified " + str(sites.at[line,"Site_mod"]))
                    #modified_sites += 1
                    #print("modified sites in window is " + str(modified_sites))
            
    
        window_total_sites[win_end]=total_sites
        window_modified_sites[win_end]=mod_sites
        win_start += jump 
        window_count = 0 
        modified_sites = 0 


    csv_headers= ['Window','Number.sites']
    modified_output_filename = "results/nanodisco/genome_wide_methylation/%s_window_%s_MS_%s_contig_%s.csv" % (window,MTase,sample,contig)
    print("writing output for contig " + str(contig))
    #print(window_modified_sites)
    with open(modified_output_filename, 'w') as csvfile:
        writer=csv.writer(csvfile)
        writer.writerows(window_modified_sites.items())

    total_output_filename = "results/nanodisco/genome_wide_methylation/%s_window_%s_TS_%s_contig_%s.csv" % (window,MTase,sample,contig)
    with open(total_output_filename, 'w') as csvfile:
        writer=csv.writer(csvfile)
        writer.writerows(window_total_sites.items())


#sites= pd.read_csv(modified_sites)










