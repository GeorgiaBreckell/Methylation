#Parsing each Natural Isolate Genome to extract the location of target motifs. 

import re 
import sys

genome = sys.argv[1]
#Trim input file name to the strain name only
strain= genome.rstrip("_polished_genome.fasta")

#Set output file names
dam_output= "DAM_%s_positions.txt" % strain
dcm_output= "DCM_%s_positions.txt" % strain

#Motifs to look for: 
dam="GATC"
dcm="CC[AT]{1}GG"
#dcm regex will find instances of CCTGG or CCAGG. {} indicate to find only 1 instance of A or T. 

with open(genome) as g:
    sequence=g.read()
    dam_matches= re.finditer(dam,sequence) #returns an interatable object, used in the subsequent list comp
    dcm_matches= re.finditer(dcm,sequence)
    dam_positions=[match.start() for match in dam_matches]#List comphrehensions 
    dcm_positions=[match.start() for match in dcm_matches]
#Returns the index where each motif starts. 

with open(dam_output, "w") as txtfile: 
    txtfile.writelines("%s\n" % position for position in dam_positions)

with open(dcm_output, "w") as txtfile: 
    txtfile.writelines("%s\n" % position for position in dam_positions)


