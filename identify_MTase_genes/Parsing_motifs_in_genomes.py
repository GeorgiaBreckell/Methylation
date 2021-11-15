#Parsing each Natural Isolate Genome to extract the location of target motifs. 

import re 
import sys
import csv

genome = sys.argv[1]
#Trim input file name to the strain name only
strain= genome.rstrip("_polished_genome.fasta").lstrip("../nanodisco_genomes/")

#Set output file names
dam_output= "DAM_%s_positions.tab" % strain
dcm_output= "DCM_%s_positions.tab" % strain
KII_output="KII_%s_positions.tab" % strain
GIX_output="GIX_%s_positions.tab" % strain
GVI_output="GVI_%s_positions.tab" % strain
HK3I_output="HK3I_%s_positions.tab" % strain
Eco31I_output="Eco31I_%s_positions.tab" % strain

#Motifs to look for: 
dam="GATC"
dcm="CC[AT]{1}GG"
KII="ATGCAT"
GIX="[GC]{1}A[CT]{1}"
GVI="ATGCAT"
HK3I="[CT]{1}GGCC[AG]{1}"
Eco31I="GGTCTC"
#dcm regex will find instances of CCTGG or CCAGG. {} indicate to find only 1 instance of A or T. 

with open(genome) as g:
    seq=g.readlines()
    for l in range(len(seq)): 
        if (l %2 == 0 ):
            #the above two lines group the lines of the fasta by header and sequence
            contig=seq[l]
            contig=contig.lstrip(">")
            contig_list=contig.split(" ")
            contig=contig_list[0]
            contig=contig.strip("")
            #making the contig name tidy for the output (probably could be less lines)
            sequence=seq[l+1]

            dam_matches= re.finditer(dam,sequence) #returns an interitable object, used in the subsequent list comp
            dcm_matches= re.finditer(dcm,sequence)
            KII_matches=re.finditer(KII,sequence)
            GIX_matches=re.finditer(GIX,sequence)
            GVI_matches=re.finditer(GVI,sequence)
            HK3I_matches=re.finditer(HK3I,sequence)
            Eco31I_matches=re.finditer(Eco31I,sequence)
            
           

            dam_positions=[match.start() for match in dam_matches]#List comprehensions 
            dcm_positions=[match.start() for match in dcm_matches]#resulting variable is a list of start positions
            KII_positions=[match.start() for match in KII_matches]
            GIX_positions=[match.start() for match in GIX_matches]
            GVI_positions=[match.start() for match in GVI_matches]
            HK3I_positions=[match.start() for match in HK3I_matches]
            Eco31I_positions=[match.start() for match in Eco31I_matches]


            dam_list=[]
            dcm_list=[]
            KII_list=[]
            GIX_list=[]
            GVI_list=[]
            HK3I_list=[]
            Eco31I_list=[]

            for position in dam_positions:
                dam_list.append(contig+'\t'+str(position)+'\n')
                #making each line of the output including the new_line as an element in a list
            with open(dam_output, "a") as tabfile: 
                tabfile.writelines(dam_list)
                #writing the output from above list 
            
            for position in dcm_positions:
                dcm_list.append(contig+'\t'+str(position)+'\n')
            with open(dcm_output, "a") as tabfile: 
                tabfile.writelines(dcm_list)

            for position in KII_positions:
                KII_list.append(contig+'\t'+str(position)+'\n')
            with open(KII_output, "a") as tabfile: 
                tabfile.writelines(KII_list)

            for position in GIX_positions:
                GIX_list.append(contig+'\t'+str(position)+'\n')
            with open(GIX_output, "a") as tabfile: 
                tabfile.writelines(GIX_list)

            for position in GVI_positions:
                GVI_list.append(contig+'\t'+str(position)+'\n')
            with open(GVI_output, "a") as tabfile: 
                tabfile.writelines(GVI_list)

            for position in HK3I_positions:
                HK3I_list.append(contig+'\t'+str(position)+'\n')
            with open(HK3I_output, "a") as tabfile: 
                tabfile.writelines(HK3I_list)

            for position in Eco31I_positions:
                Eco31I_list.append(contig+'\t'+str(position)+'\n')
            with open(Eco31I_output, "a") as tabfile: 
                tabfile.writelines(Eco31I_list)

            #I think I should be using csv writer here, but I couldnt work out how to append the contig name at the 
            #start of each line and get the new line, if i had "\n" it would literally output " at the end of one line
            # and " at the start of the next, then insert the tab. 
            # contig_1_pilon    84"
            #"contig_1_pilon    103"
            #" <- like this 