# -*- coding: utf-8 -*-



# Load the necessary modules

from Bio import Entrez, SeqIO
import urllib2

from datetime import datetime
startTime=datetime.now()



# Tell to NCBI who I am

Entrez.email = "miquelsanchezosuna@gmail.com"


# FIMO output

ids = [ ]

strand = [ ]

positions = [ ]

motif = [ ]


def get_regulated_genes(species=None):
    
    for i in range(0,len(ids)):   
        
        while True:
            
            try:
                
                if strand[i] == "1":       
                    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", seq_start=int(positions[i])-66)
                elif strand[i] == "-1":
                    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", seq_start=0, seq_stop=int(positions[i])+50)
                break
        
            except urllib2.HTTPError:
                
                print "Oops! HTTP Error...Trying again"
                pass
            
        records = SeqIO.parse(handle, "gb")
        
        current_genes, current_position, current_id , intergenic_distance = [], [], [], []
        ATG_end, organism = None, None        
        
        for record in records:
            organism =  record.annotations['organism']
            if species == None or species == organism:
                for f in record.features:
                    if f.type == "CDS":
                        if str(f.strand) == strand[i]:
                            if f.strand == 1 and "<" in str(f.location) or f.strand == -1 and ">" in str(f.location):
                                continue
                            else:
                                try:                                
                                    protein, protein_id = f.qualifiers["product"], f.qualifiers["protein_id"]
                                    current_genes.append(protein[0])
                                    current_id.append(protein_id[0])
                                    if strand[i] == "1":
                                        current_position.append(int(str(f.location).split("[")[1].split(":")[0].strip("<").strip(">")))
                                    if strand[i] == "-1":
                                        current_position.append(int(str(f.location).split("[")[1].split(":")[1].split("]")[0].strip("<").strip(">")))
    
                                    if ATG_end != None:
                                        intergenic = int(str(f.location).split("[")[1].split(":")[0].strip("<").strip(">")) - ATG_end
                                        intergenic_distance.append(intergenic)
                                    ATG_end = int(str(f.location).split("[")[1].split(":")[1].split("]")[0].strip("<").strip(">"))
                                except KeyError:
                                    continue
        
        handle.close()

        if len(current_position) > 0:
            if strand[i] == "1" and current_position[0] <= 350:
                c = 0
                for j in intergenic_distance:
                    if len(intergenic_distance) == 0:                     
                        break
                    else:
                        if j > 150:                                 
                            break
                        else:
                            c+=1
                for k in range(0,c+1):
                    print organism+'\t'+ids[i]+'\t'+current_id[k]+'\t'+current_genes[k]
            elif strand[i] == "-1" and int(positions[i]) - current_position[::-1][0] <= 266: # 266 = 250 + 16 (motif size)
                current_genes, current_position, current_id , intergenic_distance = current_genes[::-1], current_position[::-1], current_id[::-1], intergenic_distance[::-1]         
                c = 0
                for j in intergenic_distance:
                    if len(intergenic_distance) == 0:                     
                        break
                    else:
                        if j > 150:                                 
                            break
                        else:
                            c+=1
                for k in range(0,c+1):
                    print organism+'\t'+ids[i]+'\t'+current_id[k]+'\t'+current_genes[k]
       


def get_reverse_regulated_genes(species=None):
            
    for i in range(0,len(ids)):   
        
        while True:
            
            try:
                
                if strand[i] == "-1":       
                    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", seq_start=int(positions[i])-66)
                elif strand[i] == "1":
                    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", seq_start=0, seq_stop=int(positions[i])+50)
                break
        
            except urllib2.HTTPError:
                
                print "Oops! HTTP Error...Trying again"
                pass
            
        records = SeqIO.parse(handle, "gb")
        
        current_genes, current_position, current_id , intergenic_distance = [], [], [], []
        ATG_end, organism = None, None        
        
        for record in records:
            organism =  record.annotations['organism']
            if species == None or species == organism:
                for f in record.features:
                    if f.type == "CDS":
                        if str(f.strand) != strand[i]:
                            if f.strand == 1 and "<" in str(f.location) or f.strand == -1 and ">" in str(f.location):
                                continue
                            else:
                                try:                                 
                                    protein, protein_id = f.qualifiers["product"], f.qualifiers["protein_id"]
                                    current_genes.append(protein[0])
                                    current_id.append(protein_id[0])
                                    if strand[i] == "-1":
                                        current_position.append(int(str(f.location).split("[")[1].split(":")[0].strip("<").strip(">")))
                                    if strand[i] == "1":
                                        current_position.append(int(str(f.location).split("[")[1].split(":")[1].split("]")[0].strip("<").strip(">")))
    
                                    if ATG_end != None:
                                        intergenic = int(str(f.location).split("[")[1].split(":")[0].strip("<").strip(">")) - ATG_end
                                        intergenic_distance.append(intergenic)
                                    ATG_end = int(str(f.location).split("[")[1].split(":")[1].split("]")[0].strip("<").strip(">"))
                                except KeyError:
                                    continue
        handle.close()

        if len(current_position) > 0:
            if strand[i] == "-1" and current_position[0] <= 350:
                c = 0
                for j in intergenic_distance:
                    if len(intergenic_distance) == 0:                     
                        break
                    else:
                        if j > 150:                                 
                            break
                        else:
                            c+=1
                for k in range(0,c+1):
                    print organism+'\t'+ids[i]+'\t'+current_id[k]+'\t'+current_genes[k]
            elif strand[i] == "1" and int(positions[i]) - current_position[::-1][0] <= 266: # 266 = 250 + 16 (motif size)
                current_genes, current_position, current_id , intergenic_distance = current_genes[::-1], current_position[::-1], current_id[::-1], intergenic_distance[::-1]         
                c = 0
                for j in intergenic_distance:
                    if len(intergenic_distance) == 0:                     
                        break
                    else:
                        if j > 150:                                 
                            break
                        else:
                            c+=1
                for k in range(0,c+1):
                    print organism+'\t'+ids[i]+'\t'+current_id[k]+'\t'+current_genes[k]
            
        
    
# Run the functions

get_regulated_genes()
get_reverse_regulated_genes()

# Print the datetime

print datetime.now() - startTime
