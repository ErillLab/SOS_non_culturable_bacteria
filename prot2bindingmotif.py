# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules

import os, subprocess, warnings
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Emboss.Applications import NeedleCommandline
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO, SeqIO, AlignIO, Entrez
from StringIO import StringIO


def tf2bindingmotif(email_address,database,input_fasta):
        
    'Function to get the putative binding motif(s) of a self-regulated Transcription Factor (TF)'
    
    # Tell to NCBI who I am

    Entrez.email = email_address
    
    
    ### Identify the putative TFs orthologs of a given database
    
    # Declare the thresholds
    
    COVERAGE = 75
    
    # Perform the BLAST search using reference TFs, an E-value of 1e-20, a coverage greater than 75 % and a maximum number of alignments of 1000
    
    blastp_cline = NcbiblastpCommandline(query=input_fasta, db=database, evalue=1e-20, outfmt='"6 std qcovs"', num_alignments = 1000)
    stdout, sterr = blastp_cline()
     
    # Parse the output to get the unique accession of each hit that has exceeded the threshold
    
    blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",fields=['std','qcovs'])
    validated_TF = []
    for blast_record in blast_records:
        for alignment in blast_record:
            if alignment.query_coverage >= COVERAGE: 
                acc = alignment.id_all[0]
                if not "_" in acc:                
                    validated_TF.append(acc)
                    validated_TF = list(set(validated_TF))
    
    
    ### Get the putative promoter (-250) of each identified TF

    # Search the db_source and the positions of each validated TF using protein information
    
    nucleotide, positions, strand = [], [], []    
    
    for TF in validated_TF:        
        
        if "WP_" in TF:
            try:
                records = Entrez.read(Entrez.efetch(db="protein", id=TF, rettype='ipg', retmode='xml'))
                for record in records["IPGReport"]["ProteinList"]:
                    r = record["CDSList"][0]
                    nucleotide.append(r.attributes["accver"])
                    if r.attributes["strand"] == "+":
                        positions.append(r.attributes["start"])
                        strand.append(1)
                    else:
                        positions.append(r.attributes["stop"])
                        strand.append(2)
                        
                    break
            except KeyError:
                continue
        
        else:
            handle = Entrez.efetch(db="protein", id=TF, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "gb")      
            
            for record in records:
                if "db_source" in record.annotations:
                    db_source = record.annotations["db_source"].split()
                for f in record.features:
                    if f.type == "CDS":
                        if "coded_by" in f.qualifiers:
                            position = f.qualifiers["coded_by"]
                            
                            # Annotate the positions and the strand of all validated TFs
                            
                            if "<" in position[0] or ">" in position[0]:
                                continue
                            elif "complement" in position[0]:
                                nucleotide.append(db_source[-1])
                                positions.append(position[0].split(")")[0].split("..")[1])
                                strand.append(2)
                            else:
                                nucleotide.append(db_source[-1])
                                positions.append(position[0].split(":")[1].split(".")[0])
                                strand.append(1)
                                
    # Search the putative promotor (-250) of each validated TF using db_source and positions supplied by the previous step
    
    out_handle = open("putativepromotor_validatedTF_%s.fasta"%database, "w")
    
    for i in range (0,len(nucleotide)):
        
        if strand[i] == 1:
            handle = Entrez.efetch(db="nucleotide", id = nucleotide[i], rettype="fasta", strand=strand[i], seq_start=int(positions[i])-250, seq_stop=int(positions[i]))
        
        else:
            handle = Entrez.efetch(db="nucleotide", id = nucleotide[i], rettype="fasta", strand=strand[i], seq_start=int(positions[i]), seq_stop=int(positions[i])+250)
        
        out_handle.write(handle.read())
        handle.close()
    out_handle.close()
    
    
    ### Remove > 90% similar promoters using the Needleman–Wunsch algorithm
    
    # Save the sequences and ids in two different lists
    
    validated_TF_sequence, validated_TF_id = [], []
    for record in SeqIO.parse("putativepromotor_validatedTF_%s.fasta"%database, "fasta"):
        validated_TF_sequence.append(record.seq)
        validated_TF_id.append(record.id.split(" ")[0])
    
    # Create an empty list to save the ids of similar sequences
    
    remove_id = []

    # Perform global alignments of each sequence against all other using needle

    for i in range(0,len(validated_TF_sequence)-1):
        
        if len(validated_TF_sequence[i]) > 10:
        
        # Ignore the IDs that are present in remove_id
    
            if validated_TF_id[i] not in remove_id:
                for j in range (i+1, len(validated_TF_sequence)):
                    if len(validated_TF_sequence[j]) > 10:                  
                        needle_cline = NeedleCommandline(asequence="asis:"+validated_TF_sequence[i], bsequence="asis:"+validated_TF_sequence[j], gapopen=10, gapextend=0.5, outfile = 'stdout')
                        stdout, stderr = needle_cline()
                        alignment = AlignIO.parse(StringIO(stdout), "emboss")
                        for needle_records in alignment:
                            query = list(needle_records[0].seq)
                            subject = list(needle_records[1].seq)
                            matches = [h for h, k in zip(query, subject) if h == k]
                            while '-' in matches: matches.remove('-')   
                            similarity = (float(len(matches))/len(query))*100
                            
                            # Save the ids of similar sequences (similarity >= 90%)                    
                            
                            if similarity >= 90:    
                                remove_id.append(validated_TF_id[j])
                                remove_id = list(set(remove_id))
        else:
            remove_id.append(validated_TF_id[i])
            remove_id = list(set(remove_id))
                
    # print in a FASTA file all sequences that are not present in remove_id
                        
    records = (r for r in SeqIO.parse("putativepromotor_validatedTF_%s.fasta"%database, "fasta") if r.id.split(" ")[0] not in remove_id)
    SeqIO.write(records, "putativeTF.fasta", "fasta")
    os.rename("putativeTF.fasta","putativepromotor_validatedTF_%s_clean.fasta"%database)
    os.remove("putativepromotor_validatedTF_%s.fasta"%database)


    ### Motif discovery using MEME
    
    # Run the MEME commands in bash
    
    bashCommand = "meme putativepromotor_validatedTF_%s_clean.fasta -dna -nostatus -time 18000 -mod anr -nmotifs 10 -minw 14 -maxw 20 -revcomp -maxsites 200 -pal -oc Motif_Discovery_Palindrome"%database
    subprocess.call(bashCommand.split())
    bashCommand = "meme putativepromotor_validatedTF_%s_clean.fasta -dna -nostatus -time 18000 -mod anr -nmotifs 10 -minw 14 -maxw 20 -revcomp -maxsites 200 -oc Motif_Discovery_Not_Palindrome"%database
    subprocess.call(bashCommand.split())
    

# Run the function

tf2bindingmotif("","","")
