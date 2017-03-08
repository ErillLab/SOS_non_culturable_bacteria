# -*- coding: utf-8 -*-



# Load the necessary modules

import os, warnings
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Emboss.Applications import NeedleCommandline
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO, SeqIO, AlignIO, Entrez
from StringIO import StringIO
from datetime import datetime
startTime=datetime.now()



# Tell to NCBI who I am

Entrez.email = "miquelsanchezosuna@gmail.com"


 
def first_blast(database):
        
    'Function to get the putative LexAs of a given database'
    
    # Declare the thresholds
    
    COVERAGE = 85
    
    # Perform the BLAST search using reference LexAs (sequence.fasta), an E-value of 1e-30, a coverage greater than 85 % and a maximum number of alignments of 1000
    
    blastp_cline = NcbiblastpCommandline(query="sequence.fasta", db=database, evalue=1e-30, outfmt='"6 std qcovs"', num_alignments = 1000)
    stdout, sterr = blastp_cline()
     
    # Parse the output to get the unique accession of each hit that has exceeded the thresholds
    
    blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",fields=['std','qcovs'])
    array = []
    for blast_record in blast_records:
        for alignment in blast_record:
            if alignment.query_coverage >= COVERAGE: 
                acc = str(alignment.id_all).split("|")[-2]
                array.append(acc)
                array = list(set(array))
    
    # Save the putative LexAs into a single FASTA file
    
    handle = Entrez.efetch(db="protein", id=array, rettype="fasta", retmode="text")
    out_handle = open("putativeLexA_%s.fasta"%database, "w")
    out_handle.write(handle.read())
    out_handle.close()
    handle.close()
    


def reciprocal_blast(database):
    
    'Function to validate the putative LexAs by reblasting the hits to the reference genomes'    
    
    # first_blast function is used to get the putative LexAs in a FASTA file
    
    first_blast(database)    
    
    # Declare the thresholds
    
    COVERAGE = 85
    
    # Number of putative LexAs
    
    putative_LexA, validated_LexA = [], []
    for record in SeqIO.parse("putativeLexA_%s.fasta"%database, "fasta"):
        putative_LexA.append(record.id)
        
    # Perform the local BLAST in order to obtain the best hit of each putative LexAs vs the reference genomes databases using thresholds of E-value < 1e-30 and coverage >= 85 %
    
    # While the putative LexAs are been validated they are deleted to avoid redundant validations
        
    db = ["Acidobacterium", "Agrobacterium", "Bsubtilis", "Bthuringiensis", "Bdellovibrio", "Burkholderia", "Caulobacter", "Clostridium", "Corynebacterium", "Dehalococcoides", "Deinococcus", "Ecoli", "Fibrobacter", "Geobacter", "Haemophilus", "Leptospira", "Listeria", "Magnetococcus", "Mleprae", "Msmegmatis", "Mtuberculosis", "Myxococcus", "Nostoc", "Paracoccus", "Pectobacterium", "Peptoclostridium", "Petrotoga", "Prochlorococcus", "Paeruginosa", "Pputida", "Rhizobium", "Rcapsulatus", "Rsphaeroides", "Salmonella", "Serratia", "Soneidensis", "Spiezotolerans", "Sinorhizobium", "Staphylococcus", "Streptomyces", "Synechococcus", "Synechocystis", "Thermotoga", "Vcholerae", "Vparahaemolyticus", "Xaxonopodis", "Xcampestris", "Xylella"]    
    reference = ["WP_041839776.1", "WP_010971587.1", "WP_003238209.1", "EAO54512.1", "WP_011165841.1", "WP_009891323.1", "WP_010919768.1", "WP_010965137.1", "WP_003857389.1", "WP_041223489.1", "WP_010889603.1", "WP_000646078.1", "WP_015732425.1", "WP_010942262.1", "WP_010869048.1", "WP_000654116.1", "WP_009930693.1", "WP_011713713.1", "WP_010908066.1", "WP_003894124.1", "WP_003899448.1", "WP_011554445.1", "WP_010999034.1", "WP_011748429.1", "WP_044203865.1", "", "CBE04725.1", "WP_012209053.1", "WP_041710637.1", "WP_003091196.1", "WP_010953131.1", "WP_012483680.1", "WP_013068070.1", "WP_002719150.1", "WP_000646079.1", "WP_025304453.1", "WP_011074200.1", "WP_020914776.1", "WP_003528303.1", "WP_001208760.1", "WP_003973219.1", "WP_011933578.1", "WP_010872400.1", "WP_015919952.1", "WP_000803693.1", "WP_005480871.1", "WP_005915033.1", "WP_011036298.1", "WP_023906043.1"]    
    c = 0    
    
    for strain in db:
        records = (r for r in SeqIO.parse("putativeLexA_%s.fasta"%database, "fasta") if r.id.split(" ")[0] not in validated_LexA)
        SeqIO.write(records, "putativeLexA.fasta", "fasta")
        os.rename("putativeLexA.fasta","putativeLexA_%s.fasta"%database)        
        blastp_cline = NcbiblastpCommandline(query="putativeLexA_%s.fasta"%database, db=strain, evalue=1e-30, outfmt='"6 std stitle qcovs"', num_alignments = 1)
        stdout, sterr = blastp_cline()
        
        blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",fields=['std', 'stitle', 'qcovs'])
        
        # The putative LexAs are validated when their best hit is a reference LexA
        
        for blast_record in blast_records:
                for alignment in blast_record:
                    if alignment.query_coverage >= COVERAGE and reference[c] == alignment.title.split("|")[3]: 
                        acc = alignment.query_id
                        validated_LexA.append(acc)
                        validated_LexA = list(set(validated_LexA))
        
        # Break the validation if all the putative LexAs have been validated
        c += 1        
        
        if len(validated_LexA) == len(putative_LexA):
            break

    
    # Return the validated LexAs
                  
    return validated_LexA



def get_upstream_region(database):
    
    'Function to get the putative promoter (-250,+50) of each validated LexA'
    
    # reciprocal_blast function is used to get the validated LexAs

    validated_LexA = reciprocal_blast(database)
    
    # Search the db_source and the positions of each validated LexA using protein information
    
    nucleotide, positions, strand = [], [], []    
    
    for LexA in validated_LexA:        
        
        if "WP_" in LexA:
            handle = Entrez.efetch(db="protein", id=LexA, rettype="ipg", retmode="xml")
            records = Entrez.parse(handle)
            
            for record in records:
                for identical_protein in range(0,len(record["RedundantGiList"])):
                    
                    try:                    
                        # Annotate the positions and the strand of all validated LexAs
						
                        if str(record["RedundantGiList"][identical_protein]["CDSList"][0].attributes["strand"]) == "+":
                            nucleotide.append(record["RedundantGiList"][identical_protein]["CDSList"][0].attributes["accver"])
                            positions.append(record["RedundantGiList"][identical_protein]["CDSList"][0].attributes["start"])
                            strand.append(1)
                        else:
                            nucleotide.append(record["RedundantGiList"][identical_protein]["CDSList"][0].attributes["accver"])
                            positions.append(record["RedundantGiList"][identical_protein]["CDSList"][0].attributes["stop"])
                            strand.append(2)
                    
                    # KeyError appears in Swiss-Prot records (no CDS region in Nucleotide)
                    
                    except KeyError:
                        continue
        
        else:
            handle = Entrez.efetch(db="protein", id=LexA, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "gb")      
            
            for record in records:
                if "db_source" in record.annotations:
                    db_source = record.annotations["db_source"].split()
                for f in record.features:
                    if f.type == "CDS":
                        if "coded_by" in f.qualifiers:
                            position = f.qualifiers["coded_by"]
                            
                            # Annotate the positions and the strand of all validated LexAs
                            
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
                                
    # Search the putative promoter (-250, +50) of each validated LexA using db_source and positions supplied by the previous step
    
    out_handle = open("putativepromotor_validatedLexA_%s.fasta"%database, "w")
    
    for i in range (0,len(nucleotide)):
        
        if strand[i] == 1:
            handle = Entrez.efetch(db="nucleotide", id = nucleotide[i], rettype="fasta", strand=strand[i], seq_start=int(positions[i])-250, seq_stop=int(positions[i])+50)
        
        else:
            handle = Entrez.efetch(db="nucleotide", id = nucleotide[i], rettype="fasta", strand=strand[i], seq_start=int(positions[i])-50, seq_stop=int(positions[i])+250)
        
        out_handle.write(handle.read())
        handle.close()
    out_handle.close()

#    # Return the correlation between protein and nucleotide ids
#
#    correlation = dict(zip(validated_LexA, zip(nucleotide, positions)))
#    return correlation



def remove_similars(database):  
    
    'Function to remove sequences with similarity >= 90 %'
    
    # get_upstream_region function is used to get the putative promotors in a FASTA file
    
    get_upstream_region(database)
    
    # Save the sequences and ids
    
    validated_LexA_sequence, validated_LexA_id = [], []
    for record in SeqIO.parse("putativepromotor_validatedLexA_%s.fasta"%database, "fasta"):
        validated_LexA_sequence.append(record.seq)
        validated_LexA_id.append(record.id.split(":")[0])
    
    # Create an empty list to save the ids of similar sequences
    
    remove_id = []

    # Perform global alignments of each sequence against all other using the Needleman\96Wunsch algorithm

    for i in range(0,len(validated_LexA_sequence)-1):
        
        if len(validated_LexA_sequence[i]) > 10:
        
        # Ignore the IDs that are present in remove_id
    
            if validated_LexA_id[i] not in remove_id:
                for j in range (i+1, len(validated_LexA_sequence)):
                    if len(validated_LexA_sequence[j]) > 10:                  
                        needle_cline = NeedleCommandline(asequence="asis:"+validated_LexA_sequence[i], bsequence="asis:"+validated_LexA_sequence[j], gapopen=10, gapextend=0.5, outfile = 'stdout')
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
                                remove_id.append(validated_LexA_id[j])
                                remove_id = list(set(remove_id))
        else:
            remove_id.append(validated_LexA_id[i])
            remove_id = list(set(remove_id))
                
    # Print in a FASTA file all sequences that are not present in remove_id
                        
    records = (r for r in SeqIO.parse("putativepromotor_validatedLexA_%s.fasta"%database, "fasta") if record.id.split(":")[0] not in remove_id)
    SeqIO.write(records, "putativeLexA.fasta", "fasta")
    os.rename("putativeLexA.fasta","putativepromotor_validatedLexA_%s_clean.fasta"%database)

    
    
# Run the function

remove_similars("")

# Print the datetime

print datetime.now() - startTime