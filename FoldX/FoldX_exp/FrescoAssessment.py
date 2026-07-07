from pathlib import Path
from Bio import Align
import pandas as pd
import numpy as np
import subprocess
import csv
import yasara

#input can be structure of heterodimer, but only one fasta sequence 

#general
work_dir = Path("/Users/yannikmeindl/FoldX/FoldX_exp")
F1 = "DGLSGAT"

#get_fasta_seq
fasta = work_dir / "rcsb_pdb_4XHF.fasta"

#get_pdb_seq
pdb = work_dir/"4XHF.pdb"

#exp_data
csv_file = work_dir/"experimental_data.csv"

#Results from FRESCO => list with all possible mutations and their respective ddG values 
MutEnergyList = pd.read_csv(
    "MutationEnergies_CompleteListcopy.tab", 
    sep = r"\s+",
    skiprows = 1, 
    header = None, 
    names = ["mutation", "ddG", "sd"]
    )

#--------------------Preperation of fasta and pdb sequences--------------------#
def get_fasta_seq(fasta):

    records = []
    current_header = None
    current_lines = []

    with open(fasta, "r") as f:
        for line in f:
            line = line.strip()

            #saves current header from last round and the current_lines that get joined in the process
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_lines)))

                #defines line that started with ">" as the new header and opens list for lines (saved as individual strings) 
                current_header = line
                current_lines = []

            #takes all lines that don't beginn with ">" and appends them to the current_lines list
            else:
                current_lines.append(line)

    #saves content of last round if there is no more line that starts with ">"
    if current_header is not None:
        records.append((current_header, "".join(current_lines)))

    return records


def get_pdb_seq(pdb):
    
    yasara.LoadPDB(pdb)
    yasara.OligomerizeObj("1", center="Yes", instance="No") #removes non biological multimers (from crystallisation)

    pdb_seq = work_dir / "pdb_sequences.txt" #variable for the file where yasara saves the sequences in the next line

    yasara.SaveSeqMol("All", pdb_seq, join="No") #saves sequences of each molecule as individual fasta in one .txt file 

    records = []
    current_header = None
    current_lines = []

    with open(pdb_seq, "r") as f:
        for line in f:
            line = line.strip()

            #saves current header from last round and the current_lines that get joined in the process
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_lines)))

                #defines line that started with ">" as the new header and opens list for lines (saved as individual strings) 
                current_header = line
                current_lines = []

            #takes all lines that don't beginn with ">" and appends them to the current_lines list
            else:
                current_lines.append(line)

    #saves content of last round if there is no more line that starts with ">"
    if current_header is not None:
        records.append((current_header, "".join(current_lines)))

    return records         

#--------------------Alignment of fasta and pdb sequences--------------------#

def get_alignment_info(fasta_seqs, pdb_seqs):

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    results = []
    best_per_fasta = {}

    #loop searches for best alignment of each sequence from the pdb
    #result: in case of 3 sequences from pdb, after the loop there will be 3 alignments
    for fasta_seq in fasta_seqs:

        for pdb_seq in pdb_seqs:

            #list of tupels with informations about the two aligned sequences that is automatically sorted by the scores
            #between to sequences there might be multiple ways to align them (usually there is one alignment that clearly has the highest score)
            alignments = aligner.align(fasta_seq[1], pdb_seq[1])

            best_align = alignments[0]

            #results is a list of dictionaries with all the relevant inforamtion about the alignment between a fasta_seq und the pdb_seqs
            results.append({
                "fasta_id": fasta_seq[0],
                "fasta_seq": fasta_seq[1],
                "pdb_chain": pdb_seq[0],
                "pdb_seq": pdb_seq[1],
                "score": best_align.score,
                "fasta_length": len(fasta_seq[1]),
                "pdb_length": len(pdb_seq[1]),
                "aligned_areas": best_align.aligned,
            })

    #finds the best matching sequence out of a pdb file per fasta
    for result in results:

        fasta_id = result["fasta_id"]

        if fasta_id not in best_per_fasta:
            best_per_fasta[fasta_id] = result 

        elif result["score"] > best_per_fasta[fasta_id]["score"]:
            best_per_fasta[fasta_id] = result 
    
    return(best_per_fasta)

#--------------------Sorts information from alignment in tables--------------------#

def get_alignment_table(fasta_seq, pdb_seq, alignment):
    # Grundtabelle: FASTA ist das Framework
    df = pd.DataFrame({
        "fasta_seq": list(fasta_seq),
        "pdb_seq": [""] * len(fasta_seq),
        "pdb_pos": pd.Series([pd.NA] * len(fasta_seq), dtype ="Int64"), 
    })

    df.index = range(1, len(fasta_seq) + 1)
    df.index.name = "fasta_pos"

    # *_blocks is a list of blocks with matching residues e.g fasta_blocks = [[5, 54], [78, 93]]
    fasta_blocks = alignment[0]
    pdb_blocks = alignment[1]

    # loop iterates through every matching block in the *_blocks lists
    for fasta_block, pdb_block in zip(fasta_blocks, pdb_blocks):
        fasta_start, fasta_end = fasta_block #defines the two numbers per block as the beginning and end of the matching sequence
        pdb_start, pdb_end = pdb_block

        #per every block, the numbers from fasta_start to fasta_end are matched with the numbers from pdb_start to pdb_end
        #the lineup is then used to correctly place the pdb_seq values and indices correctly on the fasta_index in the dataframe
        for fasta_index, pdb_index in zip(
            range(fasta_start, fasta_end),
            range(pdb_start, pdb_end)
        ):
            df.loc[fasta_index +1, "pdb_seq"] = pdb_seq[pdb_index]
            df.loc[fasta_index +1, "pdb_pos"] = pdb_index + 1

    return df

#goes through the best_per_fasta dictioniary and aligns the sequences of the best pdb sequence with each fasta sequence
#creates dictionary with tables of the aligned sequnences (only one table if protein is not a heterodimer)
def concat_alignment_tables(best_per_fasta):

    alignment_tables = {}

    #goes through entries in the best_per_fasta dictionary
    for fasta_id in best_per_fasta.keys():

        #takes necessary input for the "get_alignment_table" function out of the dictionary entry
        fasta_s = best_per_fasta[fasta_id]["fasta_seq"]
        pdb_s = best_per_fasta[fasta_id]["pdb_seq"]
        alignment = best_per_fasta[fasta_id]["aligned_areas"]

        #executes the get_alignment_table function 
        alignment_table = get_alignment_table(fasta_s, pdb_s, alignment)

        #places the generated table in a dictionary 
        alignment_tables[fasta_id] = alignment_table

    return alignment_tables

def save_alignment_tables(alignment_tables):

    for fasta_id, alignment_table in alignment_tables.items():

        seq_name = fasta_id.split("|")[0]
        alignment_table.to_csv(f"{seq_name}_align_table.tab", sep = "\t")


#--------------------Prepares dictionary containing experimentally tested motifs and respective efficiencies--------------------# 
  
def open_exp_data(csv_file): 
    
    exp_data = pd.read_csv(csv_file, sep = ";", decimal = ",")
    
    return exp_data


def get_motifs(exp_data):

    all_motifs = {
        "DGLSGAT": 1
        }

    #takes information about mutations within F1 from exp_data list, alters F1 string and appends it to motif list together with exp. tested loading efficiencies
    for index, row in exp_data.iterrows():

        mut_pos = row["Position"]
        motif = f"{F1[:mut_pos]}{row["Mutation"]}{F1[(mut_pos+1):]}"

        all_motifs[motif] = row["Efficiency"]

    return all_motifs

#--------------------Compares motifs to sequence and calculates ddG of theoretical motif implementation at each position--------------------#

def get_ddgs(alignment_tables, all_motifs, MutEnergyList):
#output structure: 

    #dictionary with list for every 7 mer, that contains information about the mutation, ddg and sd of every individual position 
    ddg_per_mutation = {}                   
    #dictionary with lists for every pos containing total ddg and sd for every 7mer 
    ddg_per_position = {}
    #dictionary with dictionaries for every motif 
    ddg_per_motif = {}
    
    #iterates through entries of "alignment_tables" dictionary, sorted by fastas
    for fasta_id, alignment_table in alignment_tables.items():

        #makes a string out of the "pdb_seq" entries in the alignment table 
        pdb_seq = "".join(alignment_table["pdb_seq"].replace("", " ").fillna(" ").astype(str))

        #iterates through list of experimetnally testes motifs
        for motif, efficiency in all_motifs.items():

            #iterates through positions of all possible 7mers
            for pos in range(len(pdb_seq) - len(F1) + 1):

                #defines the window for the current position 
                window = pdb_seq[pos:pos+len(F1)]

                #iterates through the positions of the current window
                #counts and pairs the positions 
                for i, (pbd_aa, motif_aa) in enumerate(zip(window, F1)):

                    #brings information of the mismatch into the format of the mutation that would need to happen to implement the motif 
                    mut = f"{pbd_aa}{pos+i}{motif_aa}"

                    #tests if all the current position of the window actually contains an aa
                    if pbd_aa == " ":
                        ddg_per_mutation.pop(pos, None)
                        break
                    
                    #tests if theres a mismatch 
                    #if so it takes the respective ddg and sd values for that mismatch/mutation and appends them to a list 
                    elif pbd_aa != motif_aa:

                        #list of lines in the MutationEnergiesList whose mutation equals "mut" (=> supposed to be one)
                        line_of_mut = MutEnergyList.loc[MutEnergyList["mutation"] == mut]
                        
                        if line_of_mut.empty:
                            continue#raise value error 

                        else:

                            #defines the current row of the first mutation match tha has been found in the list (it is only one match anyways)
                            row = line_of_mut.iloc[0]

                            #creates a new list everytime it runs through a new position iteration
                            #then appends values from the MutationEnergiesList for every lower level iteration
                            ddg_per_mutation.setdefault(pos, []).append({
                                "mutation": row["mutation"], 
                                "ddG": float(row["ddG"]),
                                "sd": float(row["sd"]),
                            })

                    else:
                        continue
            
            #sums up the ddg + sd values of every 7mer and puts them into a dictionary
            #=> one dictionary for every motif
            for pos, mutations in ddg_per_mutation.items():

                ddg_per_position.setdefault(pos, []).append({
                    "ddG": round(sum(mut["ddG"] for mut in mutations), 4),
                    "sd": round(sum(mut["sd"] for mut in mutations), 4),
                }) 

            #puts the dictionaries of all the motifs into one dictionary 
            ddg_per_motif[motif] = ddg_per_position

    return(ddg_per_motif)



fasta_seqs = get_fasta_seq(fasta)
pdb_seqs = get_pdb_seq(pdb)
best_per_fasta = get_alignment_info(fasta_seqs, pdb_seqs)
alignment_tables = concat_alignment_tables(best_per_fasta)
save_alignment_tables(alignment_tables)

exp_data = open_exp_data(csv_file)
all_motifs = get_motifs(exp_data)
#print(all_motifs)

ddg_per_motif = get_ddgs(alignment_tables, all_motifs, MutEnergyList)
print(ddg_per_motif)
















