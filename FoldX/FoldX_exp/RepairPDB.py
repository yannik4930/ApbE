from pathlib import Path
from Bio import Align
import pandas as pd
import numpy as np
import subprocess
import csv
import yasara
import re

#general
work_dir = Path("/Users/yannikmeindl/FoldX/FoldX_exp")
F1 = "DGLSGAT"

#get_fasta_seq
fasta = work_dir / "rcsb_pdb_4XHF.fasta"

#get_pdb_seq
pdb = work_dir/"frescoTest"/"4XHF_cleaned.pdb"

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
    #yasara.OligomerizeObj("1", center="Yes", instance="No") #removes non biological multimers (from crystallisation)

    resultlist = yasara.ListMol("Obj 1") #creates list with all molecules of obj 1 
    for i, a in enumerate(resultlist, start = 1): # assigns index to all molecules 
        yasara.NameMol(f"Obj 1 Mol {a}, {i}") #replaces old_name with index

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
    #only one entry
    for result in results:

        fasta_id = result["fasta_id"]

        if fasta_id not in best_per_fasta:
            best_per_fasta[fasta_id] = result 

        elif result["score"] > best_per_fasta[fasta_id]["score"]:
            best_per_fasta[fasta_id] = result 
    
    return(best_per_fasta)

#--------------------Writes string for yasara object selection from fasta header of aligned best aligned sequence--------------------#

def yasara_selector_from_fasta_header(best_per_fasta):

    for fasta_id, result in best_per_fasta.items():
        header = result["pdb_chain"].strip() #gets rid of spaces etc. 

        if header.startswith(">"):
            header = header[1:].strip() #strips > and potential spaces after it

        obj_match = re.search(r"\bObject\s+(\d+)\b", header, re.IGNORECASE) #searches for "object x" in string
        if obj_match is None:
            raise ValueError(f"No object number found in header: {header}")
        mol_match = re.search(r"\bMolecule\s+(\d+)\b", header, re.IGNORECASE) #searches for "molecule x" in string
        if mol_match is None:
            raise ValueError(f"No molecule number found in header: {header}")
        
        obj_num = int(obj_match.group(1)) #takes second part of obj_num string (ususally 1 => Obj number) 
        mol_num = (mol_match.group(1)) #same for molecules

        return(f"Obj {obj_num} Mol {mol_num}")

#--------------------Reassigns residue numbers of pdb based on numbers from alignment--------------------#

def reassign_pdb_residues(pdb, best_per_fasta): #output in form of prepared .pdb file in work_dir

    yasara.Clear()
    yasara.LoadPDB(pdb)
    #yasara.OligomerizeObj("1", center="Yes", instance="No") #removes non biological multimers (from crystallisation)

    resultlist = yasara.ListMol("Obj 1") #creates list with all molecules of obj 1 
    for i, a in enumerate(resultlist, start = 1): # assigns index to list of molecules 
        yasara.NameMol(f"Obj 1 Mol {a}, {i}") #replaces old mol name with index

    selector = yasara_selector_from_fasta_header(best_per_fasta) 
    first_res_of_pdb_seq = list(best_per_fasta.values())[0]["aligned_areas"][0][0][0] + 1 #takes first residue of pdb compared to 

    yasara.NumberRes(selector, first = first_res_of_pdb_seq) #Reassigns residue numbers of pdb based on numbers from alignment
    prepared_pdb = work_dir / f"{pdb.stem}_prepared.pdb" #variable for place to save the prepared pdb
    yasara.SavePDB("Obj 1", prepared_pdb)


fasta_seqs = get_fasta_seq(fasta)

pdb_seqs = get_pdb_seq(pdb)

renamed_mol_pdb = work_dir / f"{pdb.stem}_renamed.pdb" 

best_per_fasta = get_alignment_info(fasta_seqs, pdb_seqs)

reassign_pdb_residues(pdb, best_per_fasta)



