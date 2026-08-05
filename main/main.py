from pathlib import Path
from Bio import Align
from Bio.SeqUtils import seq1
import pandas as pd
import numpy as np
import subprocess
import csv
import yasara
import re
import os
import time 
import sys
import shutil
import json 
import math
import argparse

#general
BASE_DIR = Path(__file__).resolve().parent
F1 = "DGLSGAT"
window_size = len(F1)
foldx_path = "/scratch/s6765211/FoldX/foldx5_Linux/foldx_20270131"
position = 2 #input option for output mode => gives best scores/ddg values at input position"

#get_fasta_seq
fasta = BASE_DIR / "rcsb_pdb_1G3F.fasta"

#get_pdb_seq
pdb = BASE_DIR / "1G3F.pdb"

#exp_data
csv_file = BASE_DIR/"experimental_data.csv"

###pyFRESCO
pyfresco_dir = BASE_DIR / "pyFRESCO"
pyfresco_dir.mkdir(exist_ok = True)

#non-mutable region
far_enough_res = None
far_enough_zone = None

#switch for repair function 
repair = True

#open yasara in text mode (won't work on cluster otherwise)
yasara.info.mode = "txt"


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

    if len(records) > 1:
        raise ValueError(
            f"The fasta file contains {len(records)} sequences."
            "Expected is one sequence."
        )
    
    if len(records) == 0:
        raise ValueError(
            "The fasta file does not contain a sequence."
        )

    return records


def load_substitution_matrix(matrix_file):
    matrix = {}

    with open(matrix_file, "r", encoding="utf-8") as file:
        lines = file.readlines()

    header = lines[0].split()

    for line in lines[1:]:
        parts = line.split()

        row_aa = parts[0]
        scores = parts[1:]

        for col_aa, score in zip(header, scores):
            matrix[(row_aa, col_aa)] = int(score)

    return matrix


position_scoring = {
    1: {
        "matrix": load_substitution_matrix("matrices/BLOSUM62.txt"), 
        "a": 10.505,
        "b": -7.8119,
        "R2": 0.5876
    }, 
    2: {
        "matrix": load_substitution_matrix("matrices/GRANTHAM.txt"), 
        "a": -362.82,
        "b": 384.49,
        "R2": 0.9424
    },
    3: {
        "matrix": load_substitution_matrix("matrices/BLOSUM62.txt"), 
        "a": 3.6534,
        "b": -3.2351,
        "R2": 0.1255
    },
    4: {
        "matrix": load_substitution_matrix("matrices/BLOSUM45.txt"), 
        "a": 4.1893,
        "b": -2.5435,
        "R2": 0.5426
    },
    5: {
        "matrix": load_substitution_matrix("matrices/PAM160.txt"), 
        "a": 15.515,
        "b": -14.111,
        "R2": 0.5691
    },
    6: {
        "matrix": load_substitution_matrix("matrices/BLOSUM62.txt"), 
        "a": 5.603,
        "b": -3.1112,
        "R2": 0.8162
    },
}      

#--------------------Prepares dictionary containing experimentally tested motifs and respective efficiencies--------------------# 

def get_motifs(csv_file):

    exp_data = pd.read_csv(csv_file, sep = ";", decimal = ",")

    all_motifs = {
        "DGLSGAT": 1
        }

    #takes information about mutations within F1 from exp_data list, alters F1 string and appends it to motif list together with exp. tested loading efficiencies
    for index, row in exp_data.iterrows():

        mut_pos = row["Position"]
        motif = f"{F1[:mut_pos]}{row['Mutation']}{F1[(mut_pos+1):]}"

        all_motifs[motif] = row["Efficiency"]

    return all_motifs

# sliding window mechanism splits sequence into 7mers (windows)
# creates dictionary with a list of windows for each number of mismatches
def count_mismatches(sequence, all_motifs, window_size, max_mismatches: str):

    windows = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        windows.append(window) 

    mismatch_dict = {}

    for motif in all_motifs:

        mismatches_per_motif = {
                0: [],
                1: [], 
                2: [], 
                3: [], 
                4: [],
                5: [],
                6: [],
                7: [],
            }
    
        for index, window in enumerate(windows):
            mismatches = 0
            mismatch_info = []

            for pos, (a, b) in enumerate(zip(window, motif)):
                if a != b:
                    mismatches += 1
                    mismatch_info.append(((pos + 1), a, b))
        
            if mismatches <= max_mismatches:
                mismatches_per_motif[mismatches].append(((index + 1), window, mismatch_info))

        mismatch_dict[motif] = mismatches_per_motif

    return mismatch_dict


def load_mutation_data(csv_file):
    mutation_lookup = {}

    with open(csv_file, newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file, delimiter=";")

        for row in reader:
            position = int(row["Position"])
            original = row["Original"]
            mutation = row["Mutation"]
            efficiency = float(row["Efficiency"].replace(",", "."))

            mutation_lookup[(position, original, mutation)] = efficiency

    return mutation_lookup

#calculcates efficiency of mismatch w/o exp_data and then weighs it depending on the R2 
def get_efficiency(matrix_value, a, b, R2):
    x = ((matrix_value - b)/a)*R2 
    return x

#walks through mismatch_dict and assigns score to each window, depending on number of mismatches, existance of exp_data and substition matrices 
def get_assessment_score(mismatch_dict, exp_data, position_scoring):

    a_scores_per_motif = {}

    for motif, mismatches_per_motif in mismatch_dict.items():

        scored_motifs = [] 
    
        for mismatch_count, motif_hits in mismatches_per_motif.items():
        
            for pos, window_sequence, mismatch_details in motif_hits:
                assessment_score = 0 
                assessment_score += window_size - mismatch_count
            
                for mm_pos, window_aa, original_aa in mismatch_details:
                    mismatch = (mm_pos, window_aa, original_aa)

                    if mismatch in exp_data:
                        efficiency = exp_data[mismatch]
                        assessment_score += efficiency

                    else: 
                        if mm_pos == 7:
                            assessment_score += 0
                        
                        else:
                            best_matrix = position_scoring[mm_pos]["matrix"]
                            matrix_val = best_matrix[(window_aa, original_aa)]
                            matrix_based_efficiency = get_efficiency(matrix_val, position_scoring[mm_pos]["a"], position_scoring[mm_pos]["b"], position_scoring[mm_pos]["R2"])
                            capped_efficiency = min(matrix_based_efficiency, 1) 
                            assessment_score += capped_efficiency

                normalized_assessment_score = assessment_score / len(F1)
                rounded_assessment_score = round(normalized_assessment_score, 5)
                        
                scored_motifs.append((pos, window_sequence, rounded_assessment_score, mismatch_details))

        a_scores_per_motif[motif] = scored_motifs

    return(a_scores_per_motif)


def get_confidence_score(a_scores_per_motif, exp_data, position_scoring):

    ac_scores_per_motif = {}

    for motif, scored_motifs in a_scores_per_motif.items():

        scored_motifs_2 = {} 
        
        for pos, window_sequence, assessment_score, mismatch_details in scored_motifs:
            confidence_score = 0 
            confidence_score += window_size - len(mismatch_details)
            
            for mm_pos, original_aa, window_aa in mismatch_details:
                mismatch = (mm_pos, original_aa, window_aa)

                if mismatch in exp_data:
                    confidence_score += 1

                else:
                    if mm_pos == 7:
                        confidence_score += 1

                    else:
                        R2 = position_scoring[mm_pos]["R2"]
                        confidence_score += R2

            normalized_confidence_score = confidence_score / len(F1)
            rounded_confidence_score = round(normalized_confidence_score, 5)

            scored_motifs_2[pos] = {
                "Sequence Score": assessment_score, 
                "Confidence": rounded_confidence_score, 
                "MM": len(mismatch_details),
                "Window Sequence": window_sequence,
                "Details": mismatch_details
            }

        ac_scores_per_motif[motif] = scored_motifs_2
            
    return(ac_scores_per_motif)


def make_SS_output_tables_per_motif(ac_scores_per_motif, output_mode):

    SS_output_dict = {}

    for motif, scores_per_motif in ac_scores_per_motif.items():

        motif_scores = (
            pd.DataFrame.from_dict(
                scores_per_motif,
                orient = "index"
            )
            .rename_axis("Pos")
        )

        motif_scores.index = range(1, len(motif_scores) + 1)
        motif_scores.index.name = "Pos"

        if output_mode == "verbose":

            output_df = motif_scores

        elif output_mode == "min":

            output_df = motif_scores[
                ["Sequence Score", "MM"]
            ]
            
        else:
            
            output_df = motif_scores[
                ["Sequence Score", "Confidence", "MM", "Window Sequence"]
            ]

        SS_output_dict[motif] = output_df

    return SS_output_dict


#________________________________ENERGY SCORER__________________________________#

#--------------------Preperation of fasta and pdb sequences--------------------#

def get_pdb_seq(pdb):

    yasara.LoadPDB(pdb)
    yasara.OligomerizeObj("1", center="Yes", instance="No") #removes non biological multimers (from crystallisation)

    resultlist = yasara.ListMol("Obj 1") #creates list with all molecules of obj 1 
    for i, a in enumerate(resultlist, start = 1): # assigns index to all molecules 
        yasara.NameMol(f"Obj 1 Mol {a}, {i}") #replaces old_name with index

    pdb_seq = BASE_DIR / "pdb_sequences.txt" #variable for the file where yasara saves the sequences in the next line

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

#--------------------Sorts information from alignment in tables--------------------#

def get_alignment_table(fasta_seq, pdb_seq, alignment):
    # Grundtabelle: FASTA ist das Framework
    alignment_table = pd.DataFrame({
        "fasta_seq": list(fasta_seq),
        "pdb_seq": [""] * len(fasta_seq),
        "pdb_pos": pd.Series([pd.NA] * len(fasta_seq), dtype ="Int64"), 
    })

    alignment_table.index = range(1, len(fasta_seq) + 1)
    alignment_table.index.name = "fasta_pos"

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
            alignment_table.loc[fasta_index +1, "pdb_seq"] = pdb_seq[pdb_index]
            alignment_table.loc[fasta_index +1, "pdb_pos"] = pdb_index + 1

    return alignment_table

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

#--------------------Selects best aligning object/molecule to fasta sequence in yasara--------------------#

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

#--------------------Reassigns residue numbers of pdb based on alignment--------------------#

def reassign_pdb_residues(pdb, best_per_fasta): #output in form of prepared .pdb file in work_dir

    pdb = Path(pdb)

    yasara.Clear()
    yasara.LoadPDB(pdb)
    yasara.OligomerizeObj("1", center="Yes", instance="No") #removes non biological multimers (from crystallisation)

    resultlist = yasara.ListMol("Obj 1") #creates list with all molecules of obj 1 
    for i, a in enumerate(resultlist, start = 1): # assigns index to list of molecules 
        yasara.NameMol(f"Obj 1 Mol {a}, {i}") #replaces old mol name with index

    selector = yasara_selector_from_fasta_header(best_per_fasta) 
    first_res_of_pdb_seq = list(best_per_fasta.values())[0]["aligned_areas"][0][0][0] + 1 #takes first residue of pdb compared to fasta

    yasara.NumberRes(selector, first = first_res_of_pdb_seq) #Reassigns residue numbers of pdb based on numbers from alignment
    prepared_pdb = BASE_DIR / f"{pdb.stem}_prepared.pdb" #variable for place to save the prepared pdb
    yasara.SavePDB("Obj 1", prepared_pdb)
    
    return(prepared_pdb)

#--------------------Fix protonation states--------------------#

def repair_pdb(prepared_pdb, pdb):

    pdb = Path(pdb)

    yasara.Clear()
    yasara.LoadPDB(prepared_pdb)
    yasara.CleanAll()
    yasara.OptHydAll(method = "Yasara")
    prepaired_pdb = BASE_DIR/f"{pdb.stem}_prepaired.pdb"
    yasara.SavePDB("Obj 1", prepaired_pdb)
    os.remove(prepared_pdb)

#--------------------Execute pyFRESCO to get ddG of all possible mutations and prepares resulting rawMutEnergyList--------------------#

def make_rawMutEnergyList(work_dir, pyfresco_dir, prepaired_pdb, foldx_path, far_enough_res, far_enough_zone):

    subfolder_pdb = pyfresco_dir / prepaired_pdb.name
    selection_tab = pyfresco_dir / f"{prepaired_pdb.stem}.tab"

    shutil.copy(prepaired_pdb, subfolder_pdb)
    shutil.copy(work_dir / "DistributeFoldx.py", pyfresco_dir / "DistributeFoldx.py")
    shutil.copy(work_dir / "submit_fresco.sh", pyfresco_dir / "submit_fresco.sh")
    shutil.copy(work_dir / "fresco_job.sh", pyfresco_dir / "fresco_job.sh")

    far_enough_command = [
        sys.executable, 
        "FarEnoughZone.py", 
        f"--pdb={prepaired_pdb}",
        f"--output={selection_tab}", 
    ]

    if far_enough_res is not None:
        far_enough_command.append(f"--residue={far_enough_res}")

    if far_enough_zone is not None:
        far_enough_command.append(f"--distance={far_enough_zone}")

    subprocess.run(
        far_enough_command,
        cwd = work_dir, 
        check=True,
    )

    subprocess.run([
        sys.executable,
        "DistributeFoldx.py", 
        "Phase1", 
        prepaired_pdb.name,
        f"{prepaired_pdb.stem}.tab",
        "500",
        foldx_path
    ], cwd = pyfresco_dir, check = True)

    subprocess.run(
        ["bash", "submit_fresco.sh"],
        cwd = pyfresco_dir, 
        check = True
    )

    subprocess.run([
        sys.executable, 
        "DistributeFoldx.py", 
        "Phase2",
        prepaired_pdb.name,
        "-5"
    ], cwd = pyfresco_dir, check = True)

    return pyfresco_dir / "MutationEnergies_CompleteList.tab"


def get_MutEnergyList(rawMutEnergyList):

    #makes rawMutEnergyList a pandas object
    MutEnergyList = pd.read_csv(
        rawMutEnergyList, 
        sep = r"\s+",
        skiprows = 1, 
        header = None, 
        names = ["mutation", "ddG", "sd"]
        )
    return MutEnergyList

#--------------------Compares motifs to sequence and calculates ddG of theoretical motif implementation at each position--------------------#

def get_ddgs(alignment_tables, all_motifs, MutEnergyList):

    #dictionary with dictionaries for every motif 
    ddg_per_motif = {}
    
    #iterates through entries of "alignment_tables" dictionary, sorted by fastas
    for fasta_id, alignment_table in alignment_tables.items():

        #makes a string out of the "pdb_seq" entries in the alignment table 
        pdb_seq = "".join(alignment_table["pdb_seq"].replace("", " ").fillna(" ").astype(str))

        #iterates through list of experimentally testes motifs
        for motif, efficiency in all_motifs.items():

            #dictionary with list for every 7 mer, that contains information about the mutation, ddg and sd of every individual position 
            ddg_per_mutation = {}                   
            #dictionary with lists for every pos containing total ddg and sd for every 7mer 
            ddg_per_position = {}

            #iterates through positions of all possible 7mers
            for pos in range(len(pdb_seq) - len(motif) + 1):

                #defines the window for the current position 
                window = pdb_seq[pos:pos+len(motif)]

                #iterates through the positions of the current window
                #counts and pairs the positions 
                for i, (pbd_aa, motif_aa) in enumerate(zip(window, motif)):

                    #brings information of the mismatch into the format of the mutation in the FRESCO output
                    mut = f"{pbd_aa}{pos+i+1}{motif_aa}"

                    #tests if all the current position of the window actually contains an aa
                    if pbd_aa == " ":
                        ddg_per_mutation.pop(pos, None)
                        break
                    
                    #tests if theres a mismatch 
                    #if so it takes the respective ddg and sd values for that mismatch/mutation and appends them to a list 
                    elif pbd_aa != motif_aa:

                        #list of lines in the MutationEnergiesList whose mutation equals "mut" (=> supposed to be one)
                        line_of_mut = MutEnergyList.loc[MutEnergyList["mutation"] == mut]
                        
                        #checks if current mutation is listed in mutation energies list
                        #breaks window loop if not 
                        if line_of_mut.empty:
                            ddg_per_mutation.pop(pos, None)
                            break
                            
                        else:

                            #defines the current row of the first mutation match tha has been found in the list (it is only one match anyways)
                            row = line_of_mut.iloc[0]

                            #creates a new list everytime it runs through a new position iteration
                            #then appends values from the MutationEnergiesList for the current mutation within the 7mer
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
                
                #pos+1 so that output shows biological residue number instead of position in python string
                ddg_per_position[pos+1] = {
                    "ddG": round(sum(mut["ddG"] for mut in mutations), 4),
                    "sd": round(sum(mut["sd"] for mut in mutations), 4),
                }

            #puts the dictionaries of all the motifs into one dictionary 
            ddg_per_motif[motif] = ddg_per_position

    with open("ddg_per_motif.json", "w") as file:
            json.dump(ddg_per_motif, file, indent = 4)

    return(ddg_per_motif)

#--------------------Calculates scores from ddg values--------------------# 

def get_energy_scores(ddg_per_motif):

    motif_scores = {}

    for motif, positions in ddg_per_motif.items():

        residue_scores = {}

        for residue, values in positions.items():

            score = 0

            ddg_value = values["ddG"]
            ddg_sd = values["sd"]

            #scoring of ddg_values
            if 0 < ddg_value < 5:
                score = 1 - 0.04 * ddg_value

            elif ddg_value > 5:
                score = 0.5 + 0.3 * math.exp(-0.5 * (ddg_value - 5))

            elif 0 > ddg_value > -5: 
                score = 1 + 0.04 * ddg_value

            elif ddg_value < -5: 
                score = 0.5 + 0.3 * math.exp(0.5 * (ddg_value + 5))

            score = round(score, 5)

            residue_scores[residue] = {
                "Energy Score": score, 
                "ddG": ddg_value,
                "sd": ddg_sd 
            }

        motif_scores[motif] = residue_scores    
    
    return motif_scores

#--------------------Makes ES output tables dependent on output mode--------------------# 

def make_ES_output_tables_per_motif(energy_scores_per_motif, output_mode):

    ES_output_dict = {}

    for motif, scores_per_motif in energy_scores_per_motif.items():

        motif_scores = (
            pd.DataFrame.from_dict(
                scores_per_motif,
                orient = "index"
            )
            .rename_axis("Pos")
        )
        
        motif_scores.index.name = "Pos"

        if output_mode == "verbose":

            output_df = motif_scores

        elif output_mode == "min":

            output_df = motif_scores[
                ["Energy Score"]
            ]
            
        else:
            
            output_df = motif_scores

        ES_output_dict[motif] = output_df

    return ES_output_dict


#________________________________SECONDARY STRUCTURE SCORER__________________________________#

#--------------------Extracts information about secondary structure of each position from pdb--------------------#
def getSecStr(prepaired_pdb):

    yasara.Clear()
    yasara.LoadPDB(prepaired_pdb)

    ca_atoms = yasara.ListAtom("CA Protein Obj 1")

    secstr_dict = {}

    for atom in ca_atoms:
        selection = f"Atom {atom}"

        chain = yasara.NameMol(selection)

        residue_three_letter = yasara.NameRes(selection)[0]
        residue_one_letter = seq1(residue_three_letter)

        position = yasara.ListRes(
            f"{selection},Format=RESNUM"
        )[0]

        secondary_structure = yasara.SecStrRes(selection)[0]

        secstr_dict[position] = {
                "Chain": chain,
                "Residue": residue_one_letter,
                "SecStr": secondary_structure,
            }

    secstr_df = pd.DataFrame.from_dict(
        secstr_dict, 
        orient = "index"
        )

    secstr_df.index.name = "Pos"

    return secstr_df

#--------------------Sums up secondary structures to helices, beta-sheets and coils--------------------#

def simplify_secstr_output(secstr_df):

    secstr_dict = {}

    for position, row in secstr_df.iterrows():

        secstr = row["SecStr"]

        if secstr in {"G", "H", "I"}:
            secstr = "H"

        elif secstr == "E":
            secstr = "E"

        elif secstr in {"T", "C"}:
            secstr = "C"

        else:
            raise ValueError(
            f"Unbekannter YASARA-SecStr-Code: "
            f"{secstr!r}"
        )

        secstr_dict[position] = {
            "Residue": row["Residue"],
            "SecStr": secstr,
        }

    simplified_df = pd.DataFrame.from_dict(
        secstr_dict,
        orient="index",
    )

    simplified_df.index.name = "Pos"

    return simplified_df

#--------------------Extracts information about secondary structure of each position from pdb--------------------#

def get_secstr_scores(simple_secstr_df, output_mode):

    windows = []
    window_size = 7
    scores = {}

    secstrchain = "".join(
        simple_secstr_df["SecStr"].astype(str)
    )

    for i in range(len(secstrchain) - window_size + 1):
        window = secstrchain[i:i + window_size]
        windows.append(window)

    for i, window in enumerate(windows): 

        window_positions = list(
            simple_secstr_df.index[i:i + 7]
        )

        position = window_positions[0]

        expected_positions = list(
            range(
                position,
                position + 7,
            )
        )

        if window_positions != expected_positions:
            continue
        
        else: 
            coil_weights = [1, 1, 1, 1, 1, 1, 1]
            position_7_helix_bonus = 3

            score = 0

            for secstr, weight in zip(
                window,
                coil_weights,
            ):
                if secstr == "C":
                    score += weight

            if window[6] == "H":
                score += position_7_helix_bonus

            normalized_score = (score / 9) * 1.5
            rounded_score = round(normalized_score, 5) 

            if output_mode == min:
                scores[position] = rounded_score

            else:
                scores[position] = {
                    "SecStr Score": rounded_score, 
                    "SecStr": simple_secstr_df.loc[position, "SecStr"]
                }

        scores_df = pd.DataFrame.from_dict(
            scores,
            orient="index",
        )

    return(scores_df)


def make_combined_table(SS_per_motif, ES_per_motif, secstr_scores_df):

    final_output_dict = {}

    for motif, SS_df in SS_per_motif.items():

        if motif in ES_per_motif:

            ES_df = ES_per_motif[motif]

            combined_df = SS_df.join(
                ES_df, 
                how = "left"
            )

            o_score = (combined_df["Sequence Score"] * combined_df["Energy Score"]).round(5)

        else:
            combined_df = SS_df.copy()
            o_score = float("nan")

        combined_df = combined_df.join(
            secstr_scores_df,
            how = "left"
        )

        secstr_factor = (
            combined_df["SecStr Score"].where(
                combined_df["SecStr Score"] > 0.8,
                0.8,
            )
        )

        o_score = (
            combined_df["Sequence Score"]* combined_df["Energy Score"]* secstr_factor
        )

        o_score = (
            o_score
            .clip(upper=1.0)
            .round(5)
        )

        combined_df.insert(
            0,
            "O-Score",
            o_score
        )

        final_output_dict[motif] = combined_df
        
    return final_output_dict

def sort_output(final_tables, all_motifs, output_sorting, top_k):

    #best_per_motif
    if output_sorting == "best_per_motif":

        best_per_motif = {}

        for motif, residues_df in final_tables.items():

            efficiency = all_motifs[motif]

            valid_df = residues_df.dropna(
                subset=["O-Score"]
            )

            sorted_df = valid_df.sort_values(
                by = "O-Score",
                ascending = False,
            )

            if sorted_df.empty:
                continue

            best_per_motif[(motif, efficiency)] = sorted_df.head(top_k)

        return(best_per_motif)



def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "fasta_file", 
        help = "Enter fasta file of POI"
    )

    parser.add_argument(
        "pdb_file",
        help = "Enter pdb file of POI"
    )

    parser.add_argument(
        "--output_mode", 
        default = "default",
        choices = ["default", "min", "verbose"],
    )

    parser.add_argument(
        "--output_sorting",
        default = "best_per_motif",
        choices = ["best_per_motif", "best_overall", "best_per_positions"],
        help = "If output_sorting = best_per_positions, the positions argument becomes mandatory"
    )

    parser.add_argument(
        "--top_k",
        default = 10, 
        type = int,
    )

    parser.add_argument(
        "--positions",
        type = int, 
        nargs ="+",
        help = "Name every position, of which you want the O-score"
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    #pandas settings for data frames
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_colwidth", None)
    pd.set_option("display.width", 300)
    
    args = parse_args()

    fasta_file = args.fasta_file
    pdb = args.pdb_file
    top_k = args.top_k
    output_mode = args.output_mode
    output_sorting = args.output_sorting
    relevant_positions = args.positions
    

    fasta_seqs = get_fasta_seq(fasta)
    fasta_seq = fasta_seqs[0][1]

    all_motifs = get_motifs(csv_file)

    mismatch_dict = count_mismatches(fasta_seq, all_motifs, window_size, 7)

    default_data = BASE_DIR/"experimental_data.csv"

    exp_data = load_mutation_data(default_data)

    a_scores_per_motif = get_assessment_score(mismatch_dict, exp_data, position_scoring)
    ac_scores_per_motif = get_confidence_score(a_scores_per_motif, exp_data, position_scoring)
    SS_per_motif = make_SS_output_tables_per_motif(ac_scores_per_motif, output_mode)

    #Adjustment of pdb residue numbers to fasta

    pdb_seqs = get_pdb_seq(pdb)
    best_per_fasta = get_alignment_info(fasta_seqs, pdb_seqs)
    alignment_tables = concat_alignment_tables(best_per_fasta)
    save_alignment_tables(alignment_tables)
    prepared_pdb = reassign_pdb_residues(pdb, best_per_fasta) #Path to pdb with reassigned residues
    repair_pdb(prepared_pdb, pdb)
    pdb_for_name = Path(pdb)
    prepaired_pdb = BASE_DIR/f"{pdb_for_name.stem}_prepaired.pdb"

    #Calculation of mutation energies

    rawMutEnergyList = make_rawMutEnergyList(BASE_DIR, pyfresco_dir, prepaired_pdb, foldx_path, far_enough_res, far_enough_zone)
    MutEnergyList = get_MutEnergyList(rawMutEnergyList)

    #Energy Scoring 

    ddg_per_motif = get_ddgs(alignment_tables, all_motifs, MutEnergyList)
    energy_scores_per_motif = get_energy_scores(ddg_per_motif)
    ES_per_motif = make_ES_output_tables_per_motif(energy_scores_per_motif, output_mode)

    #Secondary Structure Scoring

    secstr_df = getSecStr(prepaired_pdb)
    simple_secstr_df = simplify_secstr_output(secstr_df)
    secstr_scores_df = get_secstr_scores(simple_secstr_df)

    #Combination of all information in one table

    final_tables = make_combined_table(SS_per_motif, ES_per_motif)

    final_sorted = sort_output(final_tables, all_motifs, output_sorting, top_k)

    for motif, motif_df in final_sorted.items():
        print(f"\nMotif: {motif}, Efficiency: {all_motifs[motif]}")
        print("-" * 45)
        print(
            motif_df.to_string(
                index=False,
                float_format=lambda value: f"{value:.4f}"
            )
        )
    


   



    

        





        
        








"""
best_ddg_df = get_ddg_of_pos(ddg_per_motif, all_motifs)
    best_score_df = get_score_of_pos(energy_scores_per_motif, all_motifs)
    best_df = best_overall(energy_scores_per_motif, all_motifs)
    best_F1_df = best_for_F1(energy_scores_per_motif, ddg_per_motif, F1)

    print(best_df)
    print(best_ddg_df)
    print(best_score_df)
    print(best_F1_df)
"""
    
"""
#--------------------Creates list of best ddg values per motifs of a given (input) position--------------------# 

def get_ddg_of_pos(ddg_per_motif, all_motifs):

    ddg_of_input_pos = []

    for motif, positions in ddg_per_motif.items():

        for residue, values in positions.items():

            if int(residue) == int(position):

                efficiency = all_motifs[motif]
                ddg_of_input_pos.append((motif, values, efficiency))
    
    ddg_of_input_pos.sort(key=lambda x: abs(x[1]["ddG"]))  
    best_ddg_df = pd.DataFrame(
        ddg_of_input_pos,
        columns=["motif", "ddG", "efficiency"]
    )

    return best_ddg_df.head(top_k)

#--------------------Creates list of best scores per motifs of a given (input) position--------------------# 

def get_score_of_pos(motif_scores, all_motifs):

    scores_of_input_pos = []

    for motif, positions in motif_scores.items():

        for residue, score in positions.items():

            if int(residue) == int(position):

                efficiency = all_motifs[motif]
                scores_of_input_pos.append((motif, score, efficiency))
                scores_of_input_pos.sort(key=lambda x: x[1], reverse = True)

    best_score_df = pd.DataFrame(
        scores_of_input_pos,
        columns=["motif", "score", "efficiency"]
    )

    return best_score_df.head(top_k)

#--------------------Creates list of best scores overall--------------------# 

def best_overall(motif_scores, all_motifs):

    best_all = []

    for motif, positions in motif_scores.items():

        for residue, score in positions.items():

            efficiency = all_motifs[motif]
            best_all.append((motif, score, efficiency, residue))
    
    best_all.sort(key=lambda x: x[1], reverse = True)
    best_df = pd.DataFrame(
        best_all,
        columns=["motif", "score", "efficiency", "residue"]
    )
    
    return best_df.head(top_k)

#--------------------Creates list of best scores for F1--------------------# 

def best_for_F1(motif_scores, ddg_per_motif, F1):

    best_energy_F1 = []

    for motif, positions, in motif_scores.items():

        if motif == F1:

            for residue, score in positions.items():

                best_energy_F1.append((score, residue))
                best_energy_F1.sort(key=lambda x: x[0], reverse = True)

    best_F1 = []

    for entry in best_energy_F1:

        for motif2, positions in ddg_per_motif.items():
        
                if motif2 == F1:

                    for position, values in positions.items():

                        if position == entry[1]:

                            best_F1.append((entry[0], entry[1], values))

    best_F1_df = pd.DataFrame(
        best_F1,
        columns=["score", "residue", "ddg + sd"]
    )

    return best_F1_df.head(top_k)

"""



