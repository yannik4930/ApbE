from pathlib import Path
from Bio import Align
from Bio.SeqUtils import seq1
from datetime import datetime
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
foldx_path = "/scratch/s6765211/FoldX/foldx5_Linux/foldx_20270131"

#non-mutable region
far_enough_res = None
far_enough_zone = None


#open yasara in text mode (won't work on cluster otherwise)
yasara.info.mode = "txt"
yasara.run("Console Off")

recognition_motifs = {
    "F1": "GVDGLSGATLTS",
    "F2": "DGLSGAT",
}

f2_motif = recognition_motifs["F2"]


def define_recognition_motif(input_motif):
    return recognition_motifs[input_motif]


def get_f2_position(motif_position, recognition_motif):

    f2_offset = recognition_motif.index(f2_motif)
    f2_position = motif_position - f2_offset

    if 1 <= f2_position <= len(f2_motif):
        return f2_position

    return None


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

def get_motifs(csv_file, recognition_motif):

    exp_data = pd.read_csv(
        csv_file,
        sep=";",
        decimal=",",
    )

    all_motifs = {
        recognition_motif: 1.0
    }

    f2_offset = recognition_motif.index(f2_motif)

    for _, row in exp_data.iterrows():

        f2_position = int(row["Position"]) 

        if f2_position == 6:
            continue

        motif_index = f2_offset + f2_position 

        mutated_motif = list(recognition_motif)
        mutated_motif[motif_index] = row["Mutation"]
        mutated_motif = "".join(mutated_motif)

        if mutated_motif == recognition_motif:
            print("Referenzmotiv wird überschrieben:")
            print(row)

        all_motifs[mutated_motif] = row["Efficiency"]

    return all_motifs

# sliding window mechanism splits sequence into 7mers (windows)
# creates dictionary with a list of windows for each number of mismatches
def count_mismatches(sequence, all_motifs, window_size):

    windows = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        windows.append(window) 

    mismatch_dict = {}

    for motif in all_motifs:

        mismatches_per_motif = {
            mismatch_count: []
            for mismatch_count in range(window_size + 1)
        }
    
        for index, window in enumerate(windows):
            mismatches = 0
            mismatch_info = []

            for pos, (a, b) in enumerate(zip(window, motif)):
                if a != b:
                    mismatches += 1
                    mismatch_info.append(((pos + 1), a, b))
        
            mismatches_per_motif[mismatches].append((index + 1, window, mismatch_info))

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
def get_assessment_score(mismatch_dict, exp_data, position_scoring, recognition_motif):

    a_scores_per_motif = {}

    for motif, mismatches_per_motif in mismatch_dict.items():

        scored_motifs = {}
        motif_length = len(motif)

        for mismatch_count, motif_hits in mismatches_per_motif.items():

            for pos, window_sequence, mismatch_details in motif_hits:

                assessment_score = motif_length - mismatch_count

                for motif_position, window_aa, motif_aa in mismatch_details:

                    f2_position = get_f2_position(
                        motif_position,
                        recognition_motif,
                    )

                    if f2_position is None:
                        continue

                    if f2_position == 7:
                        continue

                    mismatch = (
                        f2_position,
                        window_aa,
                        motif_aa,
                    )

                    if mismatch in exp_data:

                        efficiency = exp_data[mismatch]
                        assessment_score += efficiency

                    else:

                        scoring_data = position_scoring[f2_position]
                        matrix = scoring_data["matrix"]

                        matrix_value = matrix[
                            (window_aa, motif_aa)
                        ]

                        matrix_based_efficiency = get_efficiency(
                            matrix_value,
                            scoring_data["a"],
                            scoring_data["b"],
                            scoring_data["R2"],
                        )

                        capped_efficiency = min(
                            matrix_based_efficiency,
                            1.0,
                        )

                        assessment_score += capped_efficiency

                normalized_assessment_score = (
                    assessment_score / motif_length
                )

                rounded_assessment_score = round(
                    normalized_assessment_score,
                    5,
                )

                scored_motifs[pos] = {
                    "Sequence Score": rounded_assessment_score,
                    "MM": len(mismatch_details),
                    "Window Sequence": window_sequence,
                    "Mismatch Details": mismatch_details
                }

        a_scores_per_motif[motif] = scored_motifs

    return a_scores_per_motif


def make_SS_output_tables_per_motif(a_scores_per_motif, output_mode):

    SS_output_dict = {}

    for motif, scores_per_motif in a_scores_per_motif.items():

        motif_scores = (
            pd.DataFrame.from_dict(
                scores_per_motif,
                orient = "index"
            )
            .rename_axis("Pos")
        )

        if output_mode == "verbose":

            output_df = motif_scores

        elif output_mode == "min":

            output_df = motif_scores[
                ["Sequence Score", "MM"]
            ]
            
        else:
            
            output_df = motif_scores[
                ["Sequence Score", "MM", "Window Sequence"]
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

    pdb_seq = output_dir / "pdb_sequences.txt" #variable for the file where yasara saves the sequences in the next line

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

        seq_name = fasta_id.lstrip(">").split("_")[0]

        output_file = output_dir / f"{seq_name}_align_table.tab"

        alignment_table.to_csv(output_file, sep = "\t")

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

    resultlist = yasara.ListMol("Protein Obj 1") #creates list with all molecules of obj 1 
    for i, a in enumerate(resultlist, start = 1): # assigns index to list of molecules 
        yasara.NameMol(f"Obj 1 Mol {a}, {i}") #replaces old mol name with index

    selector = yasara_selector_from_fasta_header(best_per_fasta) 
    first_res_of_pdb_seq = list(best_per_fasta.values())[0]["aligned_areas"][0][0][0] + 1 #takes first residue of pdb compared to fasta

    yasara.NumberRes(selector, first = first_res_of_pdb_seq) #Reassigns residue numbers of pdb based on numbers from alignment
    prepared_pdb = output_dir / f"{pdb.stem}_aligned.pdb" #variable for place to save the prepared pdb
    yasara.SavePDB(
        "Obj 1", 
        prepared_pdb,
        format="PDB3"
    )
    
    return(prepared_pdb)

#--------------------Fix protonation states--------------------#

def repair_pdb(prepared_pdb, pdb):

    pdb = Path(pdb)

    yasara.Clear()
    yasara.LoadPDB(prepared_pdb)
    yasara.CleanAll()
    yasara.OptHydAll(method = "Yasara")
    prepaired_pdb = output_dir/f"{pdb.stem}_prepared.pdb"
    yasara.SavePDB(
        "Obj 1",
        prepaired_pdb,
        format="PDB3"
    )
    os.remove(prepared_pdb)

#--------------------Execute pyFRESCO to get ddG of all possible mutations and prepares resulting rawMutEnergyList--------------------#

def make_rawMutEnergyList(work_dir, output_dir, prepaired_pdb, foldx_path, far_enough_res, far_enough_zone, reuse_mutational_data, pdb):

    energy_dir = create_energy_output_folder(output_dir)

    final_mutation_list = (energy_dir / "MutationEnergies_CompleteList.tab")

    if reuse_mutational_data:

        if (
            final_mutation_list.is_file()
            and final_mutation_list.stat().st_size > 0
        ):
            print(
                "Existing mutation energies will be reused."
            )
            return final_mutation_list

        raise FileNotFoundError(
            "Mutational data should be reused, but no "
            "valid MutationEnergies_CompleteList.tab "
            f"was found in {energy_dir}."
        )

    selection_tab = energy_dir / f"{prepaired_pdb.stem}.tab"

    shutil.copy(work_dir / "DistributeFoldx.py", energy_dir / "DistributeFoldx.py")
    shutil.copy(work_dir / "submit_fresco.sh", energy_dir / "submit_fresco.sh")
    shutil.copy(work_dir / "fresco_job.sh", energy_dir / "fresco_job.sh")

    shutil.copy(output_dir / f"{pdb.stem}_prepaired.pdb", energy_dir / f"{pdb.stem}_prepaired.pdb")

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
    ], cwd = energy_dir, check = True)

    subprocess.run(
        ["bash", 
        "submit_fresco.sh"
    ], cwd = energy_dir, check = True)

    subprocess.run([
        sys.executable, 
        "DistributeFoldx.py", 
        "Phase2",
        prepaired_pdb.name,
        "-5"
    ], cwd = energy_dir, check = True)

    if (
        not final_mutation_list.is_file()
        or final_mutation_list.stat().st_size == 0
    ):
        raise RuntimeError(
            "FoldX Phase2 did not create a valid "
            "MutationEnergies_CompleteList.tab"
        )

    """
    for path in output_dir.iterdir():

            directory_number = path.name.removeprefix(
                "Subdirectory"
            )

            if (
                path.is_dir()
                and path.name.startswith("Subdirectory")
                and directory_number.isdigit()
            ):
                shutil.rmtree(path)
                print(f"Deleted {path.name}")
    """
    
    return energy_dir / "MutationEnergies_CompleteList.tab"


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
                        ddg_per_mutation.setdefault(pos, []).append({
                            "mutation": "None",
                            "ddG": 0,
                            "sd": 0,
                        })
            
            #sums up the ddg + sd values of every 7mer and puts them into a dictionary
            #=> one dictionary for every motif
            for pos, mutations in ddg_per_mutation.items():
                
                #pos+1 so that output shows biological residue number instead of position in python string
                ddg_per_position[pos+1] = {
                    "ddG": round(sum(abs(mut["ddG"]) for mut in mutations), 4),
                    "sd": round(sum(mut["sd"] for mut in mutations), 4),
                }

            #puts the dictionaries of all the motifs into one dictionary 
            ddg_per_motif[motif] = ddg_per_position

    with open(output_dir / "ddg_per_motif.json", "w") as file:
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
            if 0 < ddg_value <= 20:
                score = 1 - 0.01 * ddg_value

            elif ddg_value > 20:
                score = 0.5 + 0.3 * math.exp(-0.04 * (ddg_value - 20))

            elif 0 > ddg_value >= -20: 
                score = 1 + 0.01 * ddg_value

            elif ddg_value < -20: 
                score = 0.5 + 0.3 * math.exp(0.04 * (ddg_value + 20))

            else:
                score = 1

            score = round(score, 5)

            residue_scores[residue] = {
                "Energy Score": score, 
                "ddG": ddg_value,
                "sd": ddg_sd 
            }

        motif_scores[motif] = residue_scores    

    with open(output_dir / "motif_scores.json", "w") as file:
        json.dump(motif_scores, file, indent = 4)
    
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

def getSecStr(prepaired_pdb, fasta_seqs, pdb_seq):

    yasara.Clear()
    yasara.LoadPDB(prepaired_pdb)

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
    
    prepaired_pdb_best_per_fasta = get_alignment_info(fasta_seqs, records)
    selector = yasara_selector_from_fasta_header(prepaired_pdb_best_per_fasta) 
    
    ca_atoms = yasara.ListAtom(f"CA Protein {selector}")

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

def get_secstr_scores(simple_secstr_df, output_mode, recognition_motif):

    f2_offset = recognition_motif.index(f2_motif)
    essential_thr_index = f2_offset + 6
    beginning_of_allowed_aH = essential_thr_index - 1
    end_of_allowed_bS = essential_thr_index - 5

    windows = []
    scores = {}

    secstrchain = "".join(simple_secstr_df["SecStr"].astype(str))

    for i in range(len(secstrchain) - window_size + 1):
        window = secstrchain[i:i + window_size]
        windows.append(window)

    for i, window in enumerate(windows): 

        window_positions = list(simple_secstr_df.index[i:i + window_size])

        position = window_positions[0]

        expected_positions = list(
            range(
                position,
                position + window_size,
            )
        )

        if window_positions != expected_positions:
            continue
        
        else: 
            coil_weights = [1] * window_size
            position_7_helix_bonus = 3

            score = 0

            for secstr, weight in zip(window, coil_weights,):

                if secstr == "C":
                    score += weight

            for index in range(end_of_allowed_bS + 1):

                if window[index] == "E":
                    score += 1

            for index in range(beginning_of_allowed_aH, len(window)):

                if (window[index] == "H" and index != essential_thr_index):
                    score += 1

            if window[essential_thr_index] == "H":
                score += position_7_helix_bonus

            maximum_score = (window_size - 1 + position_7_helix_bonus)
            normalized_score = (score / maximum_score) * 1.5
            rounded_score = round(normalized_score, 5) 

            if output_mode == min:
                scores[position] = {
                    "SecStr Score" : rounded_score
                }

            else:
                scores[position] = {
                    "SecStr Score": rounded_score, 
                    "SecStr": window,
                }

    scores_df = pd.DataFrame.from_dict(
        scores,
        orient="index",
    )

    return(scores_df)


#________________________________VOLUME SCORER__________________________________#
#--------------------Calculates occupancy for volumes in front of each residues sidechain--------------------#

def get_volumes(prepaired_pdb, fasta_seqs, pdb_seq):

    yasara.Clear()
    yasara.LoadPDB(prepaired_pdb)

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
    
    prepaired_pdb_best_per_fasta = get_alignment_info(fasta_seqs, records)
    selector = yasara_selector_from_fasta_header(prepaired_pdb_best_per_fasta)

    parts = selector.split()
    obj_num = int(parts[1])

    resnumber_list = yasara.ListRes(f"Protein {selector}, Format=RESNUM")
    residue_list = yasara.ListRes(f"Protein {selector}, Format=RESName")

    residue_dict = {
        int(resnumber): residue_name
        for resnumber, residue_name in zip(
            resnumber_list,
            residue_list,
            strict=True,
        )
    }   

    volume_dict = {}

    for pos, residue in residue_dict.items():    

        print(residue)

        if residue in ["Pro", "Gly"]:

            yasara.SwapRes(pos, "Ala", isomer = "L")

            CA_selection = f"{selector} Res {pos} Atom CA"
            CB_selection = f"{selector} Res {pos} Atom CB"

            CA_atoms = yasara.ListAtom(CA_selection)
            CB_atoms = yasara.ListAtom(CB_selection)

            CA_atom = CA_atoms[0]
            CB_atom = CB_atoms[0]

            CA = np.array(yasara.PosAtom(CA_atom, coordsys = "global"), dtype = float)
            CB = np.array(yasara.PosAtom(CB_atom, coordsys = "global"), dtype = float)


            yasara.SwapRes(pos, residue, isomer = "L")

        else:

            CA_selection = f"{selector} Res {pos} Atom CA"
            CB_selection = f"{selector} Res {pos} Atom CB"

            CA_atoms = yasara.ListAtom(CA_selection)
            CB_atoms = yasara.ListAtom(CB_selection)

            CA_atom = CA_atoms[0]
            CB_atom = CB_atoms[0]

            CA = np.array(yasara.PosAtom(CA_atom, coordsys = "global"), dtype = float)
            CB = np.array(yasara.PosAtom(CB_atom, coordsys = "global"), dtype = float)


        BA = CB - CA
        BA_length = np.linalg.norm(BA)
        unit_BA = BA / BA_length

        C = CB + unit_BA * 9
        R = 9
        
        x, y, z = C

        dummy_obj = yasara.BuildAtom("C")

        yasara.PosObj(dummy_obj, x=x, y=y, z=z)

        sphere_selection = (
            f"Protein with distance<{R} from Obj {dummy_obj}"
        )

        atoms_in_sphere = yasara.ListAtom(sphere_selection)

        atom_volumes = yasara.VolumeAtom(atoms_in_sphere, Type = "VdW")[0]

        atom_volumes = round(atom_volumes, 5)

        volume_dict[pos] = atom_volumes

        yasara.DelObj("2")

    return volume_dict


#--------------------Calculates scores out of the occupancy rates--------------------#

def score_volumes(volume_dict):

    pi = math.pi
    sphere_radius = float(9)
    sphere_volume = 4/3*pi*sphere_radius**3
    sphere_volume = round(sphere_volume, 5)

    volume_scores = {}

    for pos, atom_volume in volume_dict.items():

        occupancy = atom_volume / sphere_volume

        score = 1 - occupancy
        score = round(score, 5)

        volume_scores[pos] = {
            "Volume Score": score,
            "Occupied Volume": f"{round(atom_volume, 0)}/{round(sphere_volume, 0)} Å\u00B3"
        }


    volumes_df = pd.DataFrame.from_dict(
        volume_scores,
        orient="index",
    )    

    return volumes_df

#________________________________OUTPUT MANAGEMENT__________________________________#

def create_output_folder(pdb, input_run_name = None):

    pdb_id = pdb.stem

    if input_run_name is None:
        run_name = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    else:
        run_name = input_run_name.strip().replace(" ", "_")

    output_dir = (BASE_DIR / "outputs" / pdb_id / run_name)

    output_dir.mkdir(parents = True, exist_ok = True)

    return output_dir

def create_energy_output_folder(output_dir):

    energy_dir = (output_dir / "pyFRESCO")
    energy_dir.mkdir(parents = True, exist_ok = True)

    return energy_dir

def get_final_output_path(output_dir, reuse_mutational_data, pdb, input_run_name):

    if not reuse_mutational_data:
        return output_dir / f"{pdb.stem}_{input_run_name}_final.tsv"

    reuse_number = 2

    while True:

        output_path = (output_dir / f"{pdb.stem}_{input_run_name}_final{reuse_number}.tsv")

        if not output_path.exists():
            return output_path

        reuse_number += 1

def make_combined_table(scorer_results, all_motifs, input_motif, final_output_path):

    sequence_tables = scorer_results.get("sequence")
    energy_tables = scorer_results.get("energy")
    secstr_df = scorer_results.get("secstr")
    volumes_df = scorer_results.get("volume")

    final_tables = {}

    if volumes_df is not None:

        volume_offset = {
            "F2": 6,  
            "F1": 8,  
        }[input_motif]

        volumes_aligned = volumes_df.copy()
        volumes_aligned.index = (
            volumes_aligned.index - volume_offset
        )

    else:
        volumes_aligned = None

    for motif in all_motifs:

        dataframes = []

        if (
            sequence_tables is not None
            and motif in sequence_tables
        ):
            dataframes.append(sequence_tables[motif])
    
        if (
            energy_tables is not None
            and motif in energy_tables
        ):
            dataframes.append(energy_tables[motif])
            
        if secstr_df is not None:

            dataframes.append(secstr_df)
            
        if volumes_df is not None:

            dataframes.append(volumes_aligned)

        if not dataframes:
            continue

        combined_df = pd.concat(
            dataframes,
            axis=1,
            join="outer",
        )

        combined_df = combined_df.loc[combined_df.index >= 1].copy()

        o_score = pd.Series(1.0, index=combined_df.index, dtype=float)
        maximum_o_score = 1.0

        if (
            sequence_tables is not None
            and "Sequence Score" in combined_df.columns
        ):
            o_score *= combined_df["Sequence Score"]

        if (
            energy_tables is not None
            and "Energy Score" in combined_df.columns
        ):
            o_score *= combined_df["Energy Score"]

        if (
            secstr_df is not None
            and "SecStr Score" in combined_df.columns
        ):
            secstr_factor = combined_df["SecStr Score"].apply(
                lambda value: (
                    1.5
                    if value > 0.8
                    else 0.8
                    if pd.notna(value)
                    else np.nan
                )
            )

            o_score *= secstr_factor
            maximum_o_score *= 1.5

        if (
            volumes_aligned is not None
            and "Volume Score" in combined_df.columns
        ):

            o_score *= combined_df["Volume Score"]


        combined_df.insert(
            0,
            "O-Score", 
            (
                o_score
                .div(maximum_o_score)
                .clip(lower=0.0, upper=1.0)
                .round(5)
            ),
        )

        combined_df.index.name = "Pos"

        if "MM" in combined_df.columns:
            combined_df["MM"] = combined_df["MM"].astype("Int64")

        columns_to_end = [
            "Window Sequence",
            "Mismatch Details",
        ]

        for column in columns_to_end:
            if column in combined_df.columns:
                values = combined_df.pop(column)
                combined_df[column] = values

        efficiency = all_motifs[motif]

        final_tables[(motif, efficiency)] = combined_df

    output_tables = []

    for (motif, efficiency), df in final_tables.items():

        output_df = df.copy()
        output_df.insert(0, "Efficiency", efficiency)
        output_df.insert(0, "Motif", motif)

        output_tables.append(output_df)
    
    if output_tables:
        best_per_motif_df = pd.concat(output_tables)

        best_per_motif_df.to_csv(
            final_output_path,
            sep="\t",
            index=True,
            index_label="Pos",
        )

    return final_tables
       
def sort_output(final_tables, all_motifs, top_k):
    
    best_per_motif = {}

    for motif, residues_df in final_tables.items():

        valid_df = residues_df.dropna(
            subset=["O-Score"]
        )

        sorted_df = valid_df.sort_values(
            by = "O-Score",
            ascending = False,
        )

        if sorted_df.empty:
            continue

        best_per_motif[(motif)] = sorted_df.head(top_k)

    return best_per_motif
    

#________________________________SCORING FUNCTIONS__________________________________#

def run_sequence_scorer(fasta_seq, all_motifs, window_size, exp_data, recognition_motif, output_mode):

    mismatch_dict = count_mismatches(fasta_seq, all_motifs, window_size)
    a_scores_per_motif = get_assessment_score(mismatch_dict, exp_data, position_scoring, recognition_motif)
    SS_per_motif = make_SS_output_tables_per_motif(a_scores_per_motif, output_mode)

    return SS_per_motif

def run_energy_scorer(alignment_tables, all_motifs, prepaired_pdb, output_dir, output_mode, reuse_mutational_data):

    rawMutEnergyList = make_rawMutEnergyList(BASE_DIR, output_dir, prepaired_pdb, foldx_path, far_enough_res, far_enough_zone, reuse_mutational_data, pdb)
    MutEnergyList = get_MutEnergyList(rawMutEnergyList)
    ddg_per_motif = get_ddgs(alignment_tables, all_motifs, MutEnergyList)
    energy_scores_per_motif = get_energy_scores(ddg_per_motif)
    ES_per_motif = make_ES_output_tables_per_motif(energy_scores_per_motif, output_mode)

    return ES_per_motif

def run_secstr_scorer(prepaired_pdb, fasta_seqs, pdb_seq, output_mode, recognition_motif):

    secstr_df = getSecStr(prepaired_pdb, fasta_seqs, pdb_seq)
    simple_secstr_df = simplify_secstr_output(secstr_df)
    secstr_scores_df = get_secstr_scores(simple_secstr_df, output_mode, recognition_motif)

    return secstr_scores_df

def run_volume_scorer(prepaired_pdb, fasta_seqs, pdb_seq):

    volumes_dict = get_volumes(prepaired_pdb, fasta_seqs, pdb_seq)
    volumes_df = score_volumes(volumes_dict)

    return volumes_df

#________________________________INPUT__________________________________#

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "motif", 
        default = "F1",
        choices = ["F1", "F2"],
        help = "Specify which motif you want the code to run on. \nOptions: F1, F2"
    )

    parser.add_argument(
        "fasta_file", 
        help = "Enter fasta file of POI"
    )

    parser.add_argument(
        "pdb_file",
        help = "Enter pdb file of POI"
    )

    parser.add_argument(
        "--run_name",
        help= "Specify name of run"        
    )

    parser.add_argument(
        "--reuse_mutational_data", 
        action = "store_true",  
        help = "Pass this flag to reuse the mutational data from a previous run.\nImportant: Specify the run name whose data should be reused."
    )

    parser.add_argument(
        "--scorers",
        nargs="+",
        choices=SCORERS,
        default=SCORERS,
        help="Scorers that should be applied. Options: 'sequence', 'energy', 'secstr' 'volume'. Without further specification all are applied.",
    )

    parser.add_argument(
        "--output_mode", 
        default = "default",
        choices = ["default", "min", "verbose"],
        help = "Choose how detailed the output table is supposed to be.\nOptions: default, min, verbose"
    )

    parser.add_argument(
        "--top_k",
        default = 10, 
        type = int,
    )

    args = parser.parse_args()

    return args


#________________________________MAIN__________________________________#

if __name__ == "__main__":

    SCORERS = (
        "sequence",
        "energy",
        "secstr",
        "volume",
    )

    #Arguments
    
    args = parse_args()
    input_motif = args.motif
    fasta = Path(args.fasta_file)
    pdb = Path(args.pdb_file)
    top_k = args.top_k
    output_mode = args.output_mode
    input_run_name = args.run_name
    selected_scorers = args.scorers
    reuse_mutational_data = args.reuse_mutational_data

    #Pandas settings

    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_colwidth", None)
    pd.set_option("display.width", 1000)

    #General

    recognition_motif = define_recognition_motif(input_motif)
    window_size = len(recognition_motif)
    output_dir = create_output_folder(pdb, input_run_name)
    final_output_path = get_final_output_path(output_dir, reuse_mutational_data, pdb, input_run_name)

    #Fasta sequence preparation

    fasta_seqs = get_fasta_seq(fasta)
    fasta_seq = fasta_seqs[0][1]

    #Preparation of exprimental data and tested motifs

    default_data = BASE_DIR/"experimental_data.csv"
    exp_data = load_mutation_data(default_data)
    all_motifs = get_motifs(default_data, recognition_motif)

    #Preparation and repairation of PDB + adjustment to fasta sequence
    
    pdb_seqs = get_pdb_seq(pdb)
    pdb_seq = output_dir / "pdb_sequences.txt"
    best_per_fasta = get_alignment_info(fasta_seqs, pdb_seqs)
    alignment_tables = concat_alignment_tables(best_per_fasta)
    save_alignment_tables(alignment_tables)
    prepared_pdb = reassign_pdb_residues(pdb, best_per_fasta) #Path to pdb with reassigned residues
    repair_pdb(prepared_pdb, pdb)
    pdb_for_name = Path(pdb)
    prepaired_pdb = output_dir/f"{pdb_for_name.stem}_prepaired.pdb"

    #Modular scoring

    scorer_results = {}

    if "sequence" in selected_scorers:
        scorer_results["sequence"] = run_sequence_scorer(fasta_seq, all_motifs, window_size, exp_data, recognition_motif, output_mode)

    if "energy" in selected_scorers:
        scorer_results["energy"] = run_energy_scorer(alignment_tables, all_motifs, prepaired_pdb, output_dir, output_mode, reuse_mutational_data)

    if "secstr" in selected_scorers:
        scorer_results["secstr"] = run_secstr_scorer(prepaired_pdb, fasta_seqs, pdb_seq, output_mode, recognition_motif)

    if "volume" in selected_scorers:
        scorer_results["volume"] = run_volume_scorer(prepaired_pdb, fasta_seqs, pdb_seq)

    #Combination of all information in one table

    final_tables = make_combined_table(scorer_results, all_motifs, input_motif, final_output_path)
    final_sorted = sort_output(final_tables, all_motifs, top_k)

    for motif, motif_df in final_sorted.items():
        print(f"\nMotif: {motif[0]}, Efficiency: {motif[1]}")
        print("-" * 45)
        print(
            motif_df.to_string(
                index=True,
                float_format=lambda value: f"{value:.4f}"
            )
        )


