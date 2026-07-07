from pathlib import Path
from Bio import Align
import pandas as pd
import numpy as np
import subprocess
import csv
import yasara

work_dir = Path("/Users/yannikmeindl/FoldX/FoldX_exp")

#get_fasta_sequence
fasta = work_dir / "rcsb_pdb_4XHF.fasta"

#get individual_list
SS_output = work_dir/"sequence_scorer_output.txt"
motif = "DGLSGAT"

#repairpdb
foldx = Path("/Users/yannikmeindl/FoldX/Programm/FoldX")
pdb = work_dir/"4XHF.pdb"

#buildmodel
repaired_pdb = work_dir/"Shewa_NqrC_model_0.pdb"
mutant_file = work_dir/"individual_list.txt"

#exp_data
csv_file = work_dir/"experimental_data.csv"

#read_mutation_data
mutation_results = work_dir/"BuildModel_test2"
avg_list = next(mutation_results.glob("Average*"))


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
    yasara.OligomerizeObj("1", center="Yes", instance="No")

    pdb_seq = work_dir / "pdb_sequences.txt"

    yasara.SaveSeqMol("All", pdb_seq, join="No")

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

fasta_seqs = get_fasta_seq(fasta)
pdb_seqs = get_pdb_seq(pdb)
best_per_fasta = get_alignment_info(fasta_seqs, pdb_seqs)
alignment_tables = concat_alignment_tables(best_per_fasta)

print(alignment_tables)










        

    

#get_individual_list parses through SS_output and aligns the sequences of the hits with the motif
#then, the information about the necessary mutation and its position is brought into the correct format for FoldX
def get_individual_list(SS_output, motif):

    mut_list = []

    with SS_output.open() as f:

        for line in f:
            line = line.strip()

            if line.startswith("Pos") or line.startswith("_"):
                continue

            elements = line.split(maxsplit=5)

            pos = int(elements[0])
            seq = elements[1]
            mm = int(elements[2])
            assess = float(elements[3])
            confi = float(elements[4])
            details = elements[5]

            pp_mutlist = []

            for i, (a, b) in enumerate(zip(motif, seq)):

                if b!=a:
                    pp_mutlist.append(f"{b}A{pos+i}{a}")

                mut_line = ",".join(pp_mutlist) + ";"
            
            mut_list.append(mut_line)

            individual_list = "\n".join(mut_list)

    return(individual_list)

#output: repaired structure (is automatically saved in working directory)
def repair_pdb(foldx, pdb, work_dir):

    result = subprocess.run(
        [
            str(foldx),
            "--command=RepairPDB",
            f"--pdb={pdb.name}",
            f"--pdb-dir=./",
            f"--output-dir=./",
        ],
        cwd = work_dir, 
        capture_output = True, 
        text = True,
    )

#make output folder variable (by option in argparse)
#Attention: FoldX doesnt create a output Folder => Folder needs to be made before by hand
#models and lists again automatically saved in working directory 
def get_model(repaired_pdb, mutant_file, work_dir):

    mutants = subprocess.run(
        [
            str(foldx),
            "--command=BuildModel",
            f"--pdb={repaired_pdb.name}",
            f"--mutant-file={mutant_file}",
            f"--pdb-dir=./",
            f"--output-dir={work_dir}/BuildModel_test2",
            f"--numberOfRuns=5",
        ],
        cwd = work_dir, 
        capture_output = True, 
        text = True,
    )

    return mutants

def read_mutant_energy(avg_list):
    
    with open(avg_list, "r") as f:
        lines = f.readlines()
    
    header_line = None

    for i, line in enumerate(lines):
        
        if line.startswith("Pdb\t"):
            header_line = i
            break

    if header_line == None:
        raise ValueError(f"File doesn't contain mutant energies!")

    energy_list = pd.read_csv(avg_list, sep="\t", skiprows = header_line)

    return energy_list


def open_exp_data(csv_file):
    
    exp_data = pd.read_csv(csv_file, sep = ";", decimal = ",")
    
    return exp_data





#individual_list = (get_individual_list(SS_output, motif))
#output_file = work_dir / "individual_list.txt"
#output_file.write_text(individual_list, encoding="utf-8")


#repair_pdb(foldx, pdb, work_dir)

#result = get_model(repaired_pdb, mutant_file, work_dir)

#energy_list = (read_mutant_energy(avg_list))

#energy_output = pd.DataFrame(energy_list)

#energy_output.to_csv(work_dir / "energy_output.csv")

#print(energy_list[["Pdb", "total energy", "SD"]])




#exp_data = open_exp_data(csv_file)

#print(exp_data)



#print(load_mutation_data(csv_file))

#print("Returncode:", result.returncode)
#print("STDOUT:")
#print(result.stdout)
#print("STDERR:")
#print(result.stderr)








