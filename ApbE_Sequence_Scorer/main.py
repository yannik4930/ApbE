import csv
import argparse
from pathlib import Path

motif = "DGLSGAT"
window_size = len(motif)


def read_fasta(fasta_content):
   
    if fasta_content[0].startswith(">"):
        pass
    else:
        print("No FASTA file found!\nMake sure file has the following format:\n>POI name\nPOI sequence")

    sequence = ""

    for line in fasta_content[1:]:

        if line.startswith(">"):
            print("Textfile contains more than one FASTA file. \nMake sure that textfile contains only one FASTAfile of the following format:\n>POI name\nPOI sequence")
            raise ValueError()  
        else:
            pass

        sequence += line.strip()

    return sequence


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


# sliding window mechanism splits sequence into 7mers (windows)
# creates dictionary with a list of windows for each number of mismatches
def count_mismatches(sequence, motif, window_size, max_mismatches: str):
    windows = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        windows.append(window) 

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
        pos_mismatch = []

        for pos, (a, b) in enumerate(zip(window, motif)):
            if a != b:
                mismatches += 1
                pos_mismatch.append(((pos + 1), a, b))
    
        if mismatches <= max_mismatches:
            mismatches_per_motif[mismatches].append(((index + 1), window, pos_mismatch))

    return mismatches_per_motif


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
    scored_motifs = [] 
   
    for mismatch_count, motif_hits in mismatch_dict.items():
       
        for pos, window_sequence, mismatch_details in motif_hits:
            assessment_score = 0 
            assessment_score += window_size - mismatch_count
        
            for mm_pos, original_aa, window_aa in mismatch_details:
                mismatch = (mm_pos, original_aa, window_aa)

                if mismatch in exp_data:
                    efficiency = exp_data[mismatch]
                    assessment_score += efficiency

                else: 
                    if mm_pos == 7:
                        assessment_score += 0
                    
                    else:
                        best_matrix = position_scoring[mm_pos]["matrix"]
                        matrix_val = best_matrix[(original_aa, window_aa)]
                        matrix_based_efficiency = get_efficiency(matrix_val, position_scoring[mm_pos]["a"], position_scoring[mm_pos]["b"], position_scoring[mm_pos]["R2"])
                        capped_efficiency = min(matrix_based_efficiency, 1) 
                        assessment_score += capped_efficiency
                    
            scored_motifs.append((pos, window_sequence, assessment_score, mismatch_details))
        
    scored_motifs.sort(key=lambda motif: motif[2], reverse=True)

    return(scored_motifs)


def get_confidence_score(scored_motifs, exp_data, position_scoring):
    scored_motifs_2 = [] 
       
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

        scored_motifs_2.append((pos, window_sequence, len(mismatch_details), assessment_score, confidence_score, mismatch_details))
        
    return(scored_motifs_2)


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "fasta_file", 
        help = "Enter fasta file of POI"
    )
    parser.add_argument(
        "--output", 
        default = "default",
        choices = ["default", "min", "verbose"],
    )
    parser.add_argument(
        "--top_k",
        default = 10, 
        type = int,
    )
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_args()

    fasta_file = args.fasta_file
    top_k = args.top_k
    output_mode = args.output

    with open(fasta_file, "r", encoding="utf-8") as f: 
        fasta_content = f.readlines()   

    sequence = read_fasta(fasta_content)

    mismatch_dict = count_mismatches(sequence, motif, window_size, 7)

    BASE_DIR = Path(__file__).resolve().parent
    default_data = BASE_DIR/"experimental_data.csv"

    exp_data = load_mutation_data(default_data)

    scored_motifs = get_assessment_score(mismatch_dict, exp_data, position_scoring)

    scored_motifs_2 = get_confidence_score(scored_motifs, exp_data, position_scoring)

#don't change spaces and positions (=> will mess up formating of the output table)
    if output_mode == "verbose":

        print(f"{"Pos":<5} {"Sequence":<9} {"MM":>4} {"Assessment":>12} {"Confidence":>12}   {"Details":<12}")
        print("_" * 64)

        for motif in scored_motifs_2[:top_k]:
            print(f"{motif[0]:>4}  {motif[1]:<9} {motif[2]:>4} {f"{motif[3]:.5f}":>12} {f"{motif[4]:.5f}":>12}   {f"{motif[5]}":<}")

    elif output_mode == "min":

        print(f"{"Pos":<5}  {"Assessment":<12} {"Confidence":<12}")
        print("_" * 33)

        for motif in scored_motifs_2[:top_k]:
            print(f"{motif[0]:>4}  {f"{motif[2]:.5f}":>11} {f"{motif[3]:.5f}":>12}")
        
    else:
         
        print(f"{"Pos":<5} {"Sequence":<9} {"MM":>4} {"Assessment":>12} {"Confidence":>12}")
        print("_" * 46)

        for motif in scored_motifs_2[:top_k]:
            print(f"{motif[0]:>4}  {motif[1]:<9} {motif[2]:>4} {f"{motif[3]:.5f}":>12} {f"{motif[4]:.5f}":>12}")



        



    

                    











