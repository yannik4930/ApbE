### imports mismatch_function from project folder to provide test dictionary for scorer function###

import csv 

from mismatch_function import count_mismatches
from open_csv import load_mutation_data
from open_matrix import load_substitution_matrix

sequence = "DGLSGATDKLSGATGEDGISGETDGRSGAT"
motif = "DGLSGAT"
window_size = len(motif)

mismatch_dict = count_mismatches(sequence, motif, window_size, 7)
exp_data = load_mutation_data("experimental_data.csv")

position_scoring = {
    0: {
        "matrix": load_substitution_matrix("matrices/BLOSUM62.txt"), 
        "a": 10.505,
        "b": -7.8119,
        "R2": 0.5876
    }, 
    1: {
        "matrix": load_substitution_matrix("matrices/GRANTHAM.txt"), 
        "a": -362.82,
        "b": 384.49,
        "R2": 0.9424
    },
    2: {
        "matrix": load_substitution_matrix("matrices/BLOSUM62.txt"), 
        "a": 3.6534,
        "b": -3.2351,
        "R2": 0.1255
    },
    3: {
        "matrix": load_substitution_matrix("matrices/BLOSUM45.txt"), 
        "a": 4.1893,
        "b": -2.5435,
        "R2": 0.5426
    },
    4: {
        "matrix": load_substitution_matrix("matrices/PAM160.txt"), 
        "a": 15.515,
        "b": -14.111,
        "R2": 0.5691
    },
    5: {
        "matrix": load_substitution_matrix("matrices/BLOSUM45.txt"), 
        "a": 5.776,
        "b": -2.9128,
        "R2": 0.7573
    },
}      
        
        

#print(mismatch_dict)
#print(exp_data)
#print(position_scoring[0]["matrix"])

def get_weighted_efficiency(matrix_value, a, b, R2):
    x = (matrix_value - b)/a  
    return x * R2


def get_score(mismatch_dict, exp_data, position_scoring):

    scored_motifs = [] 
    ###.items splits dictionary into keys and values (values being the objects of each list)
    #=> first for-loop iterates through keys
    # => second for loop iterates through objects of each list  
    for mismatch_count, motif_hits in mismatch_dict.items():

        #names all elements in motif_hits for later reference
        #most importantly opens list of mismatches to "mismatch_details"
        #=> mismatch_details is a list of tupels that contain the information about the mismatches for each motif
        for pos, window_sequence, mismatch_details in motif_hits:

            score = 0 
            score += window_size - mismatch_count
        
            # iterates through mismatches in each mismatch details list
            for mm_pos, original_aa, window_aa in mismatch_details:

                mismatch = (mm_pos, original_aa, window_aa)

                if mismatch in exp_data:
                    efficiency = exp_data[mismatch]
                    score += efficiency

                else: 
                    if mm_pos == 6:
                        score += 0

                    else:
                        best_matrix = position_scoring[mm_pos]["matrix"]
                        matrix_val = best_matrix[(original_aa, window_aa)]
                        weighted_efficiency = get_weighted_efficiency(matrix_val, position_scoring[mm_pos]["a"], position_scoring[mm_pos]["b"], position_scoring[mm_pos]["R2"])
                        score += weighted_efficiency

            #has to be at same indentation as mismatch loop
            #=> loop walks through multiple mismatches and thereby adds up the efficieny values before they are appended to the list
            scored_motifs.append((pos, window_sequence, score, mismatch_details))

    # key=lambda motif defines motif as the individual element of the tupel, motif[2] then tells the sort function which is elelment is relevant for sorting 
    scored_motifs.sort(key=lambda motif: motif[2], reverse=True)

    return(scored_motifs)


scored_motifs = get_score(mismatch_dict, exp_data, position_scoring)

top_k = 20 

for motif in scored_motifs[:top_k]:
    print(motif[0], motif[1], motif[2], motif[3])


    
                    








