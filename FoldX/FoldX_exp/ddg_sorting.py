import json
import csv
import pandas as pd 
import math 
from pathlib import Path

F1 = "DGLSGAT"
position = 2

work_dir = Path("/Users/yannikmeindl/ApbE/FoldX/FoldX_exp")
csv_file = work_dir/"experimental_data.csv"

with open("ddg_per_motif.json", "r") as file:
    ddg_per_motif = json.load(file)

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
        motif = f"{F1[:mut_pos]}{row['Mutation']}{F1[(mut_pos+1):]}"

        all_motifs[motif] = row["Efficiency"]

    return all_motifs

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

            residue_scores[residue] = score

        motif_scores[motif] = residue_scores    
    
    return motif_scores

#--------------------Creates list of best ddg values per motifs of a given (input) position--------------------# 

def get_ddg_of_pos(ddg_per_motif, all_motifs):

    ddg_of_input_pos = []

    for motif, positions in ddg_per_motif.items():

        for residue, values in positions.items():

            if residue == str(position):

                efficiency = all_motifs[motif]
                ddg_of_input_pos.append((motif, values, efficiency))
    
    ddg_of_input_pos.sort(key=lambda x: abs(x[1]["ddG"]))  

    return ddg_of_input_pos

#--------------------Creates list of best ddg values per motifs of a given (input) position--------------------# 

def get_score_of_pos(motif_scores, all_motifs):

    scores_of_input_pos = []

    for motif, positions in motif_scores.items():

        for residue, score in positions.items():

            if residue == str(position):

                efficiency = all_motifs[motif]
                scores_of_input_pos.append((motif, score, efficiency))
                scores_of_input_pos.sort(key=lambda x: x[1])

    return scores_of_input_pos

#--------------------Creates list of best scores overall--------------------# 

def best_overall(motif_scores, all_motifs):

    best_all = []

    for motif, positions in motif_scores.items():

        for residue, score in positions.items():

            efficiency = all_motifs[motif]
            best_all.append((motif, score, efficiency, residue))
    
    best_all.sort(key=lambda x: x[1], reverse = True)
    
    return best_all 


exp_data = open_exp_data(csv_file)
all_motifs = get_motifs(exp_data)
motif_scores = get_energy_scores(ddg_per_motif)


best_ddg_per_motif = get_ddg_of_pos(ddg_per_motif, all_motifs)
best_score_per_motif = get_score_of_pos(motif_scores, all_motifs)
best = best_overall(motif_scores, all_motifs)

print(best)
print(best_ddg_per_motif)
print(best_score_per_motif)



        





