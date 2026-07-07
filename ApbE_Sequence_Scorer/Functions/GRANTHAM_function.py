import pandas as pd

ideal = "DGLSGAT"
motif = "WSTYLPV"
grantham = pd.read_csv("matrix.csv", index_col = "AA")

def motif_distance(grantham, ideal, motif):

    motif_score = 0

    for ideal_aa, motif_aa in zip(ideal, motif):
        motif_score += grantham.loc[ideal_aa, motif_aa]

    return motif_score

print(motif_distance(grantham, ideal, motif))


