motifs = ["ABC", "CDE", "DFE", "GHJ", "LASD"]

for index, motif in enumerate(motifs):
    print(index,":", motif)


mismatches_per_motif = {
    0: [],
    1: [], 
    2: [], 
    3: [], 
    4: [],
    5: [],
}

print(type(mismatches_per_motif))