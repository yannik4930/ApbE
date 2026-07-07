sequence = "DGLSGATMADRDSGSEQGGAALGDGLSGETDAGLSATDELSGATAGLSRATDRLWGATSGGSLGHPGQHETQEATLLLQGEEEGEED"
motif = "DGLSGAT"
window_size = len(motif)


def count_mismatches(sequence, motif, window_size, max_mismatches: str):


    windows = []

    #sliding window mechanism
    #range is defined as number of possible positions in the protein where a recognition motif could be
    #window is cut out of sequence from i to i + window size (i+window_size itself is not included)
    #window is attached to windows list 
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        windows.append(window) 

    #opens dictionary with lists for each number of mismatches
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
    #iterates through windows list while numbering the interation 
    for index, window in enumerate(windows):

        mismatches = 0
        pos_mismatch = []
    
    #assigns each character of the two strings to a/b, then checks identity
    #each pair of a's and b's get a number assigne (at which position the mismaatch is in the motif)
    # pos of the mismatch within the motif as well as the identity of the respective amino acids is appended to the mismatches_per_motif list
     
        for pos, (a, b) in enumerate(zip(motif, window)):
            if a != b:
                mismatches += 1
                pos_mismatch.append((pos, a, b,))

    #checks if sequence has less than 6 mismatches, then sorts them into the respective lists 
    
        if mismatches <= max_mismatches:
            mismatches_per_motif[mismatches].append((index, window, pos_mismatch))
    
    return mismatches_per_motif

# mismatch_dict = count_mismatches(sequence, motif, window_size, 7)

        

