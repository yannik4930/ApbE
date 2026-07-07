def load_substitution_matrix(matrix_file):
    matrix = {}

    with open(matrix_file, "r", encoding="utf-8") as file:
        lines = file.readlines()

    #lines takes whole row (first one here) of the matrix as a string, .split makes it a list, everything that is seperated by a space is then one entry of the list 
    header = lines[0].split()


    for line in lines[1:]:
        parts = line.split()

        row_aa = parts[0]
        scores = parts[1:]

        #header and scores get zipped together and thereby in the forloop internally redefined as col_aa and score (could also be named a and b)
        #the second line then then adds an entry to the matrix dictionary 
        for col_aa, score in zip(header, scores):
            matrix[(row_aa, col_aa)] = int(score)

    return matrix


#matrix = load_substitution_matrix("matrices/BLOSUM45.txt")
#print(matrix)