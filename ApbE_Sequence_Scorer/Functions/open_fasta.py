with open("fasta_test.txt", "r", encoding="utf-8") as f: #opens file and automatically closes it after end of indentation 
    fasta_content = f.readlines()                        #saves text as list, each entry represents one line


def read_fasta(fasta_content):
    #checks format of FASTA file
    if fasta_content[0].startswith(">"):
        pass
    else:
        print("No FASTA file found!\nMake sure file has the following format:\n>POI name\nPOI sequence")

    #opens empty string for sequence
    sequence = ""

    #checks if there is more than one FASTA in the text file 
    #stops skript after printing error message
    for line in fasta_content[1:]:

        if line.startswith(">"):
            print("Textfile contains more than one FASTA file. \nMake sure that textfile contains only one FASTAfile of the following format:\n>POI name\nPOI sequence")
            raise ValueError()  
        else:
            pass

        sequence += line.strip() #adds up entries of fasta_content list (without spaces, /n, etc.) in empty "sequence" string 

    return sequence

print(read_fasta(fasta_content))

