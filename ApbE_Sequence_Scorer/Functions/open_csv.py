import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
default_matrix = BASE_DIR/"experimental_data.csv"


def load_mutation_data(csv_file):
    mutation_lookup = {}

    with open(csv_file, newline="", encoding="utf-8-sig") as file:
        #opens csv, takes first row as a list of fieldnames, ";" tells it when each entry stops
        #then it takes each row as a list of strings and assigns them by number to the respective fieldname(first entry => first fieldname)
        reader = csv.DictReader(file, delimiter=";")

        for row in reader:
            #"csv.DictReader" takes every entry as a string, so some of the data types need to be changed
            position = int(row["Position"])
            original = row["Original"]
            mutation = row["Mutation"]
            efficiency = float(row["Efficiency"].replace(",", "."))

            #standard line for adding entrys to a dictionary, the part in [] is the key and "=" assigns it a value
            mutation_lookup[(position, original, mutation)] = efficiency

    return mutation_lookup

"""exp_data = load_mutation_data(default_matrix)

for i in exp_data:
    print(i)"""


