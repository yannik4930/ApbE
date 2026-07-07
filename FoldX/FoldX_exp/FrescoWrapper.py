import subprocess
from pathlib import Path


#general
work_dir = Path("/Users/yannikmeindl/FoldX/FoldX_exp")
F1 = "DGLSGAT"

#get_fasta_seq
fasta = work_dir / "rcsb_pdb_4XHF.fasta"

#get_pdb_seq
pdb = work_dir/"4XHF.pdb"
repaired_cleaned_pdb = work_dir/"frescoTest"/"4XHF_cleaned_prepared.pdb"

#exp_data
csv_file = work_dir/"experimental_data.csv"

