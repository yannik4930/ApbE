import pandas as pd

fasta_seq = "MSYYHHHHHHDYDIPTTENLYFQGAMGSLKERKLAKKRDELQRYVLMAADVNLGQGNEFRDIFAKSVKPLLINLDTGKVDSDANVLDFDERMAAINPETSSTPKKDIAKIKTRANDARVFKVFDDSGKLSSVVVPFYGKGLWSMIYGYVAVEPDFNTIKGVVVYEHGETPGIGDFVTDPHWLSLWKGKQLFDDKGKFAMRLVKGGVKEGDIHGVDAVSGATMTGRGVQRAMEFWFGVEGFQTFFNQLKASADQGELGGAK"
pdb_seq = "SLKERKLAKKRDELQRYVLMAADVNLGQGNEFRDIFAKSVKPLLINLDTGKVDSDANVLDFDERMAAINPETSSTPKKDIAKIKTRANDARVFKVFDDSGKLSSVVVPFYGKGLWSMIYGYVAVEPDFNTIKGVVVYEHGETPGIGDFVTDPHWLSLWKGKQLFDDKGKFAMRLVKGGVKEGDIHGVDAVSGATMTGRGVQRAMEFWFGVEGFQTFFNQLKAS"
alignment = ([[[ 27, 250]], [[  0, 223]]])


def alignment_table_on_fasta(fasta_seq, pdb_seq, alignment):
    # Grundtabelle: FASTA ist das Framework
    df = pd.DataFrame({
        "fasta_seq": list(fasta_seq),
        "pdb_seq": [""] * len(fasta_seq),
        "pdb_pos": pd.Series([pd.NA] * len(fasta_seq), dtype ="Int64"), 
    })

    df.index = range(1, len(fasta_seq) + 1)
    df.index.name = "fasta_pos"

    # *_blocks is a list of blocks with matching residues e.g fasta_blocks = [[5, 54], [78, 93]]
    fasta_blocks = alignment[0]
    pdb_blocks = alignment[1]

    # loop iterates through every matching block in the *_blocks lists
    for fasta_block, pdb_block in zip(fasta_blocks, pdb_blocks):
        fasta_start, fasta_end = fasta_block #defines the two numbers per block as the beginning and end of the matching sequence
        pdb_start, pdb_end = pdb_block

        #per every block, the numbers from fasta_start to fasta_end are matched with the numbers from pdb_start to pdb_end
        for fasta_index, pdb_index in zip(
            range(fasta_start, fasta_end),
            range(pdb_start, pdb_end)
        ):
            print(fasta_index, pdb_index)
            df.loc[fasta_index +1, "pdb_seq"] = pdb_seq[pdb_index]
            df.loc[fasta_index +1, "pdb_pos"] = pdb_index + 1

    return df

df = alignment_table_on_fasta(fasta_seq, pdb_seq, alignment)
print(df.to_string())


