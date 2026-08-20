from pathlib import Path 
import pandas as pd

output_dir

def load_final_tables(final_tsv):

    final_df = pd.read_csv(
        final_tsv,
        sep="\t",
        index_col="Pos",
    )

    final_tables = {
        (motif, efficiency): motif_df
        .drop(columns=["Motif", "Efficiency"])
        .copy()

        for (motif, efficiency), motif_df
        in final_df.groupby(
            ["Motif", "Efficiency"],
            sort=False,
        )
    }

    return final_tables


def sort_output(final_tables, all_motifs, output_sorting, top_k):

    #best_per_motif
    if output_sorting == "best_per_motif":

        best_per_motif = {}

        for motif, residues_df in final_tables.items():

            efficiency = all_motifs[motif]

            valid_df = residues_df.dropna(
                subset=["O-Score"]
            )

            sorted_df = valid_df.sort_values(
                by = "O-Score",
                ascending = False,
            )

            if sorted_df.empty:
                continue

            best_per_motif[(motif, efficiency)] = sorted_df.head(top_k)

        best_per_motif = pd.read_csv(output_dir / "final.tsv", sep="\t")

        return(best_per_motif)


final_tables = load_final_tables(
    output_dir / "final.tsv"
)