#!/usr/bin/env python3
import pandas as pd
import re
from functools import reduce

# -----------------------------
# Config: put your file paths here
# -----------------------------
FILES = [
    "GSE99248_raw_counts_GRCh38.p13_NCBI.tsv"
]

OUTPUT = "merged_all_counts.tsv"

# -----------------------------
# Helpers
# -----------------------------
def standardize_geneid_column(df: pd.DataFrame) -> pd.DataFrame:
    """Rename the gene id column (whatever its case/variant) to 'GeneID'."""
    id_candidates = [c for c in df.columns if c.lower() in ("geneid", "gene_id", "id_ref")]
    if not id_candidates:
        raise ValueError("Could not find a gene id column in one of the input files.")
    df = df.rename(columns={id_candidates[0]: "GeneID"})
    return df

def strip_xy_suffixes(col: str) -> str:
    """Remove trailing _x / _y (and optional .1, .2, …) that pandas adds when merging."""
    if col == "GeneID":
        return col
    return re.sub(r"(_x|_y)(\.\d+)?$", "", col)

# -----------------------------
# Main
# -----------------------------
def main():
    # Read all files as strings (safe for merging; we coerce to numeric later)
    dfs = [pd.read_csv(f, sep="\t", dtype=str) for f in FILES]

    # Ensure each has a 'GeneID' column
    dfs = [standardize_geneid_column(df) for df in dfs]

    # Merge all on GeneID (outer to keep any gene that appears in any file)
    merged = reduce(lambda l, r: pd.merge(l, r, on="GeneID", how="outer"), dfs)

    # Drop empty GeneIDs
    merged = merged[merged["GeneID"].notnull() & (merged["GeneID"] != "")]

    # Fill missing values with '0' (strings for now)
    merged = merged.fillna("0")

    # Clean column names (_x/_y)
    merged.columns = [strip_xy_suffixes(c) for c in merged.columns]

    # If stripping created duplicate column names, sum their values
    if merged.columns.duplicated().any():
        merged = merged.set_index("GeneID")
        merged = (
            merged.apply(pd.to_numeric, errors="coerce")  # turn counts numeric
                  .fillna(0)
                  .groupby(merged.columns, axis=1)        # group duplicated columns
                  .sum()
                  .astype(int)
                  .reset_index()
        )
    else:
        # Still cast to numeric ints for cleanliness
        numeric_cols = [c for c in merged.columns if c != "GeneID"]
        merged[numeric_cols] = (merged[numeric_cols]
                                .apply(pd.to_numeric, errors="coerce")
                                .fillna(0)
                                .astype(int))

    # Save
    merged.to_csv(OUTPUT, sep="\t", index=False)
    print(f"Merged file saved as {OUTPUT}")

if __name__ == "__main__":
    main()
