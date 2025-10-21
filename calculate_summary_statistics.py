"""Module for calculating summary statistics from protein binding site predictions provided in process_boltz_from_index."""
import re
import pandas as pd
from pathlib import Path
from typing import Dict, Set, Tuple


def _norm_group_str(sites: list[int]) -> str:
    "Normalize/canonicalize string representing input peptide"
    return ";".join(str(x) for x in sorted(sites))


def parse_site_groups(index_file: Path, index_sheet: str, index_header: int = 2) -> Dict[Tuple[str, int], Set[str]]:
    """
    Build a site -> peptide mapping
    """
    ind_df = pd.read_excel(index_file, sheet_name=index_sheet, header=index_header)

    # This is a dictionary of {(accession, site): set(normalized strings)}
    site_to_peptides: Dict[Tuple[str, int], Set[str]] = {}

    for _, row in ind_df.iterrows():
        accession = str(row["Uniprot ID"]).strip()

        # robust split on ; , or whitespace
        sites_str = str(row["Site"])
        sites = [int(s) for s in re.split(r"[;,\s]+", sites_str.strip()) if s]

        if not sites:
            continue

        group_str = _norm_group_str(sites)

        for s in sites:
            key = (accession, int(s))
            site_to_peptides.setdefault(key, set()).add(group_str)

    return site_to_peptides


def calculate_bulk_summary_statistics(df: pd.DataFrame, index_file: Path, index_sheet: str, index_header: int = 2) -> pd.DataFrame:
    """
    Reduce the full input dataframe to a summary dataframe representing the predicted structure that places
    the ligand closest to each cysteine on each peptide.

    Peptides are denoted as groups of sites in the input index file, listed on a single line and separated
    by semicolons.
    """
    if df.empty:
        return df

    df = df.copy()
    df["distance_to_cys"] = pd.to_numeric(df["distance_to_cys"], errors="coerce")

    # Map (accession, site) to set of group strings
    site_to_groups = parse_site_groups(index_file, index_sheet, index_header)

    # For each prediction row, get ALL groups this [accession, cys_site] belongs to.
    # If none, fall back to the singleton site label.
    def groups_for_row(row):
        acc = str(row["accession"]).strip()
        try:
            site = int(row["cys_site"])
        except Exception:
            return {str(row["cys_site"])}
        groups = site_to_groups.get((acc, site))
        if not groups:
            return {str(site)}
        return groups

    df["site_group"] = df.apply(groups_for_row, axis=1)

    # expand to one row per site_group/peptide
    df = df.explode("site_group", ignore_index=True)

    # find the closest residue per [accession, site_group]
    groupby_cols = ["accession", "site_group"]
    idx = df.groupby(groupby_cols, dropna=False)["distance_to_cys"].idxmin()
    idx = idx.dropna()

    summary_df = df.loc[idx].reset_index(drop=True)

    # Add fraction of rows < 5 per group
    frac_near_cys = df.groupby(groupby_cols)["distance_to_cys"].apply(lambda x: (x < 5).mean())
    summary_df = summary_df.merge(frac_near_cys.rename("fraction_near_cys"), left_on=groupby_cols, right_index=True)

    return summary_df
