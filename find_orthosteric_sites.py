#!/usr/bin/env python3
"""Module for fetching and processing ligand binding information from UniProt."""

import time

import pandas as pd
import requests

# Ligands to exclude from analysis
# Note that corner case criteria for inclusion/exclusion are defined in _determine_ligand_status()
EXCLUDED_LIGANDS = [
    # This should catch most cofactors
    "Zn", "Zinc", "Zn(2+)", "Zn2+", "zinc",
    "Mg", "Magnesium", "Mg(2+)", "Mg2+", "magnesium",
    "Ca", "Calcium", "Ca(2+)", "Ca2+", "calcium",
    "Fe", "Iron", "Fe(2+)", "Fe2+", "Fe(3+)", "Fe3+", "iron",
    "Cu", "Copper", "Cu(2+)", "Cu2+", "Cu(1+)", "Cu+", "copper",
    "Mn", "Manganese", "Mn(2+)", "Mn2+", "manganese",
    "Co", "Cobalt", "Co(2+)", "Co2+", "cobalt",
    "Ni", "Nickel", "Ni(2+)", "Ni2+", "nickel",
    "K", "Potassium", "K(+)", "K+", "potassium",
    "Na", "Sodium", "Na(+)", "Na+", "sodium",
    "Mo", "Molybdenum", "molybdenum",
    "V", "Vanadium", "vanadium",
    "W", "Tungsten", "tungsten",
    "Al", "Aluminum", "Aluminium", "Al(3+)", "Al3+",
    "Cd", "Cadmium", "Cd(2+)", "Cd2+",
    "Hg", "Mercury", "Hg(2+)", "Hg2+",
    "Pb", "Lead", "Pb(2+)", "Pb2+",
]


def fetch_uniprot_data(accession: str, max_retries: int = 3) -> dict | None:
    """Fetch complete UniProt data for an accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 404:
                print(f"  UniProt entry not found for {accession}")
                break
            if attempt < max_retries - 1:
                time.sleep(1)
        except requests.RequestException as e:
            print(f"  Error fetching UniProt data for {accession}: {e}")
            if attempt < max_retries - 1:
                time.sleep(1)
    return None


def get_gene_name(uniprot_data: dict) -> str | None:
    """Extract primary gene name from UniProt data."""
    genes = uniprot_data.get("genes", [])
    if genes:
        gene_names = genes[0].get("geneName", {})
        if gene_names:
            return gene_names.get("value")
    return None


def parse_binding_sites(uniprot_data: dict, gene_name: str | None) -> list[dict]:
    """Parse binding site and active site information from UniProt JSON data."""
    binding_sites = []
    features = uniprot_data.get("features", [])
    for feature in features:
        feature_type = feature.get("type")
        if feature_type not in ("Binding site", "Active site", "Mutagenesis"):
            continue

        description = feature.get("description")

        ligand_name = None
        ligand_note = None

        if feature_type == "Binding site":
            ligand_info = feature.get("ligand", {})
            ligand_name = ligand_info.get("name")
            ligand_note = ligand_info.get("note")  # Extract note field
        elif feature_type == "Active site":
            ligand_name = feature.get("description")
            if len(ligand_name) < 2:
                ligand_name = "unnamed active site"

        location = feature.get("location", {})
        # Extract position(s)
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        # Format binding site
        if start and end and start != end:
            # For sequential ranges, use the center residue
            center = (int(start) + int(end)) // 2
            binding_site = str(center)
        elif start:
            binding_site = str(start)
        else:
            print(f"No position found for {feature}")
            continue

        if feature_type == "Mutagenesis":
            # Important to put this in name
            ligand_name = f"Res {binding_site} mutation"

        evidences = feature.get("evidences", [])
        eco_codes = [ev.get("evidenceCode") for ev in evidences if ev.get("evidenceCode")]
        if eco_codes:
            print(f"  ECO codes for {feature_type} at {binding_site}: {', '.join(sorted(set(eco_codes)))}")

        binding_sites.append({
            "gene_name": gene_name,
            "orthosteric_ligand_name": ligand_name,
            "orthosteric_ligand_binding_site": binding_site,
            "orthosteric_ligand_note": ligand_note,  # Add note to output
            "description": description,
            "feature_type": feature_type,
        })

    return binding_sites


def _determine_ligand_status(
    ligand_name: str,
    ligand_note: str | None,
    feature_type: str,
    exclusions: list[str],
    exclude_mutations: bool,
) -> str:
    """Determine if a ligand should be included or excluded.

    Metals are excluded unless they have 'catalytic' in the notes field.
    Mutations are excluded if exclude_mutations is True.
    """
    if exclude_mutations and feature_type == "Mutagenesis":
        return "excluded"

    if not exclusions:
        return "included"

    if not ligand_name:
        return "included"

    exclusions_lower = [e.lower() for e in exclusions]
    ligand_lower = ligand_name.lower()

    if ligand_lower in exclusions_lower:
        # Check if it has catalytic in notes (case-insensitive)
        if ligand_note and "catalytic" in ligand_note.lower():
            return "included"
        return "excluded"

    return "included"


def get_known_ligands(
    accession: str,
    exclude_metals: bool = True,
    exclude_mutations: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get known ligands for a protein from UniProt.

    Args:
        accession: UniProt accession ID
        exclude_metals: Whether to exclude metal ions from results
        exclude_mutations: Whether to exclude mutagenesis sites from results

    Returns:
        Tuple of (filtered_df, unfiltered_df with status column)
    """
    # Fetch UniProt data
    uniprot_data = fetch_uniprot_data(accession)
    if not uniprot_data:
        return pd.DataFrame(), pd.DataFrame()

    # Get gene name and parse binding sites
    gene_name = get_gene_name(uniprot_data)
    binding_sites = parse_binding_sites(uniprot_data, gene_name)

    if not binding_sites:
        return pd.DataFrame(), pd.DataFrame()

    # Convert to DataFrame
    df = pd.DataFrame(binding_sites)
    df["accession"] = accession

    # Build exclusion list
    exclusions = []
    if exclude_metals:
        exclusions.extend(EXCLUDED_LIGANDS)

    # Add status column
    df["status"] = df.apply(
        lambda row: _determine_ligand_status(
            row["orthosteric_ligand_name"],
            row["orthosteric_ligand_note"],
            row["feature_type"],
            exclusions,
            exclude_mutations,
        ),
        axis=1,
    )

    # Reorder columns for unfiltered df
    column_order = [
        "gene_name",
        "accession",
        "orthosteric_ligand_name",
        "orthosteric_ligand_binding_site",
        "orthosteric_ligand_note",
        "description",
        "feature_type",
        "status",
    ]
    df_unfiltered = df[column_order].copy()

    # Create filtered dataframe (exclude status and new columns)
    df_filtered = df[df["status"] == "included"][column_order[:5]].copy()

    return df_filtered, df_unfiltered


def get_ligands_batch(
    accessions: list[str],
    exclude_metals: bool = True,
    exclude_mutations: bool = False,
    delay: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch ligand information for multiple proteins.

    Args:
        accessions: List of UniProt accession numbers
        exclude_metals: Whether to exclude metal ions
        exclude_mutations: Whether to exclude mutagenesis sites
        delay: Delay between API calls in seconds

    Returns:
        Tuple of (filtered combined DataFrame, unfiltered combined DataFrame with status)
    """
    all_ligands_filtered = []
    all_ligands_unfiltered = []

    for i, accession in enumerate(accessions, 1):
        print(f"Fetching annotations for {accession} ({i}/{len(accessions)})")

        df_filtered, df_unfiltered = get_known_ligands(
            accession,
            exclude_metals=exclude_metals,
            exclude_mutations=exclude_mutations,
        )

        if not df_filtered.empty:
            all_ligands_filtered.append(df_filtered)
        if not df_unfiltered.empty:
            all_ligands_unfiltered.append(df_unfiltered)

        if i < len(accessions):
            time.sleep(delay)

    # Combine results
    if all_ligands_filtered:
        combined_filtered = pd.concat(all_ligands_filtered, ignore_index=True)
    else:
        combined_filtered = pd.DataFrame()

    if all_ligands_unfiltered:
        combined_unfiltered = pd.concat(all_ligands_unfiltered, ignore_index=True)
    else:
        combined_unfiltered = pd.DataFrame()

    return combined_filtered, combined_unfiltered