#!/usr/bin/env python3
"""Process Boltz2 prediction data and extract binding site information."""

import argparse
from pathlib import Path
import re
import shutil
import sys
import time
import warnings

from Bio.PDB.MMCIFParser import MMCIFParser
import numpy as np
import pandas as pd
import requests
from scipy.spatial.distance import cdist

import calculate_summary_statistics
import classify_liganding_events
import classify_ligand_stereochemistry
import find_orthosteric_sites
import validate_structures
import format_excel


def fetch_uniprot_gene_name(accession: str, max_retries: int = 3) -> dict[str, str]:
    """Fetch gene name and protein name from UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()

                # Extract gene name using nested get calls for safety
                gene_name = ""
                if genes := data.get("genes", []):
                    gene_name = genes[0].get("geneName", {}).get("value", "")

                # Extract protein name
                protein_name = ""
                if desc := data.get("proteinDescription", {}).get("recommendedName", {}):
                    protein_name = desc.get("fullName", {}).get("value", "")

                return {"gene_name": gene_name, "protein_name": protein_name}

            elif response.status_code == 404:
                print(f"  UniProt entry not found for {accession}")
                break

            if attempt < max_retries - 1:
                time.sleep(1)

        except requests.RequestException as e:
            print(f"  Error fetching UniProt data for {accession}: {e}")
            if attempt < max_retries - 1:
                time.sleep(1)

    return {"gene_name": "", "protein_name": ""}


def get_uniprot_info_batch(accessions: list[str]) -> dict[str, dict[str, str]]:
    """Fetch UniProt gene names for multiple accessions."""
    unique_accessions = list(set(accessions))
    print(f"Fetching UniProt data for {len(unique_accessions)} unique accessions...")

    uniprot_cache = {}
    for i, accession in enumerate(unique_accessions, 1):
        if i % 10 == 0:
            print(f"  Progress: {i}/{len(unique_accessions)}")

        uniprot_cache[accession] = fetch_uniprot_gene_name(accession)

        # Rate limiting
        if i < len(unique_accessions):
            time.sleep(0.1)

    return uniprot_cache


def calculate_min_distance(coords1: np.ndarray, coords2: np.ndarray) -> float | None:
    """Calculate minimum distance between two sets of coordinates."""
    if len(coords1) == 0 or len(coords2) == 0:
        return None
    return float(np.min(cdist(coords1, coords2)))


def get_binding_site(structure, cutoff: float = 5.0) -> list[int]:
    """Find residues within cutoff distance of ligand."""
    # Get ligand atoms
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " " and residue.resname != "HOH":
                    ligand_atoms.extend([atom.coord for atom in residue])

    if not ligand_atoms:
        return []

    ligand_coords = np.array(ligand_atoms)
    binding_residues = set()

    # Check each protein residue
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Protein residue
                    res_coords = np.array([atom.coord for atom in residue])
                    if calculate_min_distance(res_coords, ligand_coords) < cutoff:
                        binding_residues.add(residue.id[1])

    return sorted(binding_residues)


def get_residues_near_target(
    structure, target_resno: int, chain_id: str = "A",
    cutoff: float = 5.0, include_target: bool = False
) -> list[int]:
    """Find residues within cutoff distance of target residue."""
    target_atoms = []
    neighbor_residues = set()

    # Find target residue atoms
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == " ":  # Protein residue
                        if residue.id[1] == target_resno:
                            target_atoms = np.array([atom.coord for atom in residue])
                            if include_target:
                                neighbor_residues.add(target_resno)
                            break
                if len(target_atoms) > 0:  # Found target, no need to continue
                    break
        if len(target_atoms) > 0:  # Found target, no need to continue
            break

    if len(target_atoms) == 0:
        return []

    # Find neighbors
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == " " and residue.id[1] != target_resno:
                        res_coords = np.array([atom.coord for atom in residue])
                        if calculate_min_distance(res_coords, target_atoms) < cutoff:
                            neighbor_residues.add(residue.id[1])

    return sorted(neighbor_residues)


def calculate_plddt_stats(structure, residue_filter=None) -> dict:
    """Calculate pLDDT statistics for specified residues."""
    all_bfactors = []
    residue_scores = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue_filter is None or residue_filter(residue):
                    bfactors = [atom.bfactor for atom in residue]
                    all_bfactors.extend(bfactors)
                    if residue.id[0] == " ":  # Store per-residue for protein
                        residue_scores[residue.id[1]] = np.mean(bfactors)

    if not all_bfactors:
        return {"mean": np.nan, "stdev": np.nan, "residue_scores": residue_scores}

    return {
        "mean": float(np.mean(all_bfactors)),
        "stdev": float(np.std(all_bfactors)),
        "residue_scores": residue_scores
    }


def get_residue_ranges(structure) -> str:
    """Get residue ranges in the structure formatted as ranges."""
    residues = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Protein residues only
                    residues.add(residue.id[1])

    if not residues:
        return ""

    # Convert to sorted list and find ranges
    sorted_res = sorted(residues)
    ranges = []
    start = end = sorted_res[0]

    for res in sorted_res[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = res

    # Add final range
    ranges.append(f"{start}-{end}" if start != end else str(start))

    return ranges[0] if len(ranges) == 1 else f"({', '.join(ranges)})"


def parse_residue_ranges(ranges_str: str) -> set[int]:
    """Parse residue range string into a set of residue numbers."""
    if not ranges_str:
        return set()

    residues = set()
    ranges_str = ranges_str.strip("()")

    for part in ranges_str.split(","):
        part = part.strip()
        if "-" in part:
            try:
                start, end = map(int, part.split("-"))
                residues.update(range(start, end + 1))
            except ValueError:
                continue
        else:
            try:
                residues.add(int(part))
            except ValueError:
                continue

    return residues


def extract_domain_id(filename: str, accession: str, ligand_name: str) -> str:
    """Extract domain identifier from filename."""
    base_name = re.sub(r"_model_\d+\.cif$", "", filename)

    # Pattern: {accession}_{domain}_{ligand} or {accession}_{ligand}
    prefix = f"{accession}_"
    suffix = f"_{ligand_name}"

    if not base_name.startswith(prefix):
        return "unknown"

    if base_name == f"{accession}_{ligand_name}":
        return "main"

    if base_name.endswith(suffix):
        middle = base_name[len(prefix):-len(suffix)]
        return middle if middle else "main"

    return "main"


def process_structure_data(structure, cys_site: int, site_cutoff: float) -> dict:
    """Process all structure-related calculations."""
    result = {}

    # Get residue ranges
    result["residue_ranges"] = get_residue_ranges(structure)
    residue_set = parse_residue_ranges(result["residue_ranges"])
    result["cys_in_range"] = "Y" if cys_site in residue_set else "N"
    if cys_site not in residue_set:
        print(f"  Warning: Cys{cys_site} not in structure (residues: {result['residue_ranges']})")


    # Calculate pLDDT statistics
    global_stats = calculate_plddt_stats(structure, lambda r: r.id[0] == " ")
    result["mean_global_pLDDT"] = round(global_stats["mean"], 2) if not np.isnan(global_stats["mean"]) else ""
    result["stdev_global_pLDDT"] = round(global_stats["stdev"], 2) if not np.isnan(global_stats["stdev"]) else ""

    ligand_stats = calculate_plddt_stats(structure, lambda r: r.id[0] != " " and r.resname != "HOH")
    result["mean_ligand_pLDDT"] = round(ligand_stats["mean"], 2) if not np.isnan(ligand_stats["mean"]) else ""

    # Get binding site
    binding_site = get_binding_site(structure, site_cutoff)

    if not binding_site:
        result["calc_complete"] = "No ligand found"
        return result

    result["calc_complete"] = "Y"
    result["calc_resnos"] = ";".join(map(str, binding_site))
    result["match_with_cys_binding_site"] = "Y" if cys_site in binding_site else "N"

    # Get cysteine neighbors
    neighbors = get_residues_near_target(structure, cys_site, "A", site_cutoff, include_target=True)
    result["resnos_around_cys_site"] = ";".join(map(str, neighbors))

    # Calculate site-specific pLDDT
    residue_scores = global_stats["residue_scores"]

    binding_plddt = [residue_scores.get(r, np.nan) for r in binding_site]
    result["mean_binding_site_plddt"] = np.nanmean(binding_plddt)

    neighbor_plddt = [residue_scores.get(r, np.nan) for r in neighbors]
    result["mean_cys_sites_plddt"] = np.nanmean(neighbor_plddt)

    # Format detailed scores
    result["plddt_scores_for_each_residue_around_binding_site"] = ", ".join(
        f"{r}:{residue_scores.get(r, 0):.1f}" for r in binding_site
    )
    result["plddt_scores_for_each_residues_around_cys_site"] = ", ".join(
        f"{r}:{residue_scores.get(r, 0):.1f}" for r in neighbors
    )

    # Calculate ligand-cysteine distance
    cys_coords = []
    ligand_coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if chain.id == "A" and residue.id[0] == " " and residue.id[1] == cys_site:
                    cys_coords = np.array([atom.coord for atom in residue])
                elif residue.id[0] != " " and residue.resname != "HOH":
                    ligand_coords.extend([atom.coord for atom in residue])

    if len(cys_coords) > 0 and ligand_coords:
        distance = calculate_min_distance(cys_coords, np.array(ligand_coords))
        result["distance_to_cys"] = round(distance, 2) if distance else ""
    else:
        result["distance_to_cys"] = ""

    return result


def process_binding_sites(ind_df: pd.DataFrame, predictions_dir: Path, site_cutoff: float, skip_uniprot: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process binding site data from CIF files.

    Returns:
        Tuple of (models_data_df, failed_matches_df) where failed_matches_df contains
        protein-ligand pairs that had no matching CIF files.
    """

    if "Site" in ind_df.columns:
        ind_df["Site"] = ind_df["Site"].astype(str).str.split(r"\s*;\s*")
        ind_df = ind_df.explode("Site", ignore_index=True)
        ind_df["Site"] = ind_df["Site"].str.strip()
        ind_df = ind_df[ind_df["Site"] != ""]
        ind_df["Site"] = ind_df["Site"].astype(int)

    # Fetch UniProt data
    uniprot_cache = {}
    if not skip_uniprot:
        accessions = ind_df["Uniprot ID"].unique().tolist()
        uniprot_cache = get_uniprot_info_batch(accessions)

    # Find all CIF files
    cif_files = {f.name: f for f in predictions_dir.rglob("*_model_*.cif")}

    models_data = []
    failed_matches = []
    parser = MMCIFParser(QUIET=True)

    print(f"Processing {len(ind_df)} ligand-protein pairs from index file...")

    for _, row in ind_df.iterrows():
        accession = row["Uniprot ID"]
        ligand_name = row["Stereoprobe"]
        cys_site = int(row["Site"])
        gene_name = row["Protein"]

        # Get UniProt data
        uniprot_data = uniprot_cache.get(accession, {"gene_name": "", "protein_name": ""})

        # Find matching CIF files
        prefix = f"{accession}_"
        suffix_pattern = f"_{ligand_name}_model_"

        matching_files = [
            (name, path) for name, path in cif_files.items()
            if name.startswith(prefix) and suffix_pattern in name
        ]

        if not matching_files:
            print(f"  No CIF files found for {gene_name} with ligand {ligand_name}")
            # Track missing pair
            missing_pair = {
                # "entry": row.get("entry", ""),
                "accession": accession,
                "gene_name": uniprot_data["gene_name"], # use UniProt name
                "protein_name": uniprot_data["protein_name"],
                "ligand_name": ligand_name,
                "cys_site": cys_site,
                "reason": "No matching CIF files found"
            }
            # Add any other columns from the original index that might be useful
            for col in ind_df.columns:
                if col not in missing_pair:
                    missing_pair[col] = row.get(col, "")
            failed_matches.append(missing_pair)
            continue

        for filename, cif_path in matching_files:
            # Extract model rank
            if model_match := re.search(r"model_(\d+)", filename):
                model_rank = int(model_match.group(1))
            else:
                continue

            # print(f"  Processing: {filename}")

            base_name = re.sub(r"_model_\d+\.cif$", "", filename)
            domain_id = extract_domain_id(filename, accession, ligand_name)

            model_info = {
                # "entry": row["entry"],
                "accession": accession,
                "gene_name": uniprot_data["gene_name"],
                "ligand_name": ligand_name,
                "cys_site": cys_site,
                "domain_identifier": domain_id,
                "complex_name": base_name,
                "model_rank": model_rank,
                "model_file_name": filename,
                "folder_name": str(cif_path.parent.resolve()),
                "model_path": str(cif_path.resolve())
            }

            try:
                # Parse structure directly
                structure = parser.get_structure("protein", str(cif_path))

                # Process all structure data
                structure_data = process_structure_data(structure, cys_site, site_cutoff)
                model_info.update(structure_data)

            except Exception as e:
                print(f"  Error processing {filename}: {e}")
                model_info.update({
                    "calc_complete": f"Error: {str(e)}",
                    "residue_ranges": "",
                    "cys_in_range": "",
                    "distance_to_cys": "",
                    "mean_global_pLDDT": "",
                    "stdev_global_pLDDT": "",
                    "mean_ligand_pLDDT": ""
                })

            models_data.append(model_info)

    return pd.DataFrame(models_data), pd.DataFrame(failed_matches)

def add_all_orthosteric_distances(
    merged_df: pd.DataFrame,
    orthosteric_df: pd.DataFrame,
    predictions_dir: Path
) -> pd.DataFrame:
    """
    Add both ligand-to-orthosteric and cysteine-to-orthosteric distance columns.

    Returns a DataFrame with four new columns:
        - distance_to_nearest_orthosteric_site: Min distance from ligand to orthosteric site
        - nearest_orthosteric_ligand: Name of nearest orthosteric ligand from ligand
        - cys_to_nearest_orthosteric_distance: Min distance from cys to orthosteric site
        - nearest_orthosteric_from_cys: Name of nearest orthosteric ligand from cys
    """
    # Initialize all new columns
    merged_df["distance_to_nearest_orthosteric_site"] = pd.NA
    merged_df["nearest_orthosteric_ligand"] = ""
    merged_df["cys_to_nearest_orthosteric_distance"] = pd.NA
    merged_df["nearest_orthosteric_from_cys"] = ""

    # Group orthosteric sites by accession
    orthosteric_grouped = orthosteric_df.groupby("accession").apply(
        lambda x: list(zip(x["orthosteric_ligand_binding_site"],
                           x["orthosteric_ligand_name"]))
    ).to_dict()

    parser = MMCIFParser(QUIET=True)

    # Track processed files to avoid redundant parsing
    processed_count = 0

    for idx, row in merged_df.iterrows():
        # Skip if missing required data
        if pd.isna(row.get("accession")) or pd.isna(row.get("model_file_name")):
            continue

        accession = row["accession"]
        if accession not in orthosteric_grouped:
            continue

        cys_site = int(row["cys_site"]) if not pd.isna(row.get("cys_site")) else None

        # Find CIF file
        model_filename = row["model_file_name"]
        cif_path = next(predictions_dir.rglob(model_filename), None)

        if not cif_path or not cif_path.exists():
            print(f'Warning: Could not find CIF file: {model_filename}')
            continue

        try:
            # Parse structure directly
            structure = parser.get_structure("protein", str(cif_path))
            processed_count += 1

            # Extract ligand coordinates
            ligand_coords = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != " " and residue.resname != "HOH":
                            ligand_coords.extend([atom.coord for atom in residue])

            ligand_coords = np.array(ligand_coords) if ligand_coords else None

            # Extract cysteine coordinates if cys_site is provided
            cys_coords = None
            if cys_site:
                for model in structure:
                    for chain in model:
                        if chain.id == "A":  # Assuming chain A
                            for residue in chain:
                                if residue.id[0] == " " and residue.id[1] == cys_site:
                                    cys_coords = np.array([atom.coord for atom in residue])
                                    break
                            if cys_coords is not None:
                                break
                    if cys_coords is not None:
                        break

            # Calculate distances to orthosteric sites
            ligand_min_distance = float("inf")
            ligand_nearest = ""
            cys_min_distance = float("inf")
            cys_nearest = ""

            for ortho_resno, ligand_name in orthosteric_grouped[accession]:
                # Get orthosteric site residue coordinates
                ortho_resno_num = int(ortho_resno)
                ortho_coords = []

                for model in structure:
                    for chain in model:
                        if chain.id == "A":  # Assuming chain A
                            for residue in chain:
                                if residue.id[0] == " " and residue.id[1] == ortho_resno_num:
                                    ortho_coords = [atom.coord for atom in residue]
                                    break
                            if ortho_coords:
                                break
                    if ortho_coords:
                        break

                if not ortho_coords:
                    continue

                ortho_coords = np.array(ortho_coords)

                # Calculate ligand-to-orthosteric distance
                if ligand_coords is not None and len(ligand_coords) > 0:
                    dist = calculate_min_distance(ligand_coords, ortho_coords)
                    if dist is not None and dist < ligand_min_distance:
                        ligand_min_distance = dist
                        ligand_nearest = ligand_name

                # Calculate cys-to-orthosteric distance
                if cys_coords is not None and len(cys_coords) > 0:
                    dist = calculate_min_distance(cys_coords, ortho_coords)
                    if dist is not None and dist < cys_min_distance:
                        cys_min_distance = dist
                        cys_nearest = ligand_name

            # Update all distance columns
            if ligand_min_distance != float("inf"):
                merged_df.at[idx, "distance_to_nearest_orthosteric_site"] = round(ligand_min_distance, 2)
                merged_df.at[idx, "nearest_orthosteric_ligand"] = ligand_nearest

            if cys_min_distance != float("inf"):
                merged_df.at[idx, "cys_to_nearest_orthosteric_distance"] = round(cys_min_distance, 2)
                merged_df.at[idx, "nearest_orthosteric_from_cys"] = cys_nearest

        except Exception as e:
            print(f"Warning: Error processing {model_filename}: {e}")
            continue

    return merged_df


def main():
    parser = argparse.ArgumentParser(description="Process Boltz-2 prediction data")
    parser.add_argument("--predictions-dir", required=True, help="Path to directory of cif files. This directory must contain cif files with {accession}_{arbitrary_field}_{ligand_name}_model_{num}.cif")
    parser.add_argument("--output", default=f"western_analogs_output_{int(time.time())}.xlsx", help="Output Excel file path")
    parser.add_argument("--index-file", required=True, help="Path to index file with ligand information")
    parser.add_argument("--index-sheet", required=True, help="sites sheet name in index file")
    parser.add_argument("--ligand-sheet", help="compounds sheet name in index file")
    parser.add_argument("--site-cutoff", type=float, default=5.0, help="Distance cutoff for binding site")
    parser.add_argument("--filter-cys-domains", action="store_true",
                        help="Only include domains that contain the cysteine site in output")
    parser.add_argument("--include-extended-fields", action="store_true",
                        help="Include all verbose fields in output")
    parser.add_argument("--skip-uniprot", action="store_true",
                        help="Skip fetching data from UniProt API")
    parser.add_argument("--skip-orthosteric", action="store_true",
                        help="Skip fetching orthosteric site data from UniProt")
    parser.add_argument("--skip-validation", action="store_true",
                        help="Skip Phenix/PoseBusters validation")
    parser.add_argument("--include-metals", action="store_true",
                        help="Include metal ions in orthosteric sites analysis")
    parser.add_argument("--include-mutations", action="store_true",
                        help="Include mutagenesis sites in orthosteric sites analysis")
    parser.add_argument("--index-header", type=int, default=2, help="Header in the index sheet")
    parser.add_argument("--ligands-header", type=int, default=2, help="Header in the ligand sheet")
    parser.add_argument("--save-raw-data", action="store_true",
                        help="Skip Excel formatting and write full output to an excel file")

    args = parser.parse_args()

    # make sure phenix is available in the system path
    if not args.skip_validation:
        print("Checking for phenix...")
        if shutil.which("phenix") is None:
            sys.exit("Error: 'phenix' is not in your system path. Please source phenix_env.sh first.")
        else:
            print("Phenix found.")
    else:
        print("Skipping Phenix check (validation disabled).")

    print("\nStarting to process data...")

    if not args.ligand_sheet:
        warnings.warn("No ligand sheet provided, skipping stereochemical evaluation."
                      )
    predictions_dir = Path(args.predictions_dir)

    # Parse output path
    output_file = Path(args.output)
    output_dir = output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    main_df = pd.DataFrame()
    index_df = pd.DataFrame()
    failed_matches_df = pd.DataFrame()
    ligand_df = pd.DataFrame()

    # Process binding sites if index file provided
    if args.index_file:
        index_path = Path(args.index_file)
        if index_path.exists():
            index_df = pd.read_excel(index_path, sheet_name=args.index_sheet, header=args.index_header)
            print("Processing binding site information...")
            site_df, failed_matches_df = process_binding_sites(
                index_df, predictions_dir, args.site_cutoff, args.skip_uniprot
            )

            if not site_df.empty:
                main_df = pd.merge(main_df, site_df, on=["complex_name", "model_rank"], how="outer") \
                    if not main_df.empty else site_df
        else:
            print("No index file provided!")


    # Apply filters
    # Remove domains that do not contain target cysteine
    if args.filter_cys_domains and "cys_in_range" in main_df.columns:
        original_count = len(main_df)
        main_df = main_df[main_df["cys_in_range"] == "Y"]
        print(f"Filtered out {original_count - len(main_df)} domains without cysteine site")

    # Sort results on name & rank
    if not main_df.empty:
        main_df = main_df.sort_values(["complex_name", "model_rank"])

    # Fetch orthosteric site info
    orthosteric_df = pd.DataFrame()
    orthosteric_df_unfiltered = pd.DataFrame()
    if not args.skip_orthosteric and "accession" in main_df.columns:
        unique_accessions = main_df["accession"].dropna().unique().tolist()
        if unique_accessions:
            print(f"\nFetching orthosteric site data for {len(unique_accessions)} unique proteins...")
            orthosteric_df, orthosteric_df_unfiltered = find_orthosteric_sites.get_ligands_batch(
                unique_accessions,
                exclude_metals=not args.include_metals,  # Exclude metals by default
                exclude_mutations=not args.include_mutations,  # Exclude mutations by default
                delay=0.1  # rate limit queries
            )

            if not orthosteric_df.empty:
                print(f"Found {len(orthosteric_df)} orthosteric binding sites")
                print(f"Excluded {len(orthosteric_df_unfiltered) - len(orthosteric_df)} orthosteric ligands, likely all metal cofactors.")
                proteins_with_sites = len(orthosteric_df.groupby("accession").size())
                print(f"Proteins with sites: {proteins_with_sites}/{len(unique_accessions)}")

    # Print distance summary
    if "distance_to_cys" in main_df.columns:
        valid_distances = pd.to_numeric(main_df["distance_to_cys"], errors="coerce").dropna()
        if not valid_distances.empty:
            print(f"\nLigand-cysteine distances: min={valid_distances.min():.2f}Å, "
                  f"max={valid_distances.max():.2f}Å, mean={valid_distances.mean():.2f}Å")

    # Calculate orthosteric site vs ligand distances
    if not orthosteric_df.empty and not main_df.empty:
        print("\nCalculating distances to orthosteric sites...")
        main_df = add_all_orthosteric_distances(
            main_df, orthosteric_df, predictions_dir
        )

    # Add stereochemistry evaluation
    if args.ligand_sheet:
        ligand_df = pd.read_excel(index_path, sheet_name=args.ligand_sheet, header=args.ligands_header)
        main_df = classify_ligand_stereochemistry.classify_df(main_df, ligand_df, enumerate_stereoisomers=False)

    # Add liganding event categories
    main_df = classify_liganding_events.classify_df(main_df)
    if not args.skip_validation:
        pb_df = validate_structures.validate(main_df, ligand_df)
    else:
        print("Skipping structure validation step.")
        pb_df = pd.DataFrame()
    summary_df = calculate_summary_statistics.calculate_bulk_summary_statistics(main_df, args.index_file, args.index_sheet, args.index_header)

    dataframes = {
        'analysis_data': main_df,
        'per_peptide_summary': summary_df,
        'orthosteric_sites': orthosteric_df_unfiltered,
        'failed_matches': failed_matches_df,
        'models_validation': pb_df,
        'input_index': index_df
    }

    # Save with formatting
    format_excel.save_formatted_excel(output_file, dataframes, save_raw_data=args.save_raw_data)


if __name__ == "__main__":
    main()
