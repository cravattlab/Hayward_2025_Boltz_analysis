"""Use PoseBusters to evaluate the physical validity and stereochemistry of each Boltz-2 prediction."""

from contextlib import contextmanager
from pathlib import Path
import pandas as pd
import re
import shlex
import subprocess
import tempfile

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO
from posebusters import PoseBusters
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import Descriptors
from rdkit.Chem import FindMolChiralCenters


def get_ligand_mol(cif_path: Path, template_mol: Chem.Mol, ligand_name: str = "LIG1") -> Chem.Mol:
    """
    Extract ligand atoms from Boltz-2 prediction.

    Returns: RDKit mol object representing the ligand
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", str(cif_path))

    atoms: list[tuple[str, float, float, float]] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == ligand_name:
                    for atom in residue:
                        elem = (atom.element or "").strip() or atom.get_name()[0]
                        x, y, z = map(float, atom.coord.tolist())
                        atoms.append((elem, x, y, z))

    lines = [str(len(atoms)), "LIG1"]
    for e, x, y, z in atoms:
        lines.append(f"{e} {x:.3f} {y:.3f} {z:.3f}")

    # XYZ-formatted intermediate
    ligand_xyz = "\n".join(lines) + "\n"

    ligand_mol = Chem.MolFromXYZBlock(ligand_xyz)
    if ligand_mol is None:
        return None

    rdDetermineBonds.DetermineConnectivity(ligand_mol)

    try:
        ligand_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, ligand_mol)
        Chem.SanitizeMol(ligand_mol)
    except Chem.KekulizeException:
        print(f"Failed to kekulize molecule in: {cif_path}")
        return "kekulization", template_mol
    except Chem.AtomValenceException as e:
        print(f"Valence error for molecule in: {cif_path} - {e}")
        return "valence", template_mol

    Chem.AssignStereochemistryFrom3D(ligand_mol)

    return ligand_mol


@contextmanager
def get_protein_pdb(cif_path: Path, ligand_name: str = "LIG1") -> None:
    """
    Extract atoms from Boltz-2 prediction representing the protein
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", str(cif_path))

    # Remove ligand atoms
    for model in structure:
        for chain in model:
            atoms_to_remove = []
            for residue in chain:
                if residue.resname == ligand_name:
                    atoms_to_remove.append(residue.id)
            for res_id in atoms_to_remove:
                chain.detach_child(res_id)

    # Save temporary PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp.name)
        tmp_path = Path(tmp.name)

    print(f"- Saved temp PDB file to {str(tmp_path)}")

    try:
        yield tmp_path
    finally:
        tmp_path.unlink()


def phenix_validate_protein(pdb_path: str, timeout: int = 600) -> dict:
    """Run phenix.model_statistics and return geometry + MolProbity stats."""

    # Run phenix
    p = subprocess.run(
        shlex.split(f"phenix.model_statistics {shlex.quote(pdb_path)}"),
        capture_output=True, text=True, timeout=timeout
    )

    result = {
        "phenix_available": True,
        "phenix_error": None,
        "bond_mean": None,
        "bond_max": None,
        "bond_rmsz": None,
        "angle_mean": None,
        "angle_max": None,
        "angle_rmsz": None,
        "chirality_mean": None,
        "planarity_mean": None,
        "dihedral_mean": None,
        "min_nonbonded_distance": None,
        "clashscore": None,
        "rama_outliers_pct": None,
        "rama_allowed_pct": None,
        "rama_favored_pct": None,
        "rotamer_outliers_pct": None,
        "rotamer_allowed_pct": None,
        "rotamer_favored_pct": None,
        "cbeta_deviations_pct": None,
    }

    if p.returncode != 0:
        result["phenix_error"] = (p.stderr or p.stdout or "phenix.model_statistics failed").strip()
        return result

    text = f"{p.stdout}\n{p.stderr}"
    flags = re.I | re.M | re.S

    # Geometry restraints
    if m := re.search(r"^\s*Bond\s*:\s*([0-9.]+)\s+([0-9.]+)\s+\d+\s+Z=\s*([0-9.]+)", text, flags):
        result["bond_mean"] = float(m.group(1))
        result["bond_max"] = float(m.group(2))
        result["bond_rmsz"] = float(m.group(3))
    if m := re.search(r"^\s*Angle\s*:\s*([0-9.]+)\s+([0-9.]+)\s+\d+\s+Z=\s*([0-9.]+)", text, flags):
        result["angle_mean"] = float(m.group(1))
        result["angle_max"] = float(m.group(2))
        result["angle_rmsz"] = float(m.group(3))
    if m := re.search(r"^\s*Chirality\s*:\s*([0-9.]+)", text, flags):
        result["chirality_mean"] = float(m.group(1))
    if m := re.search(r"^\s*Planarity\s*:\s*([0-9.]+)", text, flags):
        result["planarity_mean"] = float(m.group(1))
    if m := re.search(r"^\s*Dihedral\s*:\s*([0-9.]+)", text, flags):
        result["dihedral_mean"] = float(m.group(1))
    if m := re.search(r"^\s*Min\s+Nonbonded\s+Distance\s*:\s*([0-9.]+)", text, flags):
        result["min_nonbonded_distance"] = float(m.group(1))

    # MolProbity stats
    if m := re.search(r"All-atom\s+Clashscore\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)", text, flags):
        result["clashscore"] = float(m.group(1))
    if m := re.search(r"Ramachandran.*?Outliers\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rama_outliers_pct"] = float(m.group(1))
    if m := re.search(r"Ramachandran.*?Allowed\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rama_allowed_pct"] = float(m.group(1))
    if m := re.search(r"Ramachandran.*?Favou?red\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rama_favored_pct"] = float(m.group(1))
    if m := re.search(r"\bRotamer\b.*?Outliers\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rotamer_outliers_pct"] = float(m.group(1))
    if m := re.search(r"\bRotamer\b.*?Allowed\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rotamer_allowed_pct"] = float(m.group(1))
    if m := re.search(r"\bRotamer\b.*?Favou?red\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["rotamer_favored_pct"] = float(m.group(1))
    if m := re.search(r"C(?:-?beta|β)\s+Deviations\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags):
        result["cbeta_deviations_pct"] = float(m.group(1))

    return result

def run_validation(df_main: pd.DataFrame, df_ligand: pd.DataFrame) -> pd.DataFrame:
    """
    Run PoseBusters analysis on protein-ligand complexes from main dataframe.

    Returns:
        DataFrame with PoseBusters analysis results for each model (one row per model)
    """
    if df_main.empty or df_ligand.empty:
        return pd.DataFrame()

    results = []

    for index, row in df_main.iterrows():

        print(f"\nRunning PoseBusters + phenix analysis on {index+1}/{len(df_main)} models...\n")

        result = {
            "model_rank": row["model_rank"],
            "accession": row["accession"],
            "gene_name": row["gene_name"],
            "ligand_name": row["ligand_name"],
            "cys_site": row["cys_site"],
        }

        model_filename = row["model_file_name"]
        model_path_str = row["model_path"]

        cif_path = Path(model_path_str) if isinstance(model_path_str, str) else None
        if cif_path is not None and not cif_path.exists():
            cif_path = None

        if cif_path is None:
            folder_name = row["folder_name"]
            if isinstance(folder_name, str):
                cif_path = next(Path(folder_name).rglob(model_filename), None)
            else:
                cif_path = None

        if not cif_path or not cif_path.exists():
            print(f'Warning: Could not find CIF file: {model_filename}')
            results.append(pd.DataFrame([{**result, "pb_error": f"Missing CIF: {model_filename}", "pb_valid": False}]))
            continue

        print(f"Analyzing {cif_path}")

        try:
            # get template smiles
            template_smiles = row["template_SMILES"]
            if not template_smiles or not isinstance(template_smiles, str):
                results.append(pd.DataFrame([{**result, "pb_error": "error getting template_SMILES", "pb_valid": False}]))
                continue

            # get and condition template mol
            template_mol = Chem.MolFromSmiles(template_smiles, sanitize=True)
            params = AllChem.ETKDGv3()
            params.pruneRmsThresh = 0.5
            _ = AllChem.EmbedMolecule(template_mol, params)
            AllChem.UFFOptimizeMolecule(template_mol, maxIters=200)
            if template_mol is None or template_mol.GetNumAtoms() == 0:
                results.append(pd.DataFrame([{**result, "pb_error": "invalid or empty template mol", "pb_valid": False}]))
                continue

            # get ligand mol
            ligand_mol = get_ligand_mol(cif_path, template_mol)
            if ligand_mol is None or ligand_mol.GetNumAtoms() == 0:
                results.append(pd.DataFrame([{**result, "pb_error": "error getting ligand mol", "pb_valid": False}]))
                continue

            print(f"- Template MW: {Descriptors.HeavyAtomMolWt(template_mol)}")
            print("- Chiral centers:", FindMolChiralCenters(template_mol, includeUnassigned=True))
            print(f"- Ligand MW: {Descriptors.HeavyAtomMolWt(ligand_mol)}")
            print("- Chiral centers:", FindMolChiralCenters(ligand_mol, includeUnassigned=True))

            with get_protein_pdb(cif_path) as protein_pdb_path:
                buster = PoseBusters(config="redock_fast").bust(
                # buster = PoseBusters(config="redock").bust( # this option is much slower
                    [ligand_mol],
                    template_mol,
                    str(protein_pdb_path),
                    full_report=True
                )

                # Run phenix analysis using the temp PDB file
                phenix_summary = phenix_validate_protein(str(protein_pdb_path))

            phenix_df = pd.DataFrame([phenix_summary])
            result_df = pd.DataFrame([result])

            merged = pd.concat([result_df.reset_index(drop=True),
                                buster.reset_index(drop=True),
                                phenix_df.reset_index(drop=True)], axis=1)

            results.append(merged)

        except Exception as e:
            print(f"Error: {e}")
            results.append(pd.DataFrame([{**result, "pb_error": str(e), "pb_valid": False}]))

    if not results:
        return pd.DataFrame()

    return pd.concat(results, ignore_index=True, sort=False)


def validate(
    df_main: pd.DataFrame,
    df_ligand: pd.DataFrame,
    skip_posebusters: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate PoseBusters and phenix results for each model in the input dataframe.
    """
    if skip_posebusters or df_main.empty:
        return df_main, pd.DataFrame(), pd.DataFrame()

    print("\nRunning PoseBusters & phenix validation...")

    # Run PoseBusters and phenix on all models
    pb_df = run_validation(df_main, df_ligand)

    return pb_df
