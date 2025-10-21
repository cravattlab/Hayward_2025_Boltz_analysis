import pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def compare_stereochemistry(predicted: Chem.Mol | None | str, template: Chem.Mol) -> tuple[str, bool]:
    """
    Compare stereochemistry between two molecules.
    predicted: predicted structure (needs 3D stereochemistry assignment), None if kekulization failed, or "valence" if valence error
    template: template from SMILES (already has stereochemistry)

    Returns: (classification, is_match) where classification is "match", "enantiomer", "diastereomer", "failed to kekulize", or "failed valence"
    """
    if isinstance(predicted, str):
        return f"failed ({predicted})", False

    Chem.AssignStereochemistryFrom3D(predicted)

    predicted_chirality = dict(Chem.FindMolChiralCenters(predicted, includeUnassigned=True))
    template_chirality = dict(Chem.FindMolChiralCenters(template, includeUnassigned=True))

    if not predicted_chirality:  # No chiral centers
        return "match", True

    # Compare configurations (assuming same chiral centers)
    matches = sum(1 for idx in predicted_chirality if predicted_chirality[idx] == template_chirality.get(idx))
    total = len(predicted_chirality)

    match (matches, total):
        case (m, t) if m == t:
            return "match", True
        case (0, _):
            return "enantiomer", False
        case _:
            return "diastereomer", False


def get_ligand_structure(cif_path: Path, ligand_name: str = "LIG1") -> str:
    """
    Extract atoms from Boltz-2 prediction using the LIG1 name.

    Returns: str object representing the ligand in RDKit-compatible xyz format
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

    return "\n".join(lines) + "\n"


def get_ligand_RDKit_mol(ligand_xyz: str, smiles_template: str) -> tuple[Chem.Mol | None | str, Chem.Mol]:
    """
    Make an RDKit Mol object from xyz-formatted ligand data. Note that Hs don't get added here but this
    doesn't affect the stereochemistry evaluation.

    Returns both the generated mol and the template.
    """
    mol = Chem.MolFromXYZBlock(ligand_xyz)
    rdDetermineBonds.DetermineConnectivity(mol)

    template = Chem.MolFromSmiles(smiles_template)

    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        Chem.SanitizeMol(mol)
    except Chem.KekulizeException:
        print(f"Failed to kekulize molecule with SMILES: {smiles_template}")
        return "kekulization", template
    except Chem.AtomValenceException as e:
        print(f"Valence error for molecule with SMILES: {smiles_template} - {e}")
        return "valence", template

    return mol, template

def classify_df(
    df_main: pd.DataFrame,
    df_ligand: pd.DataFrame,
    enumerate_stereoisomers: bool = False,
) -> pd.DataFrame:
    """
    Evaluate predicted ligand stereochemistry and compare it to known stereochemistry.

    If `enumerate_stereoisomers` is False (default), compare only against the provided
    template SMILES. If True, enumerate all stereoisomers and compare against each.
    """
    df_m = df_main.copy()
    df_l = df_ligand.copy()

    all_rows = []

    for _, row in df_m.iterrows():
        model_filename = row["model_file_name"]
        model_path_str = row.get("model_path")

        cif_path = Path(model_path_str) if isinstance(model_path_str, str) else None
        if cif_path is not None and not cif_path.exists():
            cif_path = None

        if cif_path is None:
            folder_name = row.get("folder_name")
            if isinstance(folder_name, str):
                cif_path = next(Path(folder_name).rglob(model_filename), None)
            else:
                cif_path = None
        if not cif_path or not cif_path.exists():
            print(f"Warning: Could not find CIF file: {model_filename}")
            continue

        # Extract ligand structure
        ligand_xyz = get_ligand_structure(cif_path)

        # Get SMILES template
        ligand_name = row["ligand_name"]
        print(f"Evaluating {ligand_name} stereochemistry.")
        try:
            smiles_template = df_l.loc[df_l["Name"] == ligand_name, "Structure (SMILES)"].iloc[0]
        except IndexError:
            print(f"Warning: no SMILES found for ligand {ligand_name}")
            continue

        # Generate RDKit molecules
        mol, template = get_ligand_RDKit_mol(ligand_xyz, smiles_template)

        # Get stereochemistry info for mol if it didn't throw an RDKit error
        if not isinstance(mol, str):
            Chem.AssignStereochemistryFrom3D(mol)
            mol_chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            mol_smiles = Chem.MolToSmiles(mol)
        else:
            mol_chiral = None
            mol_smiles = None

        # Reference stereochemistry from the template
        template_chiral_reference = Chem.FindMolChiralCenters(template, includeUnassigned=True)

        # Choose what to compare against
        if enumerate_stereoisomers:
            opts = StereoEnumerationOptions(unique=True, onlyUnassigned=False)
            stereoisomers = list(EnumerateStereoisomers(template, options=opts))
        else:
            stereoisomers = [template]   # <- only compare to the template

        # Process each target (template or stereoisomers)
        for stereo_idx, target in enumerate(stereoisomers):
            target_chiral = Chem.FindMolChiralCenters(target, includeUnassigned=True)
            target_smiles = Chem.MolToSmiles(target)

            # Compare stereochemistry
            stereo_result, is_match = compare_stereochemistry(mol, target)

            # Create row
            new_row = row.to_dict()
            new_row['stereoisomer_index'] = stereo_idx
            if target_chiral == template_chiral_reference:
                new_row['stereoisomer_source'] = "template"
            else:
                new_row['stereoisomer_source'] = "synthetic"
            new_row['source_smiles'] = smiles_template
            new_row['template_SMILES'] = target_smiles
            new_row['predicted_ligand_SMILES'] = mol_smiles
            new_row['template_stereochemistry'] = str(target_chiral)
            if isinstance(mol, str) and mol == "valence":
                new_row['ligand_stereochemistry'] = "failed valence"
            elif isinstance(mol, str) and mol == "kekulization":
                new_row['ligand_stereochemistry'] = "failed to kekulize"
            else:
                new_row['ligand_stereochemistry'] = str(mol_chiral)
            new_row['stereochemistry_result'] = stereo_result
            new_row['is_stereochemistry_match'] = is_match

            all_rows.append(new_row)

    return pd.DataFrame(all_rows)