"""Module for classifying liganding events based on distance thresholds to orthosteric sites."""

import pandas as pd


def add_boolean_classifications(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add boolean classification columns to the dataframe based on distance thresholds.

    Args:
        df: Input dataframe with distance columns

    Returns:
        DataFrame with added boolean classification columns
    """
    # Check if distance_to_nearest_orthosteric_site column exists and has non-null values
    df["has_orthosteric_distance"] = (
        df["distance_to_nearest_orthosteric_site"].notna()
        if "distance_to_nearest_orthosteric_site" in df.columns
        else False
    )

    # Check if distance_to_cys <= 5
    df["ligand_near_cys"] = (
        pd.to_numeric(df["distance_to_cys"], errors="coerce") <= 5
        if "distance_to_cys" in df.columns
        else False
    )

    # Check if distance_to_nearest_orthosteric_site <= 5
    df["ligand_near_orthosteric"] = (
        pd.to_numeric(df["distance_to_nearest_orthosteric_site"], errors="coerce") <= 5
        if "distance_to_nearest_orthosteric_site" in df.columns
        else False
    )

    # Check if cys_to_nearest_orthosteric_distance <= 5
    df["cys_near_orthosteric"] = (
        pd.to_numeric(df["cys_to_nearest_orthosteric_distance"], errors="coerce") <= 5
        if "cys_to_nearest_orthosteric_distance" in df.columns
        else False
    )

    return df


def add_scenarios(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate named scenario annotations for input df with boolean distance annotations.
    """
    # Initialize boolean columns
    df['Ligand Bridging'] = False
    df['Direct Proximity'] = False
    df['Edge (10 Å)'] = False
    df['Intermediate'] = False
    df['Distant'] = False
    df['orthosteric event'] = False
    df['non-orthosteric event'] = False
    df['other'] = False

    # Only get rows with orthosteric data
    mask = df['has_orthosteric_distance']

    if mask.any():
        pl_c = pd.to_numeric(df.loc[mask, 'distance_to_cys'], errors='coerce')
        pl_os = pd.to_numeric(df.loc[mask, 'distance_to_nearest_orthosteric_site'], errors='coerce')
        os_c = pd.to_numeric(df.loc[mask, 'cys_to_nearest_orthosteric_distance'], errors='coerce')

        # Compute conditions for cysteine orthostery
        ligand_bridging = (pl_c < 5) & (pl_os < 5)
        direct_proximity = os_c < 5
        edge = (os_c < 10) & (pl_os < 10) & (pl_c >= 5)
        intermediate = (os_c >= 10) & (os_c < 12.5)
        distant = os_c >= 12.5

        ortho_cys = (ligand_bridging | direct_proximity) & ~edge
        nonortho_cys = distant & ~edge
        other = edge | intermediate | (~ortho_cys & ~nonortho_cys & ~intermediate)

        # Assign values to their own rows
        df.loc[mask, 'Ligand Bridging'] = ligand_bridging
        df.loc[mask, 'Direct Proximity'] = direct_proximity
        df.loc[mask, 'Edge (10 Å)'] = edge
        df.loc[mask, 'Intermediate'] = intermediate
        df.loc[mask, 'Distant'] = distant
        df.loc[mask, 'orthosteric event'] = ortho_cys
        df.loc[mask, 'non-orthosteric event'] = nonortho_cys
        df.loc[mask, 'other'] = other

    # Handle no-orthosteric-annotation cases
    df['no orthosteric annotation'] = ~mask

    # Create scenario column for primary classification
    df['scenario'] = 'Unclassified'
    df.loc[df['no orthosteric annotation'], 'scenario'] = 'no orthosteric annotation'
    df.loc[df['Ligand Bridging'], 'scenario'] = 'Ligand Bridging'
    df.loc[df['Direct Proximity'] & ~df['Ligand Bridging'], 'scenario'] = 'Direct Proximity'
    df.loc[df['Edge (10 Å)'], 'scenario'] = 'Edge (10 Å)'
    df.loc[df['Intermediate'], 'scenario'] = 'Intermediate'
    df.loc[df['Distant'] & ~df['Edge (10 Å)'] & ~df['Intermediate'], 'scenario'] = 'Distant'

    # Add cysteine orthostery classification
    df['cysteine_orthostery_classification'] = 'ambiguous'  # default for ambiguous cases
    df.loc[~mask, 'cysteine_orthostery_classification'] = 'no orthosteric annotation'
    df.loc[df['orthosteric event'], 'cysteine_orthostery_classification'] = 'orthosteric'
    df.loc[df['non-orthosteric event'], 'cysteine_orthostery_classification'] = 'non-orthosteric'

    # Add ligand orthostery classification
    df['ligand_orthostery_classification'] = None
    df.loc[~mask, 'ligand_orthostery_classification'] = 'no orthosteric annotation'

    # For rows with orthosteric data, classify based on ligand distance to orthosteric site
    if mask.any():
        ligand_near_orthosteric = pd.to_numeric(df.loc[mask, 'distance_to_nearest_orthosteric_site'], errors='coerce') < 5
        df.loc[mask & ligand_near_orthosteric, 'ligand_orthostery_classification'] = 'orthosteric'
        df.loc[mask & ~ligand_near_orthosteric, 'ligand_orthostery_classification'] = 'non-orthosteric'

    return df


def classify_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add classification columns to the dataframe based on distance thresholds. Classifications
    follow the logic refined by Bruno.

    This function orchestrates the classification process by:
    1. Adding boolean distance-based classifications
    2. Adding scenario classifications based on combinations of boolean flags

    Args:
        df: Input dataframe with distance columns

    Returns:
        DataFrame with all classification columns added
    """
    # Create a copy to avoid modifying the original
    df = df.copy()

    df = add_boolean_classifications(df)

    df = add_scenarios(df)

    # df.to_csv("test_out.csv", index=False)

    return df
