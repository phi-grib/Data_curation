"""
    Code for extracting ChEMBL datasets using
    a target.

    Target can be a string (cox2) or a ChEMBL ID (CHEMBL230).

    TODO:Might incoporate other arguments like species (e.g: Homo sapiens)

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/04/2021, 11:13 PM
"""

import pandas as pd

from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import PandasTools

def get_dataframe_from_target(chembl_id: str) -> pd.DataFrame:
    """
        Creates a client to ChEMBL and uses chembl_id to extract
        all the compounds related to that target.
        Returns a dataframe with all the information

        :param target_name: string with a valid chembl_id

        :return df_res: dataframe containing the results
    """

    # Create an activity query
    activities = new_client.activity

    # Make chembl_id uppercase
    chembl_id_upper = chembl_id.upper()
    
    # Select only activities with a pchembl_value (-log(IC50, Ki, Kd, EC50...).
    # We also use the chembl flags to remove the duplicates and the records where there is a validity comment
    response = activities.filter(target_chembl_id=chembl_id_upper, pchembl_value__isnull=False,\
                                potential_duplicate=False, data_validity_comment__isnull=True )

    # create a dataframe with the activity data
    df_activities = pd.DataFrame(response)
    assays = new_client.assay

    # select assays.
    response = assays.filter(assay_chembl_id__in=list(df_activities.assay_chembl_id.unique()))

    # create a dataframe with the assay data
    df_assays = pd.DataFrame(response)

    # keep only the assays where the link between the protein target and the assay is direct
    df_assays = df_assays[df_assays.confidence_score==9]

    df_activities = df_activities[df_activities.assay_chembl_id.isin(df_assays.assay_chembl_id)]

    # keep only the columns you need
    df_res = df_activities[['assay_description','molecule_chembl_id','molecule_pref_name', 'canonical_smiles','pchembl_value',\
                'standard_type','standard_relation','standard_value','standard_units','target_pref_name',
                'target_organism']]

    # remove canonical smiles with NaN
    curated_df = df_res.drop(df_res[df_res['canonical_smiles'].isna()].index, axis=0)

    return curated_df
