"""
    Code for extracting ChEMBL datasets using
    a ChEMBL ID.

    ChEMBL ID has to point to a target or protein (e.g: CHEMBL230 is Cox2).

    TODO:Might incoporate other arguments like species (e.g: Homo sapiens)

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/04/2021, 11:13 PM
"""

import numpy as np
import pandas as pd

from chembl_webresource_client.new_client import new_client

from curate.util import get_logger
LOG = get_logger(__name__)

def get_dataframe_from_target(chembl_id: str) -> pd.DataFrame:
    """
        Creates a client to ChEMBL and uses chembl_id to extract
        all the compounds related to that target.
        Returns a dataframe with all the information.

        :param target_name: string with a valid chembl_id

        :return curated_df: dataframe containing the results
    """

    # Create an activity query
    activities = new_client.activity
    
    # Make chembl_id uppercase
    chembl_id_upper = chembl_id.upper()
    
    # Select only activities with a pchembl_value (-log(IC50, Ki, Kd, EC50...).
    # We also use the chembl flags to remove the duplicates and the records where there is a validity comment
    response = activities.filter(target_chembl_id=chembl_id_upper, pchembl_value__isnull=False,\
                                potential_duplicate=False, data_validity_comment__isnull=True )
    
    if not response:
        curated_df = None
    else:
        curated_df = get_dataframe_from_response(response)

        if curated_df.empty:
            curated_df = None
        else:
            curated_df.loc[:,'chembl_id'] = chembl_id_upper

    return curated_df

def get_dataframe_from_response(response: list) -> pd.DataFrame:
    """
        Checks response when True and returns a curated dataframe with the information from ChEMBL.

        :param response: list of dictionaries with the information from ChEMBL
        :return curated_df: dataframe with all the information
    """

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

def process_list_of_chembl_ids(raw_df: pd.DataFrame) -> np.ndarray:
    """
        Gets the CHEMBLIDs from the column where they're stored and returns them as a numpy array.

        :param raw_df: dataframe with the CHEMBLIDs
        :return chembl_id_list: list of CHEMBLIDs
    """

    chembl_id_column = [name for name in raw_df.columns if 'chembl' in name.lower()]
    chembl_id_list = raw_df[chembl_id_column].values.flatten()

    return chembl_id_list

def concatenate_dataframes_from_different_chembl_ids(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
        Concatenates all the dataframes from different CHEMBLIDs.

        :param raw_df: dataframe with the CHEMBLIDs
        :return chembl_targets_concat: dataframe with all the information
    """

    chembl_id_list = process_list_of_chembl_ids(raw_df)
    
    chembl_targets_concat = pd.DataFrame()
    not_valid_ids = []
    for chembl_id in chembl_id_list:
        df_to_add = get_dataframe_from_target(chembl_id)
        if isinstance(df_to_add, pd.DataFrame):
            chembl_targets_concat = pd.concat([chembl_targets_concat, df_to_add], ignore_index=True)
        else:
            not_valid_ids.append(chembl_id)
    
    if not_valid_ids:
        warning = "WARNING: The following CHEMBLIDs were not processed: {}\nPlease, check that your IDs point to a target/protein that has enough compounds assayed\n".format(not_valid_ids)
    else:
        warning = 'All IDs were processed successfully\n'

    LOG.warning(warning)

    return chembl_targets_concat, warning