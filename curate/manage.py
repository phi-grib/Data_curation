"""
    Managing functions for Data curation CLI.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/02/2021, 17:49 PM
"""

import os
import pathlib
import sys
import time

def get_metadata(data: pd.DataFrame, structure_colname: str) -> list:
    """
        Checks the columns of the input data and returns all that are 
        not the SMILES. Everything but the structure is considered metadata.
        
        :param data: dataframe to be used
        :param structure_colname: name of the column containing the SMILES

        :return metadata: returns a list with the metadata, which are the 
        column names of the input data without the SMILES
    """

    metadata = data.loc[:, data.columns != structure_colname].columns

    return metadata

def convert_to_json(data: pd.DataFrame, filename: str = None) -> str:
    """
        Converts into JSON the input dataframe.
        If filename is provided, json is written as a file.

        :param data: dataframe to be converted.
        :param filename: filename where to dump the JSON.

        :return json_data: JSON string from pandas dataframe
    """

    json_data = data.to_json(path_or_buf= filename, orient='index')

    return json_data

def set_curation_repository(path: str = None):
    """
        Set the path to the curation repository

        :param path: string indicating the path of the curation repository
    """

    pass

#### Functions to represent simple statistics after the curation

def get_number_of_processed_vs_unprocessed(smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
        This function returns a dataframe with the number of total SMILES, the ones
        that have been processed by the code and the ones that haven't.

        Is it expected to receive an input from the API using the output of the previous
        function get_list_of_filtered_SMILES_and_substance_type()

        :param smiles_dataframe:

        :return smiles_stats_df:
    """

    unprocessed_smiles = smiles_dataframe[smiles_dataframe['structure_curated'].isna()]
    processed_smiles = smiles_dataframe[~smiles_dataframe['structure_curated'].isna()]

    smiles_stats_dict = {'Total SMILES':len(smiles_dataframe.index), 
                         'Processed SMILES':len(processed_smiles.index), 
                         'Unable to process':len(unprocessed_smiles.index)}
    
    smiles_stats_df = pd.DataFrame(data=smiles_stats_dict)

    return smiles_stats_df

def get_total_of_smiles_per_type_of_substance(smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
        This function returns the amount of substance types found of each kind
        after curating the SMILES.

        :param smiles_dataframe:

        :return subs_count:
    """

    subs_count = smiles_dataframe.groupby('substance_type_name')['substance_type_name'].count()

    return subs_count