"""
    Managing functions for Data curation CLI.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/02/2021, 17:49 PM
"""

import json
import os
import pathlib
import sys

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

def convert_to_json(data: pd.DataFrame) -> json:
    """
        Converts into JSON the input dataframe

        :param data:

        :return json_data:
    """

    json_data = data.to_json(orient='index')

    return json_data