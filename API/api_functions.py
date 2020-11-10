"""
    This code includes the functions that are called by the API to interact with
    SMILES curation code

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 10/11/2020, 17:16 PM
"""

import pandas as pd

from curate import curation as cur
from rdkit import Chem
from typing import Optional


#### Initialize object Curator

data_cur = cur.Curator()

##### Function to call from the API

def get_filtered_SMILES_and_substance_type(smiles: str) -> Optional[str]:
    """
        Function that converts SMILES to mol object of RDKit 
        and applies the filters. Returns substance type and a curated SMILES
        if possible

        :param smiles: SMILES string from the compound of interest

        :return substance_type: Substance type defined from the filters
        :return sanitized_smiles: Curated SMILES of the compound of interest. None if any error happened
    """

    data_cur.get_rdkit_mol(smiles)
    substance_type, sanitized_smiles = data_cur.filter_smiles()   

    return substance_type, sanitized_smiles

def get_list_of_filtered_SMILES_and_substance_type(list_of_smiles: pd.DataFrame, smiles_col_name: str) -> pd.DataFrame:
    """
        This function accepts a pandas dataframe as input with the column name where
        the SMILES are stored and returns a dataframe with the curated structures and the substance type

        :param list_of_smiles:
        :param smiles_col_name:

        :return smiles_df:
    """

    smiles_df = list_of_smiles.copy()

    for i, row in list_of_smiles.iterrows():
        smi = row[smiles_col_name]
        sub_type, san_smi = get_filtered_SMILES_and_substance_type(smi)
        smiles_df.ix[i,'structure_curated'] = san_smi
        smiles_df.ix[i,'substance_type_name'] = sub_type
    
    return smiles_df

#### Functions to represent simple statistics after the curation

def 