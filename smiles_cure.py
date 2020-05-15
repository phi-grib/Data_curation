"""
    Code for curating SMILES. To be used from a list/pandas dataframe.
    If SMILES are in a text file, first it will be processed either as a python list or a pandas dataframe.

    In principle it is written to work with CII and CR databases, but eventually it should be extended for 
    all phi projects if needed.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 15/05/2020, 16:29 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from Typing import Optional, Union

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps


class Curator(object):
    """
        Initializes the class with a SMILES input and applies standardization and other functions to curate the data
    """

    def __init__(self, smiles_dataframe: pd.DataFrame, smiles_field: str):
        """
            Initialized class with the smiles dataframe

            :param smiles_dataframe: Dataframe containing the input SMILES
            :param smiles_field: Column name in the DF containing SMILES
        """

        self.smi_df = smiles_dataframe
        self.smi_field = smiles_field
        
    def get_
