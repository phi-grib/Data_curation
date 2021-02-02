"""
    Code that handles the input file(s) that need to be curated.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 02/02/2021, 17:32 PM
"""

import numpy as np
import pandas as pd
import sys

from typing import Optional, Union, Tuple

class DataCuration(object):
    """
        Main class for handling inputs from users.
        Internally, it will curate the data in the following way:
            - Structure normalization.
            - Substance type identification from structure
            - Dataset selection (optional)
            - Dataset resampling (optional)
        
        TODO:More features will be implemented
    """
    
    def __init__(self, data_input: Union[pd.DataFrame,str]):
        """
            Initialize class getting substance types for structure curation.
        """

        self.substance_types = self.get_substance_types()
        self.input_data = self.process_input(data_input)

    def process_input(self, data_input: Union[pd.DataFrame,str]) -> pd.DataFrame:
        """
            Checks if input is an Excel file and converts it into pandas dataframe.
            If it already is a pandas dataframe, nothing changes.

            :param data_input: it can be either a pandas dataframe or an excel file

            :return i_data: input data to be curated
        """

        if isinstance(data_input,pd.DataFrame):
            i_data = data_input
        elif isinstance(data_input,str):
            if data_input.endswith('.xlsx'):
                i_data = pd.read_excel(data_input)
            else:
                sys.stderr.write('Please provide with a file with a proper Excel format (.xlsx)')
        
        return i_data

    def get_substance_types(self) -> pd.DataFrame:
        """
            Uses the dictionary of substance type to generate a pandas dataframe.

            :return substance_types:
        """

        substance_types = {1:'organic', 2:'organic_salt', 3:'organometallic',  4:'peptide', 5:'inorganic',
                            6:'inorganic_metal',7:'inorganic_salt', 8:'no_sanitizable',9:'no_sanitizable_organic',
                            10:'no_sanitizable_inorganic',11:'no_sanitizable_organometallic'}
        substance_types = pd.DataFrame(substance_types, columns=['id','type'])
        
        return substance_types
    
    def curate_data(self, structure_column: str) -> pd.DataFrame:
        """
            Check SMILES column to get a curated SMILES and the type of substance.

            :param structure_column: string with the column name that contains the SMILES

            :return curated_data: dataframe containing the curated information
        """
        
        import .curation as cur
        data_cur = cur.Curator()

        curated_data = self.input_data.copy()

        for i, row in curated_data.iterrows():
            smi = row[structure_column]
            data_cur.get_rdkit_mol(smi)
            sub_type, san_smi = data_cur.filter_smiles()
            sub_type_id = threads_.loc[threads_['type'] == sub_type, 'id'].values[0]
            curated_data.ix[i,'structure_curated'] = san_smi
            curated_data.ix[i,'substance_type_id'] = sub_type_id
            curated_data.ix[i,'substance_type_name'] = sub_type
        
        return curated_data