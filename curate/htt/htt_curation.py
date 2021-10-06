"""
    Child class of dataset curation.
    Handles specific requirements of HTT input files.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 06/10/2021, 17:18 PM
"""

import pandas as pd

from curate.dataset_curation import DataCuration

from typing import Union

class htt_curation(DataCuration):

    """
        Child class of dataset curation.
        It uses the same input for inizialization and adds specific methods
        for handling HTT files.
    """

    def __init__(self, data_input: Union[pd.DataFrame,str], molecule_identifier: str, structure_column: str, output_dir: str, 
                endpoint: str, metadata: Union[list,str], separator: str = None, remove_problematic: bool = None):
        """
            Initializes class with main arguments of Data curation
        """
        
        super().__init__(data_input, molecule_identifier, structure_column, output_dir, 
                            endpoint, metadata, separator, remove_problematic)
        
            