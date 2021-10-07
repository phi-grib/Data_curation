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
        
        self.x_matrix = self.get_x_matrix()
        self.y_matrix = pd.DataFrame()

    def get_x_matrix(self) -> pd.DataFrame:
        """
            Searches the columns that have float values and passes them to the X matrix.

            :return x_matrix:
        """

        htt_data = self.input_data.copy()
        htt_data.drop([self.identifier, self.structure_column], axis=1, inplace=True)
        x_mat_cols = []

        for column in htt_data.columns:
            if htt_data[column].dtype == 'float64':
                x_mat_cols.append(column)
        
        x_matrix = htt_data[x_mat_cols]

        return x_matrix
    
    def get_y_matrix(self) -> pd.DataFrame:
        """
        """
        pass

    def chemical_curation(self):
        """
        """

        self.input_data = self.input_data[[self.identifier, self.structure_column]]
        self.curate_data()
    
    def write_chem_output_and_x_matrix(self):
        """
        """

        self.write_output_curation_data()

        outfile_x_matrix_path = '/'.join([self.output_dir,'x_matrix.pkl'])
        self.x_matrix.to_pickle(outfile_x_matrix_path)