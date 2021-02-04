"""
    Code that handles the input file(s) that need to be curated.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 02/02/2021, 17:32 PM
"""

import numpy as np
import pandas as pd
import sys

from rdkit import Chem
from rdkit.Chem import PandasTools
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
            elif data_input.endswith('.csv'):
                i_data = pd.read_csv(data_input, sep=',')
            elif data_input.endswith('.tsv'):
                i_data = pd.read_csv(data_input, sep='\t')
            elif data_input.endswith('.sdf'):
                i_data = PandasTools.LoadSDF(data_input)
            else:
                sys.stderr.write('Please provide a file with a valid format (xlsx, csv, tsv, sdf)\n')
        
        return i_data
    
    def get_output_file(self, outfile_name: str, outfile_type: str):
        """
            Saves the curated data into a specific file format.
            Requires output file name and type of file (Excel, CSV, TSV, sdf)

            :param outfile_name: name of the output file
            :param outfile_type: format for the output file {.xlsx, .csv, .tsv, .sdf}

            :return output_file:
        """

        if 'sdf' in outfile_type.lower():
            self.write_sdf(outfile_name)
        elif 'xlsx' in outfile_type.lower() or 'excel' in outfile_type.lower():
            output_name_format = '.'.join([outfile_name.split('.')[0],'xlsx'])
            self.curated_data.to_excel(output_name_format)
        elif 'csv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_name.split('.')[0],'csv'])
            self.curated_data.to_csv(output_name_format, sep=',')
        elif 'tsv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_name.split('.')[0],'tsv'])
            self.curated_data.to_csv(output_name_format, sep='\t')

    def write_sdf(self, outfile_name: str):
        """
            Prepares curated data to be converted into sdf file using
            PandasTools. Returns non processed molecules in excel format.

            :param outfile_name: output file name
        """

        output_name_format = '.'.join([outfile_name.split('.')[0],'sdf'])
        copy_curated_data = self.curated_data.copy()

        PandasTools.AddMoleculeColumnToFrame(copy_curated_data,'structure_curated')
        no_mol = copy_curated_data[copy_curated_data['ROMol'].isna()]
        copy_curated_data.drop(no_mol.index, axis=0, inplace=True)
        copy_curated_data['ROMol'] = [Chem.AddHs(x) for x in copy_curated_data['ROMol'].values.tolist()]

        PandasTools.WriteSDF(copy_curated_data, output_name_format, molColName='ROMol', properties=list(copy_curated_data.columns), idName='name')

        if no_mol.empty is False:
            no_mol.to_excel('Non_processed_molecules.xlsx')

    def get_substance_types(self) -> pd.DataFrame:
        """
            Uses the dictionary of substance type to generate a pandas dataframe.

            :return substance_types:
        """

        substance_types = {'id':[1,2,3,4,5,6,7,8,9,10,11], 'type':['organic', 'organic_salt','organometallic',
                            'peptide','inorganic', 'inorganic_metal','inorganic_salt','no_sanitizable',
                            'no_sanitizable_organic','no_sanitizable_inorganic','no_sanitizable_organometallic']}
        substance_types = pd.DataFrame(substance_types, columns=['id','type'])
        
        return substance_types
    
    def curate_data(self, structure_column: str, remove_problematic: bool = None) -> pd.DataFrame:
        """
            Check SMILES column to get a curated SMILES and the type of substance.

            :param structure_column: string with the column name that contains the SMILES
            :param remove_problematic: it allows the user to get rid of problematic structures for QSAR modelling. 

            :return curated_data: dataframe containing the curated information
        """
        
        from curate import structure_curation as cur
        data_cur = cur.Curator()

        curated_data = self.input_data.copy()

        for i, row in curated_data.iterrows():
            smi = row[structure_column]
            data_cur.get_rdkit_mol(smi)
            sub_type, san_smi = data_cur.filter_smiles()
            sub_type_id = self.substance_types.loc[self.substance_types['type'] == sub_type, 'id'].values[0]
            curated_data.ix[i,'structure_curated'] = san_smi
            curated_data.ix[i,'substance_type_id'] = sub_type_id
            curated_data.ix[i,'substance_type_name'] = sub_type
        
        if remove_problematic:
            curated_data, problematic_structures = self.remove_problematic_structures(curated_data)
            self.problematic_structures = problematic_structures
        
        self.curated_data = curated_data
    
    def remove_problematic_structures(self, data: pd.DataFrame) -> pd.DataFrame:
        """
            Remove problematic structures from main dataset.
            Returns cleaned dataset and problematic structures a part.

            :param data: input data to be cleaned

            :return data_cleaned: data without problematic structures
            :return problematic_structures: data with the problematic structures

            TODO: add option to select specific substance types to be removed.
        """

        problem_struc_list = ['organometallic','no_sanitizable', 'inorganic_salt', 
                                'inorganic','inorganic_metal','no_sanitizable_organic',
                                'no_sanitizable_inorganic','no_sanitizable_organometallic']

        data_cleaned = data.loc[~data['type'].isin(problem_struc_list)]
        problematic_structures = data.loc[data['type'].isin(problem_struc_list)]

        return data_cleaned, problematic_structures

    def split_dataset(self):
        pass