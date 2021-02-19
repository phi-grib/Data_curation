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
    
    def __init__(self, data_input: Union[pd.DataFrame,str], molecule_identifier: str, structure_column: str):
        """
            Initialize class getting substance types for structure curation.
        """

        self.substance_types = self.get_substance_types()
        self.input_data = self.process_input(data_input)
        self.identifier = molecule_identifier
        self.structure_column = structure_column

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
    
    def get_metadata(self) -> list:
        """
            Checks the columns of the input data and returns all that are 
            not the SMILES. Everything but the structure is considered metadata.

            :return metadata: returns a list with the metadata, which are the 
            column names of the input data without the SMILES
        """

        metadata = self.input_data.loc[:, self.input_data.columns != self.structure_column].columns

        return metadata

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
        cur_data = self.prepare_data_for_sdf(copy=True)
        
        PandasTools.WriteSDF(cur_data, output_name_format, molColName='ROMol', properties=list(cur_data.columns), idName=self.identifier)

    def prepare_data_for_sdf(self, copy: bool = False) -> Optional[pd.DataFrame]:
        """
            Prepares the data to be converted to sdf.
            If copy, it copies the dataframe so it's not overwritten with new columns before being processed as sdf.
            Else, it directly uses self.curated_data. This option is used mostly in jupyter or CLI mode to keep the new columns
            in the python object so it can be manipulated directly in the backend.

            :param copy: boolean accepting True or False

            :return cur_data: dataframe with new columns added before being converted into sdf.
        """

        if copy:
            copy_curated_data = self.curated_data.copy()
            cur_data = self.add_mol_column_to_df(copy_curated_data)
            return cur_data
        else:
            self.add_mol_column_to_df(self.curated_data)

    def add_mol_column_to_df(self, data: pd.DataFrame) -> pd.DataFrame:
        """
            Applies PandasTools functionalities to process the structure into a valid format for the sdf transformation.

            :param data: dataframe to be modified

            :return data: modified data
            :return no_mol: data that hasn't been modified
        """

        PandasTools.AddMoleculeColumnToFrame(data,'structure_curated')
        no_mol = data[data['ROMol'].isna()]
        data.drop(no_mol.index, axis=0, inplace=True)
        data.loc[:,'ROMol'] = [Chem.AddHs(x) for x in data['ROMol'].values.tolist()]
        
        if no_mol.empty is False:
            no_mol.to_excel('Non_processed_molecules.xlsx')

        return data

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
    
    def curate_data(self, remove_problematic: bool = None) -> pd.DataFrame:
        """
            Check SMILES column to get a curated SMILES and the type of substance.

            :param remove_problematic: it allows the user to get rid of problematic structures for QSAR modelling. 

            :return curated_data: dataframe containing the curated information
        """
        
        from curate.chem import structure_curation as cur
        data_cur = cur.Curator()

        curated_data = self.input_data.copy()
        
        for i, row in curated_data.iterrows():
            smi = row[self.structure_column]
            data_cur.get_rdkit_mol(smi)
            sub_type, san_smi = data_cur.filter_smiles()
            curated_data.ix[i,'structure_curated'] = san_smi
            curated_data.ix[i,'substance_type_name'] = sub_type
        
        if remove_problematic:
            self.remove_problematic_structures(curated_data)
            self.problematic_structures.to_excel('Problematic_structures_removed.xlsx')
        else:
            self.curated_data = curated_data
    
    def remove_problematic_structures(self, data: pd.DataFrame = None) -> pd.DataFrame:
        """
            Remove problematic structures from main dataset.
            Returns cleaned dataset and problematic structures a part.

            :param data: input data to be cleaned

            :return data_cleaned: data without problematic structures
            :return problematic_structures: data with the problematic structures

            TODO: add option to select specific substance types to be removed.
        """

        if data.empty:
            data = self.input_data
        
        problem_struc_list =  ['organometallic', 'no_sanitizable', 'inorganic_salt', 
                              'inorganic', 'inorganic_metal', 'no_sanitizable_organic',
                              'no_sanitizable_inorganic', 'no_sanitizable_organometallic']

        self.curated_data = data.loc[~data['substance_type_name'].isin(problem_struc_list)]
        self.problematic_structures = data.loc[data['substance_type_name'].isin(problem_struc_list)]

    def split_dataset(self, train_proportion: float, test_proportion: float, activity_field: str):
        """
            Split the curated dataset into training and test sets using the proportion provided by the user.

            :param train_proportion: training set proportion of the main dataset
            :param test_proportion: test set proportion of the main dataset
            :param activity_field: column name that inclues the activity
        """

        from curate.data_handler import dataset_selection as datasel
        
        curated_data_object = datasel.Selection(self.curated_data, train_proportion, test_proportion, activity_field)
        self.train, self.test = curated_data_object.split_main_dataset()
    
    def correct_imbalance(self, dataset: Union[pd.DataFrame,str], activity_field: str, imbalance_algorithm: str) -> Optional[pd.DataFrame]:
        """
            Applies imbalance correction algorithms to the desired dataset.
            Input dataset can be the training set, the test set or an external pandas dataframe.
            Activity field and name of algorithm must be provided.

            :param dataset: it accepts two kins of inputs: a string being either 'train' or 'test' or a pandas dataframe.
            :param activity_field: string containing the name of the column with the activity data.
            :param imbalance_algorithm: string containing the name of the imbalance algorithm to use. {oversampling, subsampling, smoteenn, smotetomek}

            :return corrected_dataset: if the input is a dataframe, it returns the corrected version of the dataset.
        """

        from curate.data_handler import dataset_imbalance_correction as imb

        if 'train' in dataset:
            imb_object = imb.ImbalanceData(self.train, activity_field, imbalance_algorithm)
            self.train_corrected = imb_object.imbalance_correction()
        elif 'test' in dataset:
            imb_object = imb.ImbalanceData(self.test, activity_field, imbalance_algorithm)
            self.test_corrected = imb_object.imbalance_correction()
        else:
            imb_object = imb.ImbalanceData(dataset, activity_field, imbalance_algorithm)
            corrected_dataset = imb_object.imbalance_correction()
            return corrected_dataset