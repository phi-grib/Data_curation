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

from curate.idata import Idata
from curate.odata import Odata
from curate.util import utils

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
    
    def __init__(self, data_input: Union[pd.DataFrame,str], molecule_identifier: str, structure_column: str, output_dir: str, separator: str = None):
        """
            Initialize class getting substance types for structure curation.
        """

        self.substance_types = self.get_substance_types()
        self.separator = separator
        self.input_data = self.process_input(data_input)
        self.identifier = molecule_identifier
        self.structure_column = structure_column
        self.output_dir = output_dir
        
        ## Stores a copy of the input data in the curation endpoint directory
        self.write_input_data()
        
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
                if not self.separator:
                    self.separator = ','
                i_data = pd.read_csv(data_input, sep=self.separator)
            elif data_input.endswith('.tsv'):
                if not self.separator:
                    self.separator = '\t'
                i_data = pd.read_csv(data_input, sep=self.separator)
            elif data_input.endswith('.sdf'):
                i_data = PandasTools.LoadSDF(data_input)
            else:
                sys.stderr.write('Please provide a file with a valid format (xlsx, csv, tsv, sdf)\n')
                sys.exit()
                
        return i_data
    
    def write_input_data(self):
        """
            Uses the get_output_file function to write a copy of the input data in sdf
        """

        self.get_output_file(outfile_type='sdf', smiles_column=self.structure_column, data=self.input_data, outfile_name='input_data')

    def get_output_file(self, outfile_type: str, smiles_column: str = None, data: pd.DataFrame = None, outfile_name: str = None):
        """
            Saves the curated data into a specific file format.
            Requires output file name and type of file (Excel, CSV, TSV, sdf)

            :param outfile_type: format for the output file {.xlsx, .csv, .tsv, .sdf}
            :param smiles_column: SMILES column in the dataframe to be processed
            :param data: Dataframe to be written
            :param outfile_name: name of the output file

            :return output_file:
        """

        if data is None:
            data = self.curated_data.copy()

        if outfile_name is None:
            outfile_name = 'curated_data'

        outfile_full_path = '/'.join([self.output_dir,outfile_name])

        if 'sdf' in outfile_type.lower():
            self.write_sdf(data, outfile_full_path, smiles_column)
        elif 'xlsx' in outfile_type.lower() or 'excel' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'xlsx'])
            data.to_excel(output_name_format)
        elif 'csv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'csv'])
            data.to_csv(output_name_format, sep=',')
        elif 'tsv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'tsv'])
            data.to_csv(output_name_format, sep='\t')
        elif 'json' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'json'])
            data.to_json(path_or_buf = output_name_format, orient = 'index')
    
    def write_sdf(self, data: pd.DataFrame, outfile_name: str, smiles_column: str):
        """
            Prepares curated data to be converted into sdf file using
            PandasTools. Returns non processed molecules in excel format.

            :param data: Dataframe to be written
            :param smiles_column: SMILES column in the dataframe to be processed
            :param outfile_name: output file name
        """

        output_name_format = '.'.join([outfile_name,'sdf'])
        cur_data = self.prepare_data_for_sdf(data, smiles_column, copy=True)
        
        PandasTools.WriteSDF(cur_data, output_name_format, molColName='ROMol', properties=list(cur_data.columns), idName=self.identifier)

    def prepare_data_for_sdf(self, data: pd.DataFrame, smiles_column: str, copy: bool = False) -> Optional[pd.DataFrame]:
        """
            Prepares the data to be converted to sdf.
            If copy, it copies the dataframe so it's not overwritten with new columns before being processed as sdf.
            Else, it directly uses self.curated_data. This option is used mostly in jupyter or CLI mode to keep the new columns
            in the python object so it can be manipulated directly in the backend.

            :param data: Dataframe to be treated
            :param smiles_column: SMILES column in the dataframe to be processed
            :param copy: boolean accepting True or False

            :return cur_data: dataframe with new columns added before being converted into sdf.
        """

        if copy:
            cur_data = self.add_mol_column_to_df(data, smiles_column)
            return cur_data
        else:
            self.add_mol_column_to_df(data)

    def add_mol_column_to_df(self, data: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
        """
            Applies PandasTools functionalities to process the structure into a valid format for the sdf transformation.

            :param data: dataframe to be modified
            :param smiles_column: SMILES column in the dataframe to be processed

            :return data: modified data
            :return no_mol: data that hasn't been modified
        """

        PandasTools.AddMoleculeColumnToFrame(data, smiles_column)
        no_mol = data[data['ROMol'].isna()]
        data.drop(no_mol.index, axis=0, inplace=True)
        data.loc[:,'ROMol'] = [Chem.AddHs(x) for x in data['ROMol'].values.tolist()]
        
        if no_mol.empty is False:
            self.get_output_file(outfile_type='xlsx', data=no_mol, outfile_name='Non_processed_molecules')

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

    def get_number_of_processed_vs_unprocessed(self, smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
            This function returns a dataframe with the number of total SMILES, the ones
            that have been processed by the code and the ones that haven't.

            :param smiles_dataframe:

            :return smiles_stats_df:
        """

        smiles_stats_dict = {'SMILES':['Total SMILES','Processed SMILES', 'Unable to process'],
                             'Count':[len(smiles_dataframe.index), len(self.curated_data.index), len(self.problematic_structures.index)]}
        
        smiles_stats_df = pd.DataFrame(data=smiles_stats_dict)
    
        return smiles_stats_df

    def get_total_of_smiles_per_type_of_substance(self, smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
            This function returns the amount of substance types found of each kind
            after curating the SMILES.

            :param smiles_dataframe:

            :return subs_count:
        """

        subs_count = smiles_dataframe.groupby('substance_type_name')['substance_type_name'].count()
        
        return subs_count

    def calculate_data_stats(self, dataframe: pd.DataFrame):
        """
            Counts how many substances have been processed, how many haven't and the different
            types of substances calculated.

            :param dataframe: curated data dataframe
        """

        data_stats = self.get_number_of_processed_vs_unprocessed(dataframe)
        subs_types_stats = self.get_total_of_smiles_per_type_of_substance(dataframe)

        self.get_output_file(outfile_type='json', data=data_stats, outfile_name='curation_statistics')
        self.get_output_file(outfile_type='json', data=subs_types_stats, outfile_name='substance_type_statistics')

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
            curated_data.loc[curated_data.index == i,'structure_curated'] = san_smi
            curated_data.loc[curated_data.index == i,'substance_type_name'] = sub_type

        if remove_problematic:
            self.remove_problematic_structures(curated_data)
            self.calculate_data_stats(curated_data)
            self.get_output_file(outfile_type='xlsx', data=self.problematic_structures, outfile_name='Problematic_structures_removed')
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