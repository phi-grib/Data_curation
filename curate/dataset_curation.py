"""
    Code that handles the input file(s) that need to be curated.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 02/02/2021, 17:32 PM
"""

import json
import numpy as np
import pandas as pd
import pickle
import sys

from rdkit import Chem
from rdkit.Chem import PandasTools
from typing import Optional, Union

from curate.parameters import Parameters
from curate.util import utils

class DataCuration(object):

    """
        Main class for handling inputs from users.
        Internally, it will curate the data in the following way:
            - Structure normalization.
            - Substance type identification from structure
            - Dataset selection (optional)
            - Dataset resampling (optional)
    """
    
    def __init__(self, data_input: Union[pd.DataFrame,str], molecule_identifier: str, structure_column: str, output_dir: str, 
                        endpoint: str, metadata: Union[list,str], separator: str = None, remove_problematic: bool = None, outfile_type: str = None):
        """
            Initialize class getting substance types for structure curation.
        """

        self.substance_types = self.get_substance_types()
        self.separator = separator
        self.input_data = self.process_input(data_input)
        self.identifier = molecule_identifier
        self.structure_column = structure_column
        self.output_dir = output_dir
        self.endpoint = endpoint
        self.metadata = metadata
        self.remove_problematic = remove_problematic
        self.outfile_type = outfile_type

        ## Stores a copy of the input data in the curation endpoint directory
        self.write_input_data()
        
        ## Stores parameters in curation_parameters.yaml file
        self.param = Parameters()

        param_string = {'data_input': data_input, 
                        'molecule_identifier': self.identifier,
                        'structure_column': self.structure_column,
                        'endpoint':self.endpoint,
                        'metadata':self.metadata,
                        'separator':self.separator,
                        'remove_problematic':self.remove_problematic,
                        'outfile_type':self.outfile_type}

        param_string = json.dumps(param_string)
        
        success, message = self.param.delta_curation(endpoint, param_string, iformat='JSONS')
        
        if not success:
            sys.stderr.write('Unable to load curation parameters. {}. Aborting...\n'.format(message))
            sys.exit(1)

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
                i_data = pd.read_excel(data_input, engine='openpyxl')
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

    def get_output_file(self, outfile_type: str = None, smiles_column: str = None, data: pd.DataFrame = None, outfile_name: str = None):
        """
            Saves the curated data into a specific file format.
            Requires output file name and type of file (Excel, CSV, TSV, sdf)

            :param outfile_type: format for the output file {.xlsx, .csv, .tsv, .sdf}
            :param smiles_column: SMILES column in the dataframe to be processed
            :param data: Dataframe to be written
            :param outfile_name: name of the output file

            :return output_file:
        """

        if outfile_type is None and self.outfile_type:
            outfile_type = self.outfile_type
        
        if data is None:
            data_copy = self.curated_data.copy()
        else:
            data_copy = data.copy()

        if outfile_name is None:
            outfile_name = 'curated_data'

        outfile_full_path = '/'.join([self.output_dir,outfile_name])

        if 'sdf' in outfile_type.lower():
            self.write_sdf(data_copy, outfile_full_path, smiles_column)
        elif 'xlsx' in outfile_type.lower() or 'excel' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'xlsx'])
            data_copy.to_excel(output_name_format)
        elif 'csv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'csv'])
            data_copy.to_csv(output_name_format, sep=',')
        elif 'tsv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'tsv'])
            data_copy.to_csv(output_name_format, sep='\t')
        elif 'json' in outfile_type.lower():
            output_name_format = '.'.join([outfile_full_path,'json'])
            data_copy.to_json(path_or_buf = output_name_format, orient = 'index')

    def save_output_header(self):
        """
            Stores in a pickle a header from the curated data so it can be retrieved
            in the GUI
        """

        head_pickle_full_path = '/'.join([self.output_dir,'curated_data_head.pkl'])
        if self.metadata:
            cols = [self.identifier,self.structure_column,'structure_curated','substance_type_name']
            cols.extend(self.metadata)
        else:
            cols = [self.identifier,self.structure_column,'structure_curated','substance_type_name']

        output_header = self.curated_data[cols].head(10)
        output_header.to_pickle(head_pickle_full_path)
    
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
            self.add_mol_column_to_df(data, smiles_column)

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

        if self.remove_problematic:
            smiles_stats_dict = {'Total SMILES':len(smiles_dataframe.index),
                                'Processed SMILES':len(self.curated_data.index), 
                                'Unable to process':len(self.problematic_structures.index)}
        else:
            smiles_stats_dict = {'Total SMILES':len(smiles_dataframe.index),
                                'Processed SMILES':len(self.curated_data.index)}

        return smiles_stats_dict

    def get_total_of_smiles_per_type_of_substance(self, smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
            This function returns the amount of substance types found of each kind
            after curating the SMILES.

            :param smiles_dataframe:

            :return subs_count:
        """

        subs_count = smiles_dataframe.groupby('substance_type_name')['substance_type_name'].count().to_dict()

        return subs_count

    def calculate_data_stats(self, dataframe: pd.DataFrame):
        """
            Counts how many substances have been processed, how many haven't and the different
            types of substances calculated.

            :param dataframe: curated data dataframe
        """

        data_stats = self.get_number_of_processed_vs_unprocessed(dataframe)
        subs_types_stats = self.get_total_of_smiles_per_type_of_substance(dataframe)

        general_stats = {}
        general_stats['curation_stats'] = data_stats
        general_stats['substance_types'] = subs_types_stats
        
        stats_file = utils.curation_tree_path('/'.join([self.endpoint,'statistics.pkl']))
        
        with open(stats_file, 'wb') as fo:
            pickle.dump(general_stats, fo)

    def curate_data(self) -> pd.DataFrame:
        """
            Check SMILES column to get a curated SMILES and the type of substance.

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

        if self.remove_problematic:
            self.remove_problematic_structures(curated_data)
            self.get_output_file(outfile_type='xlsx', data=self.problematic_structures, outfile_name='Problematic_structures_removed')
        else:
            self.curated_data = curated_data

        self.calculate_data_stats(self.curated_data)
        self.save_output_header()
        
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