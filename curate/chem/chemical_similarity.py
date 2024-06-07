"""
    Code that calculates fingerprints from compound structures and calculates similarities
    using Fingerprints.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/11/2022, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs

from .chemical_description import Description
from curate.util import get_logger

from typing import Union

LOG = get_logger(__name__)

class Similarity(object):
    """
        Calculates the similarity of all the compounds in the input dataframe.
    """

    def __init__(self, dataframe: pd.DataFrame, molecule_id: str ,smiles_column: str, similarity_threshold: float, activity:str = None) -> None:
        
        self.compound_dataframe = dataframe.copy()
        self.molecule_id = molecule_id
        self.activity = activity
        self.smiles_column = smiles_column
        self.threshold = similarity_threshold

    def describe_compound_dataframe(self) -> pd.DataFrame:
        """
            Initializes the descriptor object and obtains RDKit descriptors and FPs

            :return compound_dataframe_described:
        """

        description = Description(self.compound_dataframe, self.molecule_id, self.smiles_column)
        
        description.add_descriptors_and_fingerprints()

        self.compound_dataframe = description.compound_dataframe
        self.descriptor_dataframe = description.descriptor_dataframe
        self.morgan_dataframe = description.morgan_dataframe

    def get_similarities_between_all_compounds(self) -> pd.DataFrame:
        """
            This function calculates the similarity between compounds using Tanimoto index.

            :return comp_df:
        """

        comparison_dict = {'name':[],'name_structure':[], 'target_name':[],'target_structure':[],'activity':[],'target_activity':[],'similarity':[]}

        index_to_avoid = []
        for i, row in self.compound_dataframe.iterrows():
            name = row[self.molecule_id]
            struc = row['canon_smiles']
            activity = row[self.activity]
            fps = row['fps']
            index_to_avoid.append(i)
            fps_to_compare = self.compound_dataframe.loc[~self.compound_dataframe.index.isin(index_to_avoid),'fps'].values
            index_to_consider = self.compound_dataframe.loc[~self.compound_dataframe.index.isin(index_to_avoid),:].index
            try:
                s = DataStructs.BulkTanimotoSimilarity(fps, fps_to_compare)
            except TypeError:
                LOG.error(self.compound_dataframe[self.compound_dataframe.index.isin([i])])
                raise
            for sim,idx in zip(s,index_to_consider):
                comparison_dict['name'].append(name)
                comparison_dict['name_structure'].append(struc)
                comparison_dict['target_name'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),self.molecule_id].values[0])
                comparison_dict['target_structure'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),'canon_smiles'].values[0])
                comparison_dict['activity'].append(activity)
                comparison_dict['target_activity'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),self.activity].values[0])
                comparison_dict['similarity'].append(sim)
        
        comp_df = pd.DataFrame(data=comparison_dict)
        comp_df = comp_df.sort_values(by=['name','similarity'], ascending=False)
        comp_df.drop_duplicates(inplace=True)

        return comp_df

    def add_similarity_class_to_dataframe(self) -> pd.DataFrame:
        """
            Creates a similarity class based on a threshold and concatenates activity and similarity so it can be passed to 
            stratify argument of train test split.
            Class concatenation based on https://stackoverflow.com/questions/45516424/sklearn-train-test-split-on-pandas-stratify-by-multiple-columns
        """

        comp_df = self.prepare_dataframe_for_similarity()

        comp_df.loc[comp_df['similarity'] >= self.threshold, 'similar'] = 1
        comp_df.loc[comp_df['similarity'] < self.threshold, 'similar'] = 0

        similars = comp_df.loc[comp_df['similar'] == 1, 'name'].unique()
        similars = np.append(similars, comp_df.loc[comp_df['similar'] == 1, 'target_name'].unique())

        self.compound_dataframe.loc[self.compound_dataframe[self.molecule_id].isin(similars), 'similar'] = 1
        self.compound_dataframe.loc[~self.compound_dataframe[self.molecule_id].isin(similars), 'similar'] = 0
        self.compound_dataframe['similar'] = self.compound_dataframe['similar'].astype(int)

        self.compound_dataframe['new_class'] = self.compound_dataframe[self.activity].astype(str) + "_" + self.compound_dataframe['similar'].astype(str)

    def calculate_tanimoto_similarity_manual_index(self, fp1: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray], 
                                                   fp2: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray]) -> float:
        """
            Calculate Tanimoto similarity between two Morgan fingerprints given as bit vectors.
            
            :param fp1: First Morgan fingerprint as an array of bits.
            :param fp2: Second Morgan fingerprint as an array of bits.
                
            :returns similarity: Tanimoto similarity between the two fingerprints.
        """
        
        # Convert the fingerprints to sets of indexes where the fingerprint is 1
        fp1_set = set(i for i, bit in enumerate(fp1) if bit == 1)
        fp2_set = set(i for i, bit in enumerate(fp2) if bit == 1)

        # Calculate the intersection and union sizes
        intersection_size = len(fp1_set.intersection(fp2_set))
        union_size = len(fp1_set.union(fp2_set))

        # Calculate Tanimoto similarity
        similarity = intersection_size / union_size if union_size != 0 else 0.0
        
        return similarity

    def calculate_similarity_between_two_datasets(self, df1: pd.DataFrame, df2: pd.DataFrame, smiles1: str, smiles2:str, 
                                                  molid1:str, molid2:str, activity1: str, activity2: str, descriptors:str, top_feat= None):
        """
            Calculates the tanimoto similarity between two datasets using Morgan FPs obtained from chemical_description.py.

            :param df1: dataframe 1
            :param df2: dataframe 2
            :param smiles1: SMILES column name of dataframe 1
            :param smiles2: SMILES column name of dataframe 2
            :param molid1: Molecule ID column name of dataframe 1
            :param molid2: Molecule ID column name of dataframe 2
            :param descriptors: descriptors column name to use for calculating similarity
            :param top_feat: optional feature. If True, it means the dataframes include the top_10 Morgan FPs 
                            that most contribute to the modelling of the dataframe previously done.

            :return similarity_df:
        """
        # Initialize a list to store the similarities
        similarity_data = []

        # Iterate over each compound in the first dataframe
        for idx1, row1 in df1.iterrows():
            # Get the molecule object from SMILES notation
            mol1 = Chem.MolFromSmiles(row1[smiles1])
            act1 = row1[activity1]

            if top_feat:
                fp1 = row1['top_10_features']
            else:
                fp1 = row1[descriptors]

            # Iterate over each compound in the second dataframe
            for idx2, row2 in df2.iterrows():
                # Get the molecule object from SMILES notation
                mol2 = Chem.MolFromSmiles(row2[smiles2])
                act2 = row2[activity2]

                if top_feat:
                    fp2 = row2['top_10_features']
                else:
                    fp2 = row2[descriptors]
                
                if top_feat:
                    similarity = self.calculate_tanimoto_similarity_manual_index(fp1,fp2)
                else:
                    # Calculate the Tanimoto similarity between the fingerprints
                    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                
                # Append the similarity data to the list
                similarity_data.append([row1[molid1], act1, row2[molid2], act2, similarity, mol1, mol2])

        # Create a new dataframe from the similarity data
        similarity_df = pd.DataFrame(similarity_data, columns=['Compound_ID_1','Activity_1', 'Compound_ID_2', 'Activity_2','Similarity', 'Structure_1', 'Structure_2'])
        
        return similarity_df