"""
    Code that calculates fingerprints from compound structures and calculates similarities
    using Fingerprints.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/11/2022, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import DataStructs

from .chemical_description import Description
from curate.util import get_logger

LOG = get_logger(__name__)

class Similarity(object):
    """
        Calculates the similarity of all the compounds in the input dataframe.
    """

    def __init__(self, dataframe: pd.DataFrame, molecule_id: str ,smiles_column: str, activity: str, similarity_threshold: float) -> None:
        
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