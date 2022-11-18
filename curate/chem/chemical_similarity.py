"""
    Code that calculates fingerprints from compound structures and calculates similarities
    using Fingerprints.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/11/2022, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from typing import Optional, Union

class Similarity(object):
    """
        Calculates the similarity of all the compounds in the input dataframe.
    """

    def __init__(self, dataframe: pd.DataFrame, smiles_column: str) -> None:
        
        self.compound_dataframe = dataframe
        self.smiles_column = smiles_column

        self.comparison_dict = {'name':[],'name_structure':[], 'target_name':[],'target_structure':[],'activity':[],'target_activity':[],'similarity':[]}

    def prepare_dataframe_for_similarity(self):
        """
        """

        # Adds canonical SMILES
        self.compound_dataframe.loc[:, 'canon_smiles'] = self.compound_dataframe.loc[:,self.smiles_column].apply(lambda x: self.get_canonical_smiles(x))

        # Adds mol object from canonical SMILES
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'canon_smiles'].apply(lambda x: Chem.MolFromSmiles(x))

        # Adds Fingerprints
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: FingerprintMols.FingerprintMol(x))

        # Get Dataframe with similarities
        comparison_dict = self.get_similarities_between_all_compounds()

        return comparison_dict

    def get_canonical_smiles(self, smiles: str) -> Optional[str]:
        """
            Gets the canonical SMILES from the structure in the dataframe

            :return canonical_smiles: canonical smiles
        """

        try:
            canonical_smiles = Chem.CanonSmiles(smiles)
        except:
            canonical_smiles = None

        return canonical_smiles

    def get_similarities_between_all_compounds(self) -> pd.DataFrame:
        """
        """

        index_to_avoid = []
        for i, row in self.compound_dataframe.iterrows():
            name = row['name']
            struc = row['canon_smiles']
            activity = row['activity']
            fps = row['fps']
            index_to_avoid.append(i)
            fps_to_compare = self.compound_dataframe.loc[~self.compound_dataframe.index.isin(index_to_avoid),'fps'].values
            index_to_consider = self.compound_dataframe.loc[~self.compound_dataframe.index.isin(index_to_avoid),:].index
            try:
                s = DataStructs.BulkTanimotoSimilarity(fps, fps_to_compare)
            except TypeError:
                print(self.compound_dataframe[self.compound_dataframe.index.isin([i])])
                raise
            for sim,idx in zip(s,index_to_consider):
                self.comparison_dict['name'].append(name)
                self.comparison_dict['name_structure'].append(struc)
                self.comparison_dict['target_name'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),'name'].values[0])
                self.comparison_dict['target_structure'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),'canon_smiles'].values[0])
                self.comparison_dict['activity'].append(activity)
                self.comparison_dict['target_activity'].append(self.compound_dataframe.loc[self.compound_dataframe.index.isin([idx]),'activity'].values[0])
                self.comparison_dict['similarity'].append(sim)
        
        comp_df = pd.DataFrame(data=self.comparison_dict)
        comp_df = comp_df.sort_values(by=['name','similarity'], ascending=False)

        return comp_df