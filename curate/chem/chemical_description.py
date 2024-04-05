"""
    Code that adds RDKit descriptors and Fingerprints into a pandas dataframe of chemical compounds

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 5/04/2024, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

from typing import Optional

from curate.util import get_logger

LOG = get_logger(__name__)

class Description(object):
    """
        Calculates RDKit descriptors and Fingerprints for all compounds in dataframe
    """

    def __init__(self, dataframe: pd.DataFrame, molecule_id: str, smiles_column: str):
        
        self.compound_dataframe = dataframe.copy()
        self.molecule_id = molecule_id
        self.smiles_column = smiles_column
    
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

    def get_rdkit_descriptors(self, dataframe: pd.DataFrame, id_col: str, mol_col: str) -> pd.DataFrame:
        """
            Gets a dataframe with RDKit descriptors

            :return df_descriptor:
        """

        names=[x[0] for x in Chem.Descriptors._descList]
        calculator=MoleculeDescriptors.MolecularDescriptorCalculator(names)
        desc = [calculator.CalcDescriptors(mol) for mol in dataframe[mol_col].values]
        
        df_descriptor=pd.DataFrame(desc,columns=names,index=dataframe[id_col].values)
        df_descriptor=df_descriptor.drop(columns='Ipc')

        idx, idy = np.where(pd.isnull(df_descriptor))

        not_procesed=list(np.unique(df_descriptor.index[idx]))
        df_descriptor.drop(index=not_procesed,inplace=True)

        df_descriptor.reset_index(inplace=True)

        df_descriptor.rename(columns={'index':id_col},inplace=True)

        rob=StandardScaler().fit(df_descriptor.iloc[:,1:])
        X_trans_rdk_sc=rob.transform(df_descriptor.iloc[:,1:])
        X_trans_rdk_sc=pd.DataFrame(np.c_[df_descriptor.iloc[:,0],X_trans_rdk_sc],index=df_descriptor.index.values,columns=df_descriptor.columns)
        
        final_df = dataframe.merge(X_trans_rdk_sc, how='left', on=id_col)

        return final_df, X_trans_rdk_sc
    
    def add_descriptors_and_fingerprints(self):
        """
            Adds RDKit descriptors and FPs into the dataframe

            :return comparison_df:
        """

        # Adds canonical SMILES
        self.compound_dataframe.loc[:, 'canon_smiles'] = self.compound_dataframe.loc[:,self.smiles_column].apply(lambda x: self.get_canonical_smiles(x))

        # Adds mol object from canonical SMILES
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'canon_smiles'].apply(lambda x: Chem.MolFromSmiles(x))

        # Adds Fingerprints
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: FingerprintMols.FingerprintMol(x))

        # Adds RDKit Descriptors
        self.compound_dataframe, self.descriptor_dataframe = self.get_rdkit_descriptors(self.compound_dataframe, self.molecule_id, 'mols_rdkit')