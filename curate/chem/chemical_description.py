"""
    Code that adds RDKit descriptors and Fingerprints into a pandas dataframe of chemical compounds

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 5/04/2024, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

from typing import Optional, Union

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

            :param smiles: SMILES to be processed

            :return canonical_smiles: canonical SMILES
        """

        try:
            canonical_smiles = Chem.CanonSmiles(smiles)
        except:
            canonical_smiles = None

        return canonical_smiles
    
    def get_morgan_fingerprints_ecfp6(self, dataframe: pd.DataFrame, id_col: str, nbits: int, radius: int, mol_col: str) -> pd.DataFrame:
        """
            Adds Morgan fingerprints in the dataframe with its proper name in each column

            :param dataframe: Pandas dataframe to work with
            :param nbits: number of bits of Morgan FPs.
            :param radius: radius to select the Morgan FPs.
            :param mol_col: column in the df containing the RDKit mol object

            :return final_df: intial dataframe with the FPs inlcuded
            :return df_morgan: df only including FPs and chemical identifier
        """

        ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=radius, nBits=nbits) for x in dataframe[mol_col]]

        ecfp6_name = [f'Bit_{i}' for i in range(nbits)]
        ecfp6_bits = [list(l) for l in ECFP6]
        df_morgan = pd.DataFrame(ecfp6_bits, index = dataframe[id_col], columns=ecfp6_name)
        df_morgan = df_morgan.reset_index()
        
        final_df = dataframe.merge(df_morgan, how='left', on=id_col)

        return final_df, df_morgan

    def select_specific_fingerprints(self, dataframe: pd.DataFrame, list_of_selected_fps: Union[list,str], number_of_top_features: int, sep: str = None) -> np.ndarray:
        """
            This function accepts a list of indexes of selected fingerprints so they can be retrieved from
            the whole set of bits and used to calculate similarities.
            Since the bits are an RDKit object they need to be converted into an array in order to perform the selection.
            
            :param dataframe: input dataframe to obtain and add the bit selection
            :param list_of_selected_fps: A list of the selected bits.
            :param top_features: integer indicating the number of features to select
            :param sep: separator to take into account from the input file

            :return final_df: Copy from the original dataframe inlcuding the selected fingerprints as a numpy array in one column
        """

        data_copy = dataframe.copy()

        if isinstance(list_of_selected_fps, list):
            """
                TODO:  need to process a list as argument
            """
            pass
        elif isinstance(list_of_selected_fps, str):
            features = pd.read_csv(list_of_selected_fps, sep=sep)
            top_features = features.sort_values(by='value', ascending=False).head(number_of_top_features).index

        feature_column = 'top_{}_features'.format(number_of_top_features)

        train_df = pd.DataFrame(columns=['CAS',feature_column])

        for i, row in data_copy.iterrows():
            small_list = []
            cas = row[self.molecule_id]
            morgan = row['morgan_fps']
            for j, bit in enumerate(list(morgan)):
                if j in top_features:
                    small_list.append(bit)
            train_df = train_df.append({'CAS': cas, feature_column: np.array(small_list)}, ignore_index=True)
        
        final_df = data_copy.merge(train_df, on='CAS')

        return final_df

    def get_rdkit_descriptors(self, dataframe: pd.DataFrame, id_col: str, mol_col: str) -> pd.DataFrame:
        """
            Gets a dataframe with RDKit descriptors

            :param dataframe: Pandas dataframe input
            :param id_col: column name in the dataframe that has de chemical identifier
            :param mol_col: column name in the dataframe that contains the RDKit mol object

            :return final_df: initial dataframe with descriptors included
            :return df_descriptor_complete: dataframe only containing rdkit descriptors and chemical identifier
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
        df_descriptor_complete=rob.transform(df_descriptor.iloc[:,1:])
        df_descriptor_complete=pd.DataFrame(np.c_[df_descriptor.iloc[:,0],df_descriptor_complete],index=df_descriptor.index.values,columns=df_descriptor.columns)
        
        final_df = dataframe.merge(df_descriptor_complete, how='left', on=id_col)

        return final_df, df_descriptor_complete
    
    def add_descriptors_and_fingerprints(self, morgan_bits: int = None, morgan_radius: int = None):
        """
            Adds RDKit descriptors and FPs into the dataframe

            :param morgan_bits: number of bits of Morgan FPs. If None, 2048 bits are selected.
            :param morgan_radius: radius to select the Morgan FPs. If None, 3 radius is selected

            :return compound_dataframe: dataframe containing all the information plus RDKit descriptors and Morgan FPs
            :return descriptor_dataframe: dataframe containing only the RDKit descriptors plus the chemical identifier
            :return morgan_dataframe: dataframe containing only the Morgan FPs plus the chemical identifier
        """

        # Checks if morgan_bits and morgan_radius are None
        if morgan_bits is None:
            morgan_bits = 2048
        
        if morgan_radius is None:
            morgan_radius = 3

        # Adds canonical SMILES
        self.compound_dataframe.loc[:, 'canon_smiles'] = self.compound_dataframe.loc[:,self.smiles_column].apply(lambda x: self.get_canonical_smiles(x))

        # Adds mol object from canonical SMILES
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'canon_smiles'].apply(lambda x: Chem.MolFromSmiles(x))

        # Adds RDKit Fingerprints for direct bulk similarity
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: FingerprintMols.FingerprintMol(x))

        # Adds Morgan Fingerprints for direct bulk similarity
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(), 'morgan_fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x,radius=morgan_radius, nBits=morgan_bits))
        
        # Adds RDKit Descriptors for modelling
        self.compound_dataframe, self.descriptor_dataframe = self.get_rdkit_descriptors(self.compound_dataframe, self.molecule_id, 'mols_rdkit')

        # Adds Morgan Fingerprints for modelling
        self.compound_dataframe, self.morgan_dataframe = self.get_morgan_fingerprints_ecfp6(self.compound_dataframe, id_col= self.molecule_id, nbits=morgan_bits, radius=morgan_radius, mol_col='mols_rdkit')