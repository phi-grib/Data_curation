"""
    Code that adds RDKit descriptors and Fingerprints into a pandas dataframe of chemical compounds

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 5/04/2024, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

from typing import Optional, Union

from curate.util import get_logger

try:
    from rdkit.Chem.rdFingerprintGenerator import MorganGenerator
    _use_morgan_generator = True
except ImportError:
    _use_morgan_generator = False
    print("INFO: Could not import MorganGenerator. Falling back to older RDKit fingerprint method. "
          "Consider updating RDKit for better performance.")

LOG = get_logger(__name__)

class Description(object):
    """
        Calculates RDKit descriptors and Fingerprints for all compounds in dataframe
    """

    def __init__(self, dataframe: pd.DataFrame, molecule_id: str, smiles_column: str):
        
        self.compound_dataframe = dataframe.copy()
        self.molecule_id = molecule_id
        self.smiles_column = smiles_column
        self.descriptor_dataframe = None
        self.morgan_dataframe = None
        self.nan_mols = None

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
    
    # def get_morgan_fingerprints_ecfp6(dataframe: pd.DataFrame, id_col: str, nbits: int, radius: int, mol_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    #     """
    #     Adds Morgan fingerprints to the dataframe using the modern MorganGenerator
    #     with a fallback to the older method for backward compatibility.

    #     :param dataframe: Pandas dataframe to work with.
    #     :param id_col: Column name for the chemical identifier.
    #     :param nbits: Number of bits for the Morgan Fingerprints.
    #     :param radius: The radius for the Morgan Fingerprints (e.g., radius=2 for ECFP4, radius=3 for ECFP6).
    #     :param mol_col: Column in the df containing the RDKit mol object.

    #     :return final_df: Initial dataframe with the fingerprints included.
    #     :return df_morgan: DataFrame containing only the fingerprints and chemical identifier.
    #     """
    #     ECFP6 = []
    #     if _use_morgan_generator:
    #         # --- Modern, Efficient Method ---
    #         fp_generator = MorganGenerator(radius=radius, nBits=nbits)
    #         ECFP6 = [fp_generator.GetFingerprintAsBitVect(x) for x in dataframe[mol_col]]
    #     else:
    #         # --- Fallback for Older RDKit Versions ---
    #         ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x, radius=radius, nBits=nbits) for x in dataframe[mol_col]]

    #     # The rest of the function remains the same
    #     ecfp6_name = [f'Bit_{i}' for i in range(nbits)]
    #     ecfp6_bits = [list(l) for l in ECFP6]
    #     df_morgan = pd.DataFrame(ecfp6_bits, index=dataframe[id_col], columns=ecfp6_name)
    #     df_morgan = df_morgan.reset_index()
        
    #     final_df = dataframe.merge(df_morgan, how='left', on=id_col)

    #     return final_df, df_morgan
    
    def get_morgan_fingerprints_ecfp6(self, dataframe: pd.DataFrame, id_col: str, nbits: int, radius: int, mol_col: str):
        """
        Calculates Morgan fingerprints using a highly memory-efficient
        np.frombuffer approach.
        """
        mol_list = dataframe[mol_col].tolist()
        
        fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=radius, nBits=nbits) for m in mol_list]

        # --- Highly Memory-Efficient Conversion using np.frombuffer ---
        # 1. Pre-allocate the final NumPy array.
        fp_array = np.zeros((len(fps), nbits), dtype=np.int8)

        # 2. Loop through each fingerprint, convert to a byte string, and
        #    use frombuffer to create a row without intermediate Python objects.
        for i, fp in enumerate(fps):
            byte_string = fp.ToBitString().encode('utf-8')
            # frombuffer reads the raw bytes, then we convert ASCII '0'/'1' to int 0/1
            fp_array[i] = np.frombuffer(byte_string, dtype=np.uint8) - ord('0')

        ecfp6_names = [f'Bit_{i}' for i in range(nbits)]
        df_morgan = pd.DataFrame(fp_array, index=dataframe[id_col].values, columns=ecfp6_names)
        
        df_morgan.reset_index(inplace=True)
        df_morgan.rename(columns={'index': id_col}, inplace=True)

        self.morgan_dataframe = df_morgan

    def select_specific_fingerprints(self, list_of_selected_fps: Union[list,str], number_of_top_features: int, sep: str = None) -> np.ndarray:
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

        data_copy = self.compound_dataframe.copy()

        if isinstance(list_of_selected_fps, list):
            """
                TODO:  need to process a list as argument
            """
            pass
        elif isinstance(list_of_selected_fps, str):
            features = pd.read_csv(list_of_selected_fps, sep=sep)
            top_features = features.sort_values(by='value', ascending=False).head(number_of_top_features).index

        feature_column = 'top_{}_features'.format(number_of_top_features)

        train_df = pd.DataFrame(columns=[self.molecule_id,feature_column])

        for i, row in data_copy.iterrows():
            small_list = []
            cas = row[self.molecule_id]
            morgan = row['morgan_fps']
            for j, bit in enumerate(list(morgan)):
                if j in top_features:
                    small_list.append(bit)
            train_df = train_df.append({self.molecule_id: cas, feature_column: np.array(small_list)}, ignore_index=True)
        
        final_df = data_copy.merge(train_df, on=self.molecule_id)

        self.compound_dataframe = final_df

    # def get_rdkit_descriptors(self, dataframe: pd.DataFrame, id_col: str, mol_col: str) -> pd.DataFrame:
    #     """
    #         Gets a dataframe with RDKit descriptors

    #         :param dataframe: Pandas dataframe input
    #         :param id_col: column name in the dataframe that has de chemical identifier
    #         :param mol_col: column name in the dataframe that contains the RDKit mol object

    #         :return final_df: initial dataframe with descriptors included
    #         :return df_descriptor_complete: dataframe only containing rdkit descriptors and chemical identifier
    #     """

    #     names=[x[0] for x in Chem.Descriptors._descList]
    #     calculator=MoleculeDescriptors.MolecularDescriptorCalculator(names)
                
    #     desc = [calculator.CalcDescriptors(mol) for mol in dataframe[mol_col].values]

    #     #### FIX for mixed data types in id col and future df_descriptor index
    #     dataframe[id_col] = dataframe[id_col].astype(str)
        
    #     df_descriptor=pd.DataFrame(desc,columns=names,index=dataframe[id_col].values)
    #     df_descriptor=df_descriptor.drop(columns='Ipc')

    #     idx, idy = np.where(pd.isnull(df_descriptor))

    #     not_procesed=list(np.unique(df_descriptor.index[idx]))
    #     df_descriptor.drop(index=not_procesed,inplace=True)

    #     df_descriptor.reset_index(inplace=True)

    #     df_descriptor.rename(columns={'index':id_col},inplace=True)

    #     ### Handling of problematic rows (infinite or very large numbers)

    #     problematic_rows_partial = self.find_problematic_rows(df_descriptor.iloc[:, 1:])

    #     if not problematic_rows_partial.empty:
    #         # Create a new DataFrame with the original rows
    #         self.problematic_rows_full = df_descriptor.loc[problematic_rows_partial.index].copy() # .copy() to avoid SettingWithCopyWarning

    #         # Remove problematic rows from df_descriptor
    #         df_descriptor = df_descriptor.drop(problematic_rows_partial.index)

    #     # Separate the ID column from the feature columns
    #     id_series = df_descriptor.iloc[:, 0]
    #     feature_df = df_descriptor.iloc[:, 1:]

    #     # Fit and transform the numeric features ONLY
    #     scaler = StandardScaler().fit(feature_df)
    #     scaled_features = scaler.transform(feature_df)

    #     # Create the final DataFrame with the scaled numeric data and correct columns/index
    #     df_descriptor_complete = pd.DataFrame(scaled_features, index=feature_df.index, columns=feature_df.columns)

    #     # Re-insert the ID column at the beginning, preserving its type
    #     df_descriptor_complete.insert(0, id_col, id_series)
        
    #     final_df = dataframe.merge(df_descriptor_complete, how='left', on=id_col)

    #     return final_df, df_descriptor_complete

    def get_rdkit_descriptors(self, dataframe: pd.DataFrame, id_col: str, mol_col: str):
        """
        Calculates a dataframe with RDKit descriptors.
        """
        mol_list = dataframe[mol_col].tolist()
        names = [x[0] for x in Chem.Descriptors._descList]
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(names)
        
        desc = [calculator.CalcDescriptors(mol) for mol in mol_list]

        df_descriptor = pd.DataFrame(desc, columns=names, index=dataframe[id_col].values)
        df_descriptor.drop(columns='Ipc', inplace=True) # Drop descriptor known to cause issues

        # Handle NaNs and Infs produced by RDKit calculation
        df_descriptor.replace([np.inf, -np.inf], np.nan, inplace=True)
        rows_with_nan = df_descriptor[df_descriptor.isnull().any(axis=1)].index
        if not rows_with_nan.empty:
            LOG.warning(f"Dropping {len(rows_with_nan)} rows due to NaN/Inf values in RDKit descriptors.")
            df_descriptor.drop(index=rows_with_nan, inplace=True)

        df_descriptor.reset_index(inplace=True)
        df_descriptor.rename(columns={'index': id_col}, inplace=True)
        
        self.descriptor_dataframe = df_descriptor

    def find_problematic_rows(self, df, threshold=1e15):
        """
        Finds rows containing infinity or values exceeding a threshold.

        Args:
            df (pd.DataFrame): The DataFrame to check.
            threshold (float): The threshold for large values.

        Returns:
            pd.DataFrame: A DataFrame containing the problematic rows.
        """

        inf_mask = np.isinf(df)
        large_value_mask = abs(df) > threshold
        combined_mask = inf_mask | large_value_mask
        problematic_rows = df[combined_mask.any(axis=1)]
        
        return problematic_rows

    # def add_descriptors_and_fingerprints(self, morgan_bits: int = None, morgan_radius: int = None):
    #     """
    #         Adds RDKit descriptors and FPs into the dataframe

    #         :param morgan_bits: number of bits of Morgan FPs. If None, 2048 bits are selected.
    #         :param morgan_radius: radius to select the Morgan FPs. If None, 3 radius is selected

    #         :return compound_dataframe: dataframe containing all the information plus RDKit descriptors and Morgan FPs
    #         :return descriptor_dataframe: dataframe containing only the RDKit descriptors plus the chemical identifier
    #         :return morgan_dataframe: dataframe containing only the Morgan FPs plus the chemical identifier
    #     """

    #     # Checks if morgan_bits and morgan_radius are None
    #     if morgan_bits is None:
    #         morgan_bits = 2048
        
    #     if morgan_radius is None:
    #         morgan_radius = 3

    #     # Adds canonical SMILES
    #     self.compound_dataframe.loc[:, 'canon_smiles'] = self.compound_dataframe.loc[:,self.smiles_column].apply(lambda x: self.get_canonical_smiles(x))

    #     # Adds mol object from canonical SMILES
    #     self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'canon_smiles'].apply(lambda x: Chem.MolFromSmiles(x))

    #     # Adds a check for mols_rdkit NaN values
    #     self.nan_mols = self.compound_dataframe.loc[self.compound_dataframe['mols_rdkit'].isna()]
    #     self.compound_dataframe = self.compound_dataframe.loc[~self.compound_dataframe['mols_rdkit'].isna()]
        
    #     # Adds RDKit Fingerprints for direct bulk similarity
    #     self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: FingerprintMols.FingerprintMol(x))

    #     # Adds Morgan Fingerprints for direct bulk similarity
    #     self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(), 'morgan_fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x,radius=morgan_radius, nBits=morgan_bits))
        
    #     # Adds RDKit Descriptors for modelling
    #     self.compound_dataframe, self.descriptor_dataframe = self.get_rdkit_descriptors(self.compound_dataframe, self.molecule_id, 'mols_rdkit')

    #     # Adds Morgan Fingerprints for modelling
    #     self.compound_dataframe, self.morgan_dataframe = self.get_morgan_fingerprints_ecfp6(self.compound_dataframe, id_col= self.molecule_id, nbits=morgan_bits, radius=morgan_radius, mol_col='mols_rdkit')
    
    def add_descriptors_and_fingerprints(self, morgan_bits: int = 2048, morgan_radius: int = 3):
        """
        Main orchestration method to calculate all descriptors and fingerprints.
        This workflow is optimized to minimize iterations over the dataframe.
        """
        LOG.info("Canonicalizing SMILES and generating RDKit Mol objects...")
        # --- Step 1: Canonicalize SMILES and create Mol objects in a single pass ---
        mols = []
        canon_smiles_list = []
        for smi in self.compound_dataframe[self.smiles_column]:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    canon_smiles_list.append(Chem.MolToSmiles(mol, canonical=True))
                    mols.append(mol)
                else:
                    canon_smiles_list.append(None)
                    mols.append(None)
            except Exception:
                canon_smiles_list.append(None)
                mols.append(None)

        self.compound_dataframe['canon_smiles'] = canon_smiles_list
        self.compound_dataframe['mols_rdkit'] = mols

        # --- Step 2: Handle molecules that failed to process ---
        self.nan_mols = self.compound_dataframe[self.compound_dataframe['mols_rdkit'].isna()]
        LOG.warning(f"Found {len(self.nan_mols)} molecules that could not be processed.")
        # Keep only valid molecules for further processing
        self.compound_dataframe.dropna(subset=['mols_rdkit'], inplace=True)
        
        # --- Step 3: Calculate Descriptors and Fingerprints ---
        LOG.info("Calculating RDKit descriptors...")
        self.get_rdkit_descriptors(self.compound_dataframe, self.molecule_id, 'mols_rdkit')
        LOG.info("Calculating Morgan fingerprints...")
        self.get_morgan_fingerprints_ecfp6(self.compound_dataframe, self.molecule_id, morgan_bits, morgan_radius, 'mols_rdkit')