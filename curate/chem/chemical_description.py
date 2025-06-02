"""
    Code that adds RDKit descriptors and Fingerprints into a pandas dataframe of chemical compounds

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 5/04/2024, 11:35 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

from typing import Optional, Union, List, Tuple

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

        ### Handling of problematic rows (infinite or very large numbers)

        problematic_rows_partial = self.find_problematic_rows(df_descriptor.iloc[:, 1:])

        if not problematic_rows_partial.empty:
            # Create a new DataFrame with the original rows
            self.problematic_rows_full = df_descriptor.loc[problematic_rows_partial.index].copy() # .copy() to avoid SettingWithCopyWarning

            # Remove problematic rows from df_descriptor
            df_descriptor = df_descriptor.drop(problematic_rows_partial.index)

        rob=StandardScaler().fit(df_descriptor.iloc[:,1:])
    
        df_descriptor_complete=rob.transform(df_descriptor.iloc[:,1:])
        df_descriptor_complete=pd.DataFrame(np.c_[df_descriptor.iloc[:,0],df_descriptor_complete],index=df_descriptor.index.values,columns=df_descriptor.columns)
        
        final_df = dataframe.merge(df_descriptor_complete, how='left', on=id_col)

        return final_df, df_descriptor_complete
    
    # def calculate_descriptors_with_error_handling(self, mols: List[Optional[Chem.Mol]], 
    # descriptor_names: List[str]) -> Tuple[List[List[Optional[float]]], List[Optional[Chem.Mol]]]:
    #  NOT WORKNG. IT DOESN'T CAPTURE THE EXCEPTION, WHICH IS THROWN BY RDKIT AND NOT CAPTURED BY THE PYTHON WRAPPER
    #     """
    #         Calculates molecular descriptors for a list of RDKit molecules, handling potential ZeroDivisionErrors 
    #         that can occur when calculating FpDensityMorgan descriptors for molecules with zero heavy atoms.

    #         Args:
    #             mols: A list of RDKit Mol objects. Molecules may be None or represent empty molecules.
    #             descriptor_names: A list of descriptor names (strings) to calculate.

    #         Returns:
    #             A tuple containing:
    #                 - all_results: A list of lists, where each inner list contains the calculated descriptor values 
    #                             for a molecule. If a ZeroDivisionError (or another error) occurred during calculation,
    #                             the corresponding inner list will contain None values.
    #                 - problematic_mols: A list of RDKit Mol objects that caused a ZeroDivisionError (or other issue)
    #                                     during descriptor calculation.

    #         Raises:
    #             (No exceptions are raised by this function itself. Exceptions during descriptor calculation are caught and handled)

    #         Example:
    #             >>> from rdkit import Chem
    #             >>> from rdkit.Chem import Descriptors
    #             >>> mols = [Chem.MolFromSmiles("C"), Chem.MolFromSmiles(""), Chem.MolFromSmiles("[H][H]")] #Example Molecules
    #             >>> descriptor_names = ["MolWt", "FpDensityMorgan1", "NumHAcceptors"]
    #             >>> results, problematic = calculate_descriptors_with_error_handling(mols, descriptor_names)
    #             >>> print(results)
    #             [[16.043, 1.0, 0], [None, None, None], [2.016, None, 0]]
    #             >>> print([Chem.MolToSmiles(x) if x else None for x in problematic])
    #             [None, '[H][H]']
    #     """
        
    #     descriptor_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
        
    #     all_results: List[List[Optional[float]]] = []
    #     problematic_mols: List[Optional[Chem.Mol]] = []

    #     for mol in mols:
    #         try:
    #             if mol is None:
    #                 raise ValueError("Mol object is None") #Explicitly handle None molecules.
    #             results = descriptor_calculator.CalcDescriptors(mol)
    #             all_results.append(results)
    #         except ZeroDivisionError:
    #             print('why not here?')
    #             problematic_mols.append(mol)
    #             print(f"Error: ZeroDivisionError encountered for molecule: {Chem.MolToSmiles(mol) if mol else 'Empty Molecule'}")

    #             # Find the descriptor that caused the problem.
    #             for desc_name in descriptor_names:
    #                 try:
    #                     desc_func = getattr(Descriptors, desc_name)
    #                     desc_func(mol) # Attempt calculation to see which one fails
    #                 except ZeroDivisionError:
    #                     print(f"   Problematic Descriptor: {desc_name}")
    #                     break # Break after finding the problem descriptor.
    #                 except Exception:
    #                     pass # Skip other exceptions during descriptor identification.
    #             all_results.append([None] * len(descriptor_names)) #Append None to keep data structure consistent.

    #         except (ValueError, AttributeError) as e:
    #             problematic_mols.append(mol)
    #             print(f"Error: {e} encountered for molecule: {Chem.MolToSmiles(mol) if mol else 'Empty Molecule'}")
    #             all_results.append([None] * len(descriptor_names))

    #         except Exception as e:
    #             problematic_mols.append(mol)
    #             print(f"Unexpected Error: {e} encountered for molecule: {Chem.MolToSmiles(mol) if mol else 'Empty Molecule'}")
    #             all_results.append([None] * len(descriptor_names))

    #     return all_results, problematic_mols

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

        # Adds a check for mols_rdkit NaN values
        self.nan_mols = self.compound_dataframe.loc[self.compound_dataframe['mols_rdkit'].isna()]
        self.compound_dataframe = self.compound_dataframe.loc[~self.compound_dataframe['mols_rdkit'].isna()]
        
        # Adds RDKit Fingerprints for direct bulk similarity
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: FingerprintMols.FingerprintMol(x))

        # Adds Morgan Fingerprints for direct bulk similarity
        self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(), 'morgan_fps'] = self.compound_dataframe.loc[~self.compound_dataframe['canon_smiles'].isna(),'mols_rdkit'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x,radius=morgan_radius, nBits=morgan_bits))
        
        # Adds RDKit Descriptors for modelling
        self.compound_dataframe, self.descriptor_dataframe = self.get_rdkit_descriptors(self.compound_dataframe, self.molecule_id, 'mols_rdkit')

        # Adds Morgan Fingerprints for modelling
        self.compound_dataframe, self.morgan_dataframe = self.get_morgan_fingerprints_ecfp6(self.compound_dataframe, id_col= self.molecule_id, nbits=morgan_bits, radius=morgan_radius, mol_col='mols_rdkit')