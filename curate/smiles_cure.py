"""
    Code for curating SMILES. To be used from a list/pandas dataframe.
    If SMILES are in a text file, first it will be processed either as a python list or a pandas dataframe.

    In principle it is written to work with CII and CR databases, but eventually it should be extended for 
    all phi projects if needed.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 15/05/2020, 16:29 AM
"""

import numpy as np
import pandas as pd
import rdkit

from rdkit import Chem
from typing import Optional, Union, Tuple

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps


class Curator(object):
    """
        Initializes the class with a SMILES input and applies standardization and other functions to curate the data
    """

    def __init__(self):
        """
            Initialized class with SMILES string

            :param smiles_dataframe: Dataframe containing the input SMILES
            :param smiles_field: Column name in the DF containing SMILES
        """

        self.threads_ = {1:'non salt', 2:'salt', 3:'metal_ion', 4:'inorganic', 5:'peptide'}

    def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        self.smiles = smiles
        smiles_mol = Chem.MolFromSmiles(self.smiles)
    
        return smiles_mol

    def get_sub_type_id(self, type_of_sub: str) -> int:
        """
            Gets the substance type id according to the dictionary. If there's no coincidence, returns 1 (normal)

            :param type_:

            :return sub_type_id:
        """

        for id_, type_ in self.type_subs.items():
            if type_of_sub == type_:
                sub_type_id = id_
            else:
                sub_type_id = 1
        
        return sub_type_id

    def standardise_mol(self, smiles_mol: rdkit.Chem.rdchem.Mol) -> Union[dict,tuple]:
        """
            Standardises SMILES. If not, returns failed.

            :param smiles_mol:

            :return std_mol:
        """

        try:
            std_mol = ps.std(smiles_mol, returnMetals=True)
        except:
            std_mol = None

        return std_mol

    def process_smiles_mol(self, smiles_mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
        """
            Standardises SMILES.

            :param smiles_mol:
            
            :return std_mol:

            TODO: check loops. Also, change this nested code for different function blocks that identify the kind of substance 
            is given as input and also properly process the SMILES.
        """

        std_mol = self.standardise_mol(smiles_mol)
    
        if not std_mol:
            # smi_ = Chem.MolToSmiles(smiles_mol)
            # for metal in ps._metals:
            #     if metal in smi_:
            #         smiles_curated = smi_
            #         break
            # else:
            #     smiles_curated = smi_
            smiles_curated = None

        elif len(std_mol) > 1:
            print(Chem.MolToSmiles(smiles_mol), std_mol)
            for key in std_mol.keys():
                checker = self.check_organometallic(key)

        else:
            for smiles in std_mol:
                for metal in ps._metals:
                    if metal in smiles:
                        smiles_curated = smiles
                        continue
                else:
                    smiles_curated = smiles

        return smiles_curated

    def return_curate_smiles(self, smiles: str) -> Union[str,Tuple[str,bool]]:
        """
            Handles all the curation process.

            :return smiles, substance_type_id
        """

        smi_mol = self.get_rdkit_mol(smiles)

        if smi_mol is None:
            smiles = 'failed_structure_rdkit'
        else:
            smiles = self.process_smiles_mol(smi_mol)

        return smiles

    #### Substance type filters
    
    # def check_molecular_weight(self, molecule: rdkit.Chem.rdchem.Mol) -> bool:
    #     """
    #         Checks MW of the molecule.
    #         Range betwn 75 kDa to 150 kDa.
    #         More than 150 kDa is considered too big and removed.

    #         :param molecule:

    #         :return bool:
    #     """

    #     mol_weight = Chem.Descriptors.MolWt(molecule)
        
    #     if mol_weight > 150:
    #         bool_ = False
    #     else:
    #         bool_ = True
        
    #     return bool_
    
    def check_organometallic(self, molecule: str) -> Optional[str]:
        """
            Checks if the molecule has metal ions.

            :param molecule:

            :return metal_molecule:
        """

        for metal in ps.metals_:
            if metal in molecule:
                metal_molecule = True
        
        if metal_molecule:
            bool_ = True
        else:
            bool_ = False

        return bool_

    def check_peptide(self, molecule: str) -> Optional[str]:
        """
            Checks if molecule is a peptide.

            :param molecule:

            :return peptide_molecule:
        """

        pass
    
    def check_salt(self, molecule: str) -> Optional[str]:
        """
            Check if molecule is a salt

            :param molecule:

            :return salt_molecule:
        """

        pass