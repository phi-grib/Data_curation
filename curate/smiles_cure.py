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

        self.type_subs = {1:'normal', 2:'metal_ion', 3:'No non-salt/solvate components', 4:'Multiple non-salt/solvate components'}
        self.threads_ = {1:'non salt', 2:'salt', 3:'metal_ion', 4:'inorganic', 5:'peptide'}

    def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        self.smiles = smiles
        smiles_mol = Chem.MolFromSmiles(self.smiles)
    
        return smiles_mol

    def get_type_of_sub(self, std_mol: rdkit.Chem.rdchem.Mol) -> Optional[str]:
        """
            Gets the type of substance obtained by process smiles

            :param std_mol:

            :return type_of_sub:
        """

        if isinstance(std_mol,dict):
            for values in std_mol.values():
                type_of_sub = values[-1]
        else:
            type_of_sub = 'normal'
    
        return type_of_sub

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
            std_mol = ('failed_standardization', None)

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
        type_of_sub = self.get_type_of_sub(std_mol)

        if isinstance(std_mol, tuple):
            smiles_curated, subs_type = std_mol

        elif not std_mol:
            smi_ = Chem.MolToSmiles(smiles_mol)
            for metal in ps._metals:
                if metal in smi_:
                    smiles_curated = smi_
                    subs_type = 2
                    break
            else:
                smiles_curated = smi_
                subs_type = self.get_sub_type_id(type_of_sub)

        elif len(std_mol) > 1:
            new_mol = ps.std(smiles_mol)
            smiles = list(new_mol.keys())[0]
            for key in std_mol.keys():
                if key in ps._metals:
                    subs_type = 2
                    smiles_curated = self.smiles
                    break
                else:
                    subs_type = self.get_sub_type_id(type_of_sub)
                    smiles_curated = self.smiles
                    break

        else:
            for smiles in std_mol:
                for metal in ps._metals:
                    if metal in smiles:
                        smiles_curated = smiles
                        subs_type = 2
                        continue
                else:
                    smiles_curated = smiles
                    subs_type = self.get_sub_type_id(type_of_sub)

        return smiles_curated, subs_type

    def return_curate_smiles(self, smiles: str) -> Union[str,Tuple[str,bool]]:
        """
            Handles all the curation process.

            :return smiles, substance_type_id
        """

        smi_mol = self.get_rdkit_mol(smiles)

        if smi_mol is None:
            smiles = 'failed_structure_rdkit'
            subs_type = None
        else:
            smiles, subs_type = self.process_smiles_mol(smi_mol)

        return smiles, subs_type

    #### Substance type filters
    
    def check_molecular_weight(self, molecule: str) -> bool:
        """
            Checks MW of the molecule.
            Range betwn 75 kDa to 150 kDa.
            More than 150 kDa is considered too big and removed.

            :param molecule:

            :return bool:
        """

        pass
    
    def check_organometallic(self, molecule: str) -> Optional[str]:
        """
            Checks if the molecule has metal ions.

            :param molecule:

            :return metal_molecule:
        """
        pass

    def check_peptide(self, molecule: str) -> Optional[str]:
        """
            Checks if molecule is a peptide.

            :param molecule:

            :return peptide_molecule:
        """

        pass