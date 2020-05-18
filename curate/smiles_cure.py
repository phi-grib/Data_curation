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

        self.type_subs = {1:'normal',2:'metal_ion',3:'No non-salt/solvate components',4:'Multiple non-salt/solvate components'}

    def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        smiles_mol = Chem.MolFromSmiles(smiles)
    
        return smiles_mol

    def standardise_mol(self, smiles_mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
        """
            Standardises SMILES.

            :param smiles_mol:
            
            :return std_mol:
        """

        try:
            std_mol = ps.std(smiles_mol, returnMetals=True)
            for values in std_mol.values():
                type_of_sub = values[-1]
            
            if len(std_mol) > 1:
                new_mol = ps.std(mol)
                smiles = list(new_mol.keys())[0]
                for key in std_mol.keys():
                    if key in ps._metals:
                        subs_type = 2
                        smiles_curated = struc
                        break
                    else:
                        for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                        subs_type = id_to_add
                        smiles_curated = struc
                        break

            elif not std_mol:
                smi_ = Chem.MolToSmiles(mol)
                for metal in ps._metals:
                    if metal in smi_:
                        smiles_curated = smi_
                        subs_type = 2
                        break
                else:
                    for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                    smiles_curated = smi_
                    subs_type = id_to_add

            else:
                for smiles in std_mol:
                    for metal in ps._metals:
                        if metal in smiles:
                            smiles_curated = smiles
                            subs_type = 2
                            continue
                    else:
                        for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                        smiles_curated = smiles
                        subs_type = id_to_add
        except:
            smiles_curated = 'failed_standardization'
            subs_type = None

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
            smiles, subs_type = self.standardise_mol(smi_mol)

        return smiles, subs_type
