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

from rdkit import Chem
from typing import Optional, Union

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps


class Curator(object):
    """
        Initializes the class with a SMILES input and applies standardization and other functions to curate the data
    """

    def __init__(self, smiles: str):
        """
            Initialized class with SMILES string

            :param smiles_dataframe: Dataframe containing the input SMILES
            :param smiles_field: Column name in the DF containing SMILES
        """

        self.smiles_str = smiles
        self.type_subs = {1:'normal',2:'metal_ion',3:'No non-salt/solvate components',4:'Multiple non-salt/solvate components'}
        
    def get_rdkit_mol(self) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        smiles_mol = Chem.MolFromSmiles(self.smiles_str)
    
        return smiles_mol

    def return_curate_smiles(self) -> Union[str,Tuple[str,bool]]:
        """
            Standardises SMILES. Handles exceptions

            :return smiles, substance_type_id
        """

        smi_mol = self.get_rdkit_mol()

        if smi_mol is None:
            smiles = 'failed_structure_rdkit'
            subs_type = None
        
        else:
            try:
                std_mol = ps.std(smi_mol, returnMetals=True)
                type_of_sub = ''
                for values in std_mol.values():
                    type_of_sub = values[-1]
            except:
                smiles = 'failed_standardization'
                subs_type = None
            
            if len(std_mol) > 1:
                new_mol = ps.std(mol)
                smiles = list(new_mol.keys())[0]
                for key in std_mol.keys():
                    if key in ps._metals:
                        subs_type = 2
                        smiles = struc
                        break
                    else:
                        for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                        substance_structures.loc[substance_structures.index == idx, 'substance_type_id'] = id_to_add
                        substance_structures.loc[substance_structures.index == idx, 'structure_curated'] = struc
                        break

            elif not std_mol:
                smi_ = Chem.MolToSmiles(mol)
                for metal in ps._metals:
                    if metal in smi_:
                        substance_structures.loc[substance_structures.index == idx, 'structure_curated'] = smi_
                        substance_structures.loc[substance_structures.index == idx, 'substance_type_id'] = 2
                        break
                else:
                    for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                    substance_structures.loc[substance_structures.index == idx, 'structure_curated'] = smi_
                    substance_structures.loc[substance_structures.index == idx, 'substance_type_id'] = id_to_add

            else:
                for smiles in std_mol:
                    for metal in ps._metals:
                        if metal in smiles:
                            substance_structures.loc[substance_structures.index == idx, 'structure_curated'] = smiles
                            substance_structures.loc[substance_structures.index == idx, 'substance_type_id'] = 2
                            continue
                    else:
                        for id_, type_ in type_subs.items():
                            if type_of_sub == type_:
                                id_to_add = id_
                            else:
                                id_to_add = 1
                        substance_structures.loc[substance_structures.index == idx, 'structure_curated'] = smiles
                        substance_structures.loc[substance_structures.index == idx, 'substance_type_id'] = id_to_add

        return smiles, subs_type
