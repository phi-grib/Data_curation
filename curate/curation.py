"""
    Code for curating SMILES. To be used from a list/pandas dataframe.
    If SMILES are in a text file, first it will be processed either as a python list or a pandas dataframe.

    In principle it is written to work with CII and CR databases, but eventually it should be extended for 
    all phi projects if needed.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 27/05/2020, 12:16 AM
"""

import numpy as np
import pandas as pd
import re
import rdkit

from rdkit import Chem
from typing import Optional, Union, Tuple

from phitools import moleculeHelper as mh
from standardiser import process_smiles as ps

class Curator(object):
    """
        Initializes the class with a SMILES input and applies 
        standardization and other functions to curate the data
    """

    def __init__(self):
        """
            Initialized class with SMILES string

            :param smiles_dataframe: Dataframe containing the input SMILES
            :param smiles_field: Column name in the DF containing SMILES
        """

        self.threads_ = {1:'organic', 2:'salt', 3:'organometallic', 4:'inorganic',
        5:'peptide',6:'isomeric_mixture',7:'related_to_mixtures'}

    def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        self.smiles = smiles
        self.smiles_mol = Chem.MolFromSmiles(self.smiles)

    def filter_smiles(self) -> Union[rdkit.Chem.rdchem.Mol, str]:
        """
            Filters SMILES by checking different aspects:
                - organometallic
                - inorganic
                - isomeric mixture
                - related to mixtures
            
            If SMILES passes these filters (meaning it's NOT any of those above)
            curation process continues.
            If not, we keep that SMILES as the curated one and get the 
            substance type to store it in the database.

            :return filtered_smiles:
        """

        sub_type = self.check_organic_inorganic(self.smiles)
        if sub_type == 'organic':
            checker = self.check_organometallic(self.smiles)
        
        elif sub_type == 'inorganic':
            checker = self.check_isomeric_mixture(self.smiles)
            if not checker:
                checker = self.check_related_mixture(self.smiles)
        


    #### Checkers

    def check_organic_inorganic(self, molecule: str) -> str:
        """
            Checks if there's a carbon atom in the molecule by filtering which elements that include a C or a c are not carbon 
            but others such as Ca (calcium), Cu (copper) or Cs (cessium).

            :param molecule:

            :return substance_type:
        """

        C_upper = re.compile(r'C[^saeroudnfl]')
        c_lower = re.compile(r'[^SATM]c')

        if re.search(C_upper, molecule) or re.search(c_lower, molecule):
            substance_type = 'organic'
        else:
            substance_type = 'inorganic'

        return substance_type

    def check_organometallic(self, molecule: str) -> bool:
        """
            Checks if the molecule has metal ions.

            :param molecule:

            :return metal_molecule:
        """
        
        metal_molecule = False
        metals = None
        try:
            mol, metals = ps.disconnect(self.smiles_mol)
        except:
            print('errooooooooor',self.smiles, self.smiles_mol)
        if metals:
            metal_molecule = 'organometallic'
            print('organometallic',self.smiles, Chem.MolToSmiles(mol), metals)

        return metal_molecule

    def check_isomeric_mixture(self, molecule: str) -> bool:
        """
            Checks if the molecule is an isomeric mixture.

            :param molecule:

            :return isomeric_mixture_molecule:
        """

        isomeric_mixture_molecule = False

        pass

    def check_related_mixture(self, molecule: str) -> bool:
        """
            Checks if the molecule is related to mixtures.

            :param molecule:

            :return related_to_mixture_molecule:
        """

        related_to_mixture_molecule = False

        pass