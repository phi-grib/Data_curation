"""
    Code for curating SMILES. To be used from a list/pandas dataframe.
    If SMILES are in a text file, first it will be processed either as a python list or a pandas dataframe.

    In principle it is written to work with CII and CR databases, but eventually it should be extended for 
    all phi projects if needed.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 27/05/2020, 12:16 PM
"""

import numpy as np
import pandas as pd
import re
import rdkit
import sys

from chembl_structure_pipeline import standardizer
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover
from typing import Optional, Union, Tuple

from curate.chem import process_smiles as ps
# from phitools import moleculeHelper as mh

class Curator(object):
    """
        Initializes the class with a SMILES input and applies 
        standardization and other functions to curate the data
    """

    def __init__(self):
        """
            Emtpy class initialization
        """

    def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
        """
            Returns mol object from rdkit

            :return smiles_mol:
        """

        self.smiles = smiles
        self.smiles_mol = self.smiles_to_rdkit_mol(smiles)

    def smiles_to_rdkit_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """
            Converts a SMILES string to a RDKit molecule.

            :param smiles: SMILES string of the molecule

            :returns mol: RDKit Mol, None if the SMILES string is invalid
        """

        self.no_san = False
        try:
            mol = Chem.MolFromSmiles(smiles)
        except TypeError:
            sys.stderr.write('Please check your structures and remove the NaNs\n')
            raise
        
        #  Sanitization check (detects invalid valence)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            self.no_san = True

        return mol 
        
    def filter_smiles(self) -> str:
        """
            Filters SMILES by checking different aspects:
                - organic
                - organometallic
                - peptide
                - inorganic
                - inorganic metal
                - organic salt
                - inorganic salt
            
            If SMILES passes these filters (meaning it's NOT any of those above)
            curation process continues.
            If not, we keep that SMILES as the curated one and get the 
            substance type to store it in the database.

            :return substance_type, final_smi:
        """

        if self.no_san:
            sub_type = 'no_sanitizable'
            sub_type_ns = self.check_organic_inorganic(self.smiles, self.smiles_mol, no_sanitizable=True)
            checker = self.check_organometallic(self.smiles)
            if checker:
                checker = '_'.join([sub_type,checker])
            else:
                sub_type = '_'.join([sub_type,sub_type_ns])
        else:
            sub_type = self.check_organic_inorganic(self.smiles, self.smiles_mol)
        
        checker = None

        if sub_type == 'organic':
            checker = self.check_organometallic(self.smiles)
            if not checker:
                checker = self.check_peptide(self.smiles)
            if not checker:
                checker = self.check_salt(self.smiles, sub_type)
        elif sub_type == 'inorganic':
            checker = self.check_inorganic_metal(self.smiles)
            if not checker:
                checker = self.check_salt(self.smiles, sub_type)
        
        if not checker:
            substance_type = sub_type
        else:
            substance_type = checker
        
        final_smi = self.check_errors(self.smiles_mol)

        return substance_type, final_smi

    def check_errors(self, smi: str) -> str:
        """
            This function processes the SMILES in order to canonicalize it and detect any errors.
            If errors are detected, the returned SMILES is not sanitized.

            :param smi: SMILES string of the compound

            :return final_smi: canonicalized SMILES
        """

        try:
            final_smi = standardizer.standardize_mol(smi)
        except:
            final_smi = smi
        final_smi = Chem.MolToSmiles(final_smi)
        final_smi = self.salt_remover(final_smi)

        return final_smi

    #### Checkers

    def check_organic_inorganic(self, molecule: str, mol_object: Chem.Mol, no_sanitizable=False) -> str:
        """
            Checks if there's a carbon atom in the molecule by filtering which elements that include a C or a c are not carbon 
            but others such as Ca (calcium), Cu (copper) or Cs (cessium). Also, a compound is considered to be organic when it has 
            a Carbon atom bound to a Hydrogen. This means that grafite, CO2, CO, NaCN etc... are considered inorganic.
            https://www.britannica.com/science/inorganic-compound

            :param molecule:

            :return substance_type:
        """

        C_upper = re.compile(r'C[^saeroudnfl]')
        c_lower = re.compile(r'[^SATM]c')

        if re.search(C_upper, molecule) or re.search(c_lower, molecule):
            if no_sanitizable:
                substance_type = 'organic'
            else:
                h_ = self.hydrogen_check(mol_object)
                if h_ == 'organic':
                    substance_type = h_
                else:
                    hal_check = self.halogen_check(molecule)
                    substance_type = hal_check
        else:
            substance_type = 'inorganic'
        
        return substance_type

    def hydrogen_check(self, molecule_object: Chem.Mol) -> Union[bool,str]:
        """
            This function checks the presence of Hydrogen atoms in a molecule with
            Carbons.

            :param molecule_object:

            :return h_check:
        """

        # Hydrogen check
        hs_pattern_1 = re.compile(r'\(\[H\]\).?[Cc]')
        hs_pattern_2 = re.compile(r'\[H\].?[Cc]')
        hs_pattern_3 = re.compile(r'[Cc]\(\[H\]\)')

        mol_hs = Chem.AddHs(molecule_object)
        smi_hs = Chem.MolToSmiles(mol_hs)
        if re.search(hs_pattern_1, smi_hs) or re.search(hs_pattern_2, smi_hs) or re.search(hs_pattern_3, smi_hs):
            # This pattern here will always be true since C or c will be present even if there is no H surrounding them. 
            # Need to check that at least ONE H is near the C
            h_check = 'organic'
        else:
            h_check = False
        
        return h_check

    def halogen_check(self, molecule_string: str) -> str:
        """
            This function checks the presence of halogen elements in the substance
            once the Carbon atom has been identified.

            :param molecule_string: SMILES string of the substance

            :return hal_check
        """

        # halogen check
        hal_check = 'inorganic'
        halogen_elements = ['Cl','Br']
        i_pattern = re.compile(r'I[^rn]')
        f_pattern = re.compile(r'F[^rem]')
        
        if re.search(i_pattern, molecule_string) or re.search(f_pattern, molecule_string):
            hal_check = 'organic'
        
        if hal_check == 'inorganic':
            for element in halogen_elements:
                if element in molecule_string:
                    hal_check = 'organic'

        return hal_check

    def check_organometallic(self, molecule: str) -> Optional[str]:
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
            for metal in ps._metals:
                if metal in self.smiles:
                    metals = metal

        if metals:
            metal_molecule = 'organometallic'

        return metal_molecule

    def check_peptide(self, molecule: str) -> Optional[str]:
        """
            Checks if molecule is a peptide

            :param molecule:

            :return peptide_mol:
        """

        peptide_mol = None

        l_peptide_pattern = re.compile(r'N\[C@@H\].+C\(=O\)O')
        r_peptide_pattern = re.compile(r'N\[C@H\].+C\(=O\)O')
        
        if re.search(l_peptide_pattern, molecule) or re.search(r_peptide_pattern, molecule):
            peptide_mol = 'peptide'
        
        return peptide_mol

    def check_inorganic_metal(self, molecule: str) -> Optional[str]:
        """
            Checks if molecule is only an elemental metal.

            :param molecule:

            :return metal:
        """

        metal = None
        
        if molecule in ps._metals:
            metal = 'inorganic_metal'
        
        return metal

    def check_salt(self, molecule: str, subType: str) -> str:
        """
            Checks if the molecule is salt.

            :param molecule:

            :return salt:
        """

        remover = SaltRemover()
        salt = None
        
        res, deleted = remover.StripMolWithDeleted(self.smiles_mol)
        
        if len(deleted) >= 1:
            salt = '_'.join([subType,'salt'])

        return salt
    
    def salt_remover(self, smiles: str) -> str:
        """
            Removes salts and counterions. Non sanitizable molecules can't be processed

            :param smiles: smiles string

            :return cleaned_smiles: 
        """
        
        rmv = rdMolStandardize.LargestFragmentChooser()
        
        if "." in smiles and Chem.MolFromSmiles(smiles):
            cleaned_smiles = Chem.MolToSmiles(rmv.choose(self.smiles_mol))
        else:
            cleaned_smiles = smiles

        return cleaned_smiles