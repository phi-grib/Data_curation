"""
    Code for curating SMILES. To be used from a list/pandas dataframe.
    If SMILES are in a text file, first it will be processed either as a python list or a pandas dataframe.

    In principle it is written to work with CII and CR databases, but eventually it should be extended for 
    all phi projects if needed.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 27/05/2020, 12:16 PM
"""

import re
import rdkit

from chembl_structure_pipeline import standardizer
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover
from typing import Optional, Union

from curate.chem import process_smiles as ps
from curate.util import get_logger
LOG = get_logger(__name__)

class Curator(object):
    """
        Initializes the class with a SMILES input and applies 
        standardization and other functions to curate the data
    """

    def __init__(self):
        """
            Emtpy class initialization
        """

    # def get_rdkit_mol(self, smiles: str) -> rdkit.Chem.rdchem.Mol:
    #     """
    #         Returns mol object from rdkit

    #         :return smiles_mol:
    #     """

    #     self.smiles = smiles
    #     self.smiles_mol = self.smiles_to_rdkit_mol(smiles)
        
    def smiles_to_rdkit_mol(self, smiles: str):
        """
            Converts a SMILES string to a RDKit molecule and stores the SMILES and mol object into class variables.
            If the molecule is not sanitizable it raises the flag no_san for proper further processing.

            :param smiles: SMILES string of the molecule

        """

        self.smiles = smiles

        try:
            mol = Chem.MolFromSmiles(smiles)
            self.no_san = False
        except TypeError:
            LOG.error('Please check your structures and remove the NaNs\n')
            raise
        
        #  Sanitization check (detects invalid valence)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            self.no_san = True

        self.smiles_mol = mol
        
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
            sub_type_ns = self.check_organic_inorganic(self.smiles_mol, self.smiles, no_sanitizable=True)
            checker = self.check_organometallic()
            if checker:
                checker = '_'.join([sub_type,checker])
            else:
                sub_type = '_'.join([sub_type,sub_type_ns])
        else:
            sub_type = self.check_organic_inorganic(self.smiles_mol, self.smiles)
        
        checker = None

        if sub_type == 'organic':
            checker = self.check_organometallic()
            if not checker:
                checker = self.check_peptide(self.smiles)
            if not checker:
                checker = self.check_salt(sub_type)
        elif sub_type == 'inorganic':
            checker = self.check_inorganic_metal(self.smiles)
            if not checker:
                checker = self.check_salt(sub_type)
        
        if not checker:
            substance_type = sub_type
        else:
            substance_type = checker
        
        final_smi = self.check_errors(self.smiles_mol)

        return substance_type, final_smi

    def check_errors(self, mol_object: Chem.Mol) -> str:
        """
            This function processes the SMILES in order to canonicalize it and detect any errors.
            If errors are detected, the returned SMILES is not sanitized.

            :param mol_object: RDKit mol object of the compound

            :return final_smi: isomeric SMILES
        """

        if mol_object is None:
            final_smi = self.smiles
        else:
            try:
                final_mol_object = standardizer.standardize_mol(mol_object)
            except:
                final_mol_object = mol_object
            # final_smi = Chem.MolToSmiles(final_smi, isomericSmiles=True)
            final_smi = self.salt_remover(final_mol_object)

        return final_smi

    #### Checkers

    # def check_organic_inorganic(self, molecule: str, mol_object: Chem.Mol, no_sanitizable=False) -> str:
    #     """
    #         Checks if there's a carbon atom in the molecule by filtering which elements that include a C or a c are not carbon 
    #         but others such as Ca (calcium), Cu (copper) or Cs (cessium). Also, a compound is considered to be organic when it has 
    #         a Carbon atom bound to a Hydrogen. This means that grafite, CO2, CO, NaCN etc... are considered inorganic.
    #         https://www.britannica.com/science/inorganic-compound

    #         :param molecule:

    #         :return substance_type:
    #     """

    #     C_upper = re.compile(r'C[^saeroudnfl]')
    #     c_lower = re.compile(r'[^SATM]c')

    #     if re.search(C_upper, molecule) or re.search(c_lower, molecule):
    #         if no_sanitizable:
    #             substance_type = 'organic'
    #         else:
    #             h_ = self.hydrogen_check(mol_object)
    #             if h_ == 'organic':
    #                 substance_type = h_
    #             else:
    #                 hal_check = self.halogen_check(molecule)
    #                 substance_type = hal_check
    #     else:
    #         substance_type = 'inorganic'
        
    #     return substance_type

    def check_organic_inorganic(self, mol: Chem.Mol, smiles: str, no_sanitizable: bool = False) -> str:
        """
            Determines whether a given molecule is organic or inorganic based on its structure and composition.

            This function uses a combination of RDKit Mol object analysis and SMILES string pattern matching
            to classify molecules. It checks for the presence of carbon atoms and specific bond types
            characteristic of organic compounds.

            Parameters:
            -----------
            mol : Chem.Mol
                An RDKit Mol object representing the molecule. This should be a valid, sanitized
                molecular structure. If the molecule cannot be sanitized, set no_sanitizable to True
                and provide the SMILES string instead.

            smiles : str
                A string representing the molecule in SMILES (Simplified Molecular Input Line Entry System) format.
                This is used for pattern matching when the RDKit Mol object is not available or cannot be sanitized.

            no_sanitizable : bool, optional (default=False)
                A flag indicating whether the molecule can be sanitized by RDKit. If True, the function
                will rely solely on SMILES string analysis for classification.

            Returns:
            --------
            str
                Either 'organic' or 'inorganic', indicating the classification of the molecule.

            Notes:
            ------
            The function considers a molecule organic if it contains carbon atoms bonded to hydrogen,
            halogens, or other organic elements (N, P, S). Simple inorganic carbon compounds (e.g., CO2)
            are classified as inorganic.

            The classification process follows these steps:
            1. If no_sanitizable is True, it uses SMILES-based methods for classification.
            2. It checks if the molecule is a simple inorganic compound.
            3. If the molecule contains carbon, it checks for specific organic bonds.
            4. If none of the above conditions are met, it classifies the molecule as inorganic.

            Examples:
            ---------
            >>> mol = Chem.MolFromSmiles("CCO")
            >>> smiles = "CCO"
            >>> classifier.check_organic_inorganic(mol, smiles)
            'organic'

            >>> mol = Chem.MolFromSmiles("NaCl")
            >>> smiles = "NaCl"
            >>> classifier.check_organic_inorganic(mol, smiles)
            'inorganic'

            >>> smiles = "C1=CC=CC=C1"
            >>> classifier.check_organic_inorganic(None, smiles, no_sanitizable=True)
            'organic'
        """

        if no_sanitizable:
            if self.has_carbon_smiles(smiles):
                if self.has_carbon_hydrogen_bond_smiles(smiles) or self.has_carbon_halogen_bond_smiles(smiles) or self.has_carbon_organic_element_bond_smiles(smiles):
                    return 'organic'
        
        elif self.is_simple_inorganic(smiles):
            return 'inorganic'

        elif self.has_carbon(mol):
            if self.has_carbon_hydrogen_bond(mol, smiles) or self.has_carbon_halogen_bond(mol) or self.has_carbon_organic_element_bond(mol):
                return 'organic'
        
        return 'inorganic'
    
    def is_simple_inorganic(self, smiles: str) -> bool:
        """
            Checks if the molecule is a simple inorganic carbon compound like CO2 or carbonate.
            
            :param smiles: SMILES string
            :return: True if the molecule is a simple inorganic carbon compound, False otherwise
        """

        simple_inorganic =  ['O=C=O', 'C(=O)=O', '[C+]#[O+]', '[O-]C(=O)[O-]', 'CO', 'CO2']

        return smiles in simple_inorganic

    def has_carbon(self, mol: Chem.Mol) -> bool:
        """
            Checks if the molecule contains carbon atoms.
            
            :param mol: RDKit Mol object
            :return bool: True if carbon is present, False otherwise
        """
        
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    def has_carbon_smiles(self, smiles: str) -> bool:
        """
            Checks if a SMILES string contains carbon atoms.

            Uses regex to detect both aliphatic ('C') and aromatic ('c') carbons,
            while avoiding false positives from other elements containing 'C' or 'c'.

            :param smiles: SMILES representation of the molecule.
            :return bool: True if carbon is present, False otherwise.

            Example:
            --------
            >>> has_carbon_smiles("CCO")
            True
            >>> has_carbon_smiles("NaCl")
            False
        """

        C_upper = re.compile(r'C[^saeroudnfl]')
        c_lower = re.compile(r'[^SATM]c')

        C_patterns = [C_upper, c_lower]

        return any(re.search(pattern, smiles) for pattern in C_patterns)


    def has_carbon_organic_element_bond(self, mol: Chem.Mol) -> bool:
        """
            Checks if the molecule contains a bond between carbon and other organic elements (N, P, S).

            :param mol: RDKit Mol object
            :return: True if a carbon-organic element bond is present, False otherwise
        """

        organic_elements = {7, 15, 16}  # N, P, S

        return any(bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() in organic_elements
                   or bond.GetEndAtom().GetAtomicNum() == 6 and bond.GetBeginAtom().GetAtomicNum() in organic_elements
                   for bond in mol.GetBonds())

    def has_carbon_organic_element_bond_smiles(self, smiles: str) -> bool:
        """
            Checks if the molecule contains a bond between carbon and other organic elements (N, P, S) using SMILES string.

            :param smiles: SMILES string of the molecule
            :return: True if a carbon-organic element bond is present in an organic compound, False otherwise
        """

        # Match C-N, C-P, or C-S bonds, but only in organic context
        
        n_pattern = re.compile(r'N[^aepbdmoi]')  # Matches N not followed by a, e, p, b, d, m, o, i
        p_pattern = re.compile(r'P[^dabmourt]')  # Matches P not followed by d, a, b, m, o, u, r, t
        s_pattern = re.compile(r'S[^icnmrebf]')  # Matches S not followed by i, c, n, m, r, e, b, f

        patterns = [n_pattern, p_pattern, s_pattern]

        return any(re.search(pattern, smiles) for pattern in patterns)

    def has_carbon_hydrogen_bond(self, mol: Chem.Mol, smiles: str) -> bool:
        """
            Checks for the presence of carbon-hydrogen bonds in a molecule.
            
            :param mol: RDKit Mol object
            :return: True if C-H bonds are found, False otherwise
        """

        mol_with_h = Chem.AddHs(mol)

        return any(bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 1
                or bond.GetEndAtom().GetAtomicNum() == 6 and bond.GetBeginAtom().GetAtomicNum() == 1
                for bond in mol_with_h.GetBonds())

    def has_carbon_hydrogen_bond_smiles(self, smiles: str) -> bool:
        """
            Checks for the presence of carbon-hydrogen bonds in a molecule using SMILES string.
            
            :param smiles: SMILES string of the molecule
            :return: True if C-H bonds are found, False otherwise
        """

        # Check for explicit hydrogens
        if re.search(r'C[H]', smiles):
            return True
        
        # Hydrogen check
        hs_pattern_1 = re.compile(r'\(\[H\]\).?[Cc]')
        hs_pattern_2 = re.compile(r'\[H\].?[Cc]')
        hs_pattern_3 = re.compile(r'[Cc]\(\[H\]\)')
        patterns = [hs_pattern_1, hs_pattern_2, hs_pattern_3]

        # # Check for implicit hydrogens (carbons with fewer than 4 bonds)
        # carbon_pattern = r'C(?![H])(?:(?![\(\)\[\]\.#=:])[^\d]){0,3}(?=[\(\)\[\]\.#=:]|$)'

        return any(re.search(pattern, smiles) for pattern in patterns)

    def has_carbon_halogen_bond(self, mol: Chem.Mol) -> bool:
        """
            Checks if the molecule contains a carbon-halogen bond.

            :param mol: RDKit Mol object
            :return: True if a carbon-halogen bond is present, False otherwise
        """

        halogen_elements = {9, 17, 35, 53} # F, Cl, Br, I

        return any(bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() in halogen_elements
                   or bond.GetEndAtom().GetAtomicNum() == 6 and bond.GetBeginAtom().GetAtomicNum() in halogen_elements
                   for bond in mol.GetBonds())
    
    def has_carbon_halogen_bond_smiles(self, smiles: str) -> bool:
        """
            Checks if the molecule contains a carbon-halogen bond using SMILES string.

            :param smiles: SMILES string of the molecule
            :return: True if a carbon-halogen bond is present, False otherwise
        """

        halogen_elements = ['Cl','Br']
        i_pattern = re.compile(r'I[^rn]')
        f_pattern = re.compile(r'F[^rem]')

        return any(element in smiles for element in halogen_elements or re.search(i_pattern, smiles) or re.search(f_pattern, smiles))

    # def has_carbon_hydrogen_bond(self, molecule_object: Chem.Mol) -> Union[bool,str]:
    #     """
    #         This function checks the presence of Hydrogen atoms in a molecule with
    #         Carbons.

    #         :param molecule_object:

    #         :return h_check:
    #     """

    #     # Hydrogen check
    #     hs_pattern_1 = re.compile(r'\(\[H\]\).?[Cc]')
    #     hs_pattern_2 = re.compile(r'\[H\].?[Cc]')
    #     hs_pattern_3 = re.compile(r'[Cc]\(\[H\]\)')

    #     mol_hs = Chem.AddHs(molecule_object)
    #     smi_hs = Chem.MolToSmiles(mol_hs)
    #     if re.search(hs_pattern_1, smi_hs) or re.search(hs_pattern_2, smi_hs) or re.search(hs_pattern_3, smi_hs):
    #         # This pattern here will always be true since C or c will be present even if there is no H surrounding them. 
    #         # Need to check that at least ONE H is near the C
    #         h_check = 'organic'
    #     else:
    #         h_check = False
        
    #     return h_check

    # def halogen_check(self, molecule_string: str) -> str:
    #     """
    #         This function checks the presence of halogen elements in the substance
    #         once the Carbon atom has been identified.

    #         :param molecule_string: SMILES string of the substance

    #         :return hal_check
    #     """

    #     # halogen check
    #     hal_check = 'inorganic'
    #     halogen_elements = ['Cl','Br']
    #     i_pattern = re.compile(r'I[^rn]')
    #     f_pattern = re.compile(r'F[^rem]')
        
    #     if re.search(i_pattern, molecule_string) or re.search(f_pattern, molecule_string):
    #         hal_check = 'organic'
        
    #     if hal_check == 'inorganic':
    #         for element in halogen_elements:
    #             if element in molecule_string:
    #                 hal_check = 'organic'

    #     return hal_check

    def check_organometallic(self) -> Optional[str]:
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

            :param molecule: SMILES string

            :return peptide_mol: string returning peptide. None otherwise.
        """

        peptide_mol = None

        l_peptide_pattern = re.compile(r'N\[C@@H\].+C\(=O\)O')
        r_peptide_pattern = re.compile(r'N\[C@H\].+C\(=O\)O')
        patterns = [l_peptide_pattern, r_peptide_pattern]

        if any(re.search(pattern, molecule) for pattern in patterns):
            peptide_mol = 'peptide'
        
        return peptide_mol

    def check_inorganic_metal(self, molecule: str) -> Optional[str]:
        """
            Checks if molecule is only an elemental metal.

            :param molecule: SMILES string

            :return metal: returns inorganic_metal if SMILES is only a metal. None otherwise.
        """

        metal = None
        
        if molecule in ps._metals:
            metal = 'inorganic_metal'
        
        return metal

    def check_salt(self, subType: str) -> str:
        """
            Checks if the molecule is salt.

            :param molecule:

            :return salt:
        """

        remover = SaltRemover()
        salt = None
        
        res, deleted = remover.StripMolWithDeleted(self.smiles_mol)

        if len(deleted) >= 1:
            # print(Chem.MolToSmiles(res))
            # print([Chem.MolToSmiles(residue) for residue in deleted])
            salt = '_'.join([subType,'salt'])

        return salt
    
    def salt_remover(self, mol_object: Chem.Mol) -> str:
        """
            Removes salts and counterions. Non sanitizable molecules can't be processed

            :param smiles: smiles string

            :return cleaned_smiles: 
        """
        
        # Always apply LargestFragmentChooser, even if there's no "."
        remover = rdMolStandardize.LargestFragmentChooser()
        largest_mol = remover.choose(mol_object)

        # Standardize the molecule
        standardized_mol = standardizer.standardize_mol(largest_mol)

         # Convert back to SMILES
        cleaned_smiles = Chem.MolToSmiles(standardized_mol, isomericSmiles=True)

        # if "." in smiles and Chem.MolFromSmiles(smiles):
        #     cleaned_smiles = Chem.MolToSmiles(rmv.choose(self.smiles_mol))
        # else:
        #     cleaned_smiles = smiles

        return cleaned_smiles