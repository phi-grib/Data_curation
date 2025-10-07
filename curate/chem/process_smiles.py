import os

from rdkit import Chem
from standardiser2 import standardise
timeout = -1

_metal_nof = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,n,O,o,F]')
_metal_non = Chem.MolFromSmarts('[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,c,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]')
_metals = ['[Al]','[Sc]','[Ti]','[V]','[Cr]','[Mn]','[Fe]','[Co]','[Ni]','[Cu]','[Zn]','[Y]','[Zr]','[Nb]','[Mo]','[Tc]',
'[Ru]','[Rh]','[Pd]','[Pd++]','[Pd+2]','[Ag]','[Cd]','[Cd++]','[Hf]','[Ta]','[W]','[Re]','[Os]','[Ir]','[Pt]','[Au]','[Sn]',
'[Pb]','[Hg]','[Cd+2]','[Cr+3]','[Cr+6]','[Sn+3]','[Mg++]','[Sb+3]','[Al+3]','[Ba++]','[Ba+2]','[Fe+3]','[Mg+2]','[B+3]','[B]',
'[Pb++]','[Pb+2]','[Zr-2]','[Mn+2]','[Mn+3]','[Sb]','[Ti+4]','[Fe+2]','[Fe++]','[Rb+]','[Mg]','[Ca++]','[Ca+2]','[As]','[Si]','[Ge]',
'[Te]','[Bi]','[Cu+]','[Se]','[se]']

def disconnect(mol):
    """
        Adapated from molVS standardizer module. Now it returns the list of metals it has disconnected
    """

    metals = set([])
    for metal_atom in _metals:
        rwmol = Chem.RWMol(mol)
        smarts = Chem.MolFromSmarts(metal_atom)
        pairs = rwmol.GetSubstructMatches(smarts)
        for i, in reversed(pairs):
            metalSymbol = rwmol.GetAtomWithIdx(i).GetSmarts()
            metals.add(metalSymbol)
            rwmol.RemoveAtom(i)
        mol = rwmol.GetMol()

    for smarts in [_metal_nof, _metal_non]:
        pairs = mol.GetSubstructMatches(smarts)
        rwmol = Chem.RWMol(mol)
        orders = []
        for i, j in reversed(pairs):
            metalSymbol = mol.GetAtomWithIdx(i).GetSmarts()
            metals.add(metalSymbol)
            orders.append(int(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble()))
            rwmol.RemoveBond(i, j)
        # Adjust neighbouring charges accordingly
        mol = rwmol.GetMol()
        for n, (i, j) in enumerate(pairs):
            chg = orders[n]
            atom1 = mol.GetAtomWithIdx(i)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + chg)
            atom2 = mol.GetAtomWithIdx(j)
            atom2.SetFormalCharge(atom2.GetFormalCharge() - chg)
    
    return mol, metals

def normalize(inF, outF, singleF, failedF, remove_salts= True, keep_nonorganic= False, verbose=False, pH=7.4) :
      
    count = 0        ## count for the whole dataset
    count_inc = 0    ## count for only included molecules
    count_exc = 0    ## count for only excluded molecules
    all_salts = 0    ## count for entries with only salts / solvent
    fail_sanity = 0  ## count for entries that fail sanity check 
    fail_mol = 0     ## count for entries that fail to create mol object 
    fail_prot = 0    ## count for entries that fail protonation

    header = '%s\n' %('\t'.join(['CAS', 'Component', 'Original smiles', 'smiles']))
    fail_header = '%s\n' %('\t'.join(['CAS', 'Original smiles', 'Error']))

    outF.write(header)
    singleF.write(header)
    failedF.write(fail_header)
    
    for line in inF:
        count += 1
        try:
            cas, smi = line.rstrip().split('\t')
        except:
            print ('Failed parsing line:')
            print (line)
            failedF.write(line.rstrip()+'\tFailed parsing line\n')
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            count_exc += 1
            fail_mol += 1
            failedF.write(line.rstrip()+'\tFailed to create molecule object\n')
            continue

        try:
            #mol = standardise.run(mol, keep_nonorganic= keep_nonorganic, remove_salts= remove_salts)
            succ, mol, err = standardise.run(mol, keep_nonorganic= keep_nonorganic)
        except Exception as err:
            err = '{}'.format(err)
            count_exc += 1
            fail_sanity += 1
            failedF.write('{}\t{}\t{}\n'.format(cas, smi, err))
            continue

        i = 1
        if succ:
            count_inc += 1
            nHA = mol.GetNumHeavyAtoms()
            if nHA < 2:
                singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(mol, isomericSmiles=True)))
            else:
                outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(mol, isomericSmiles=True)))
                #prot, protMol = protonate(mol, pH)
                #if prot:
                #    outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(protMol, isomericSmiles=True)))
                #else:
                #    failedF.write('{}\t{}\t{}\n'.format(cas, smi, protMol))
                #    fail_prot += 1
        else:
            smis = set([Chem.MolToSmiles(moli, isomericSmiles=True) for moli in mol])
            if err == 'Multiple non-salt/solvate components':
                for smii in smis:
                    moli = Chem.MolFromSmiles(smii)
                    nHA = moli.GetNumHeavyAtoms()
                    if nHA < 2:
                        singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                    else:
                        outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                        #prot, protMol = protonate(Chem.MolFromSmiles(smii), pH)
                        #if prot:
                        #    outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(protMol, isomericSmiles=True)))
                        #else:
                        #    failedF.write('{}\t{}\t{}\n'.format(cas, smi, protMol))
                        #    fail_prot += 1
                    i += 1
                count_inc += 1
            elif err == 'No non-salt/solvate components':
                metal = False
                for smii in smis:
                    moli = Chem.MolFromSmiles(smii)
                    nHA = moli.GetNumHeavyAtoms()
                    if nHA == 1 and moli.GetAtomWithIdx(0).GetSymbol() in _metals:
                        singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                        metal = True
                        i += 1
                if metal:
                    count_inc += 1
                else:
                    count_exc += 1
                    all_salts += 1
                    failedF.write('{}\t{}\t{}\n'.format(cas, smi, err))
    
    os.system('rm in.sdf out.sdf')
    print ('the full dataset = {}'.format(count))
    print ('Molecules normalized = {}'.format(count_inc))
    print ('Molecules excluded = {}'.format(count_exc))
    print ('   Fail RDkit mol object = {}'.format(fail_mol))
    print ('   Fail protonation = {}'.format(fail_prot))
    print ('   Fail sanity check = {}'.format(fail_sanity))
    print ('   Only salts / solvent = {}'.format(all_salts))

def std(mol, returnMetals=False):
    # Standardize and return a dictionary with the smiles as keys
    # and the molecule object and whether it's a metal ion as values
    stdD = {}

    # Check single atom compounds, to see if they are metal ions
    if mol.GetNumAtoms() == 1:
        at = mol.GetAtoms()[0].GetSymbol()
        symbol = '[%s]' %at
        if at in _metals and returnMetals:
            cmpd = Chem.MolFromSmiles(symbol)
            stdD[symbol] = (cmpd, True, True, '')
        else:
            (passed, std_cmpd, errmessage) = standardise.run(mol)
            if passed:
                stdD[symbol] = (std_cmpd, False, passed, errmessage)
    else:  
        # Extract metal ions from complex compounds
        comp_mol, metals = disconnect(mol)
        if returnMetals:
            for metal in metals:
                metalmol = Chem.MolFromSmiles(metal)
                metal = '[%s]' %metalmol.GetAtoms()[0].GetSymbol()
                cmpd = Chem.MolFromSmiles(metal)
                stdD[metal] = (cmpd, True, True, '')

        # For the rest of the molecule, standardize and add
        standardise.run(comp_mol)
        try:
            (passed, std_cmpds, errmessage) = standardise.run(comp_mol)
        except:
            passed = False
            errmessage = 'Failed'

        if passed:
            stdD[Chem.MolToSmiles(std_cmpds, isomericSmiles=True)] = (std_cmpds, False, True, '')
        elif errmessage == 'Multiple non-salt/solvate components':
            cmpdD = {}
            for cmpd in std_cmpds:
                inchi = Chem.MolToInchi(cmpd)
                cmpdD[inchi] = cmpd

            for inchi in cmpdD:
                cmpd = cmpdD[inchi]
                stdD[Chem.MolToSmiles(cmpd, isomericSmiles=True)] = (cmpd, False, True, 'Multiple non-salt/solvate components')
        else:
            stdD[Chem.MolToSmiles(mol, isomericSmiles=True)] = (mol, False, False, errmessage)
                
    return stdD