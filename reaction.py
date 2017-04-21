
# coding: utf-8

# In[1]:

### -----  Oxidation of primary and secondary alcohols to caboxylic acid ----  ####

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Take molecule

def mol_input(mol_str):
    mol = Chem.MolFromSmiles(mol_str)
    if mol is None:
        raise ValueError('Molecule is None from ', mol_str)
    return mol

# Check substructure
def collection_reagents():
    return {'[CX4;!H0][OH1]':['[CX4;!H0:1][OH1:2]>>[C:1]=[O:2]','[CX4;!H0:1][OH1:3]>>[OXH:3][CX4:1]=[O:2]'],
            #'[C,c][SH1]':['[C,c:1][SH1:2].[CX4;!H0:3][Br:4]>>[C:1][S:2][C:3].[Br:4]','[C,c:1][SH1:2].[CX4;H2,H1:4][SH1:3]>>[C,c:1][S:2][S:3][C:4]'],
            '[SX2][C]': ['[SX2:1][C:2]>>[C:2][S:1]=[O:3]']}


def check_for_substructure(mol):
    react_dict = collection_reagents()
    outlist = []
    for substructure in react_dict:
        mol_substr = Chem.MolFromSmarts(substructure)
        if mol.HasSubstructMatch(mol_substr):
            outlist.extend(react_dict[substructure])
    return outlist

# Reaction
def react(outlist, mol):
    prod_list = []
    for react_str in outlist:
        rxn = AllChem.ReactionFromSmarts(react_str)
        product = rxn.RunReactants((mol,))
        prod_list.append(Chem.MolToSmiles(product[0][0]))
    return prod_list

def react_str(mol_str):
    mol =  mol_input(mol_str)
    outlist = check_for_substructure(mol)
    prod_list = react(outlist, mol)
    return prod_list


if __name__ == "__main__":
    mol_str = sys.argv[1]
    prod_list = react_str(mol_str)
    for prod in prod_list:
        print prod