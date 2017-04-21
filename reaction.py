
# coding: utf-8

# In[1]:

### -----  Oxidation of primary and secondary alcohols to caboxylic acid ----  ####

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('CCCCS') #input molecule ----butanol
sub_str = Chem.MolFromSmarts('[CX4;!H0][OH1]')
if m.HasSubstructMatch(sub_str):#substructure
    oxidation_state = input("What is the desired oxidation state? Type 1 for aldehyde/ketone and 2 for carboxylic acid.")

    if oxidation_state == 1:
        rxn = AllChem.ReactionFromSmarts('[CX4;!H0:1][OH1:2]>>[C:1]=[O:2]') #aldehyde or keton
        product = rxn.RunReactants((Chem.MolFromSmiles('CCCCS'),))
        print Chem.MolToSmiles(product[0][0])
        
    elif oxidation_state == 2:
        rxn = AllChem.ReactionFromSmarts('[CX4;!H0:1][OH1:3]>>[OXH:3][CX4:1]=[O:2]') #carboxylic acid
        product = rxn.RunReactants((Chem.MolFromSmiles('CCCCS'),))
        print Chem.MolToSmiles(product[0][0])
    else:
        print "Error. Please pick 1 or 2."
    
else:
    print "Error! This is not the right substructure. No further oxidation possible. How is your chemistry?"


# In[2]:

### -----  Thioether formation ----  ####

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('CCS') #input molecule --- thioether
sub_str = Chem.MolFromSmarts('[C,c][SH1]')
if m.HasSubstructMatch(sub_str):#substructure
    rxn = AllChem.ReactionFromSmarts('[C,c:1][SH1:2].[CX4;!H0:3][Br:4]>>[C:1][S:2][C:3].[Br:4]') #thioether
    product = rxn.RunReactants((Chem.MolFromSmiles('CCS'), Chem.MolFromSmiles('CBr')))
    print Chem.MolToSmiles(product[0][0])

else:
    print "Error! This is not the right substructure. No thioether formation is possible. How is your chemistry?"


# In[3]:

### -----  Disulfide formation ----  ####

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('CCCCS') #input molecule --- thioether
sub_str = Chem.MolFromSmarts('[C,c][SH1]')
if m.HasSubstructMatch(sub_str):#substructure
    rxn = AllChem.ReactionFromSmarts('[C,c:1][SH1:2].[CX4;H2,H1:4][SH1:3]>>[C,c:1][S:2][S:3][C:4]') #disulfide
    product = rxn.RunReactants((Chem.MolFromSmiles('CCCS'), Chem.MolFromSmiles('CCCS')))
    print Chem.MolToSmiles(product[0][0])

else:
    print "Error! This is not the right substructure. No disulfide formation is possible. How is your chemistry?"


# In[4]:

### -----  Sulfoxide formation ----  ####

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('CSC') #input molecule --- thioether
sub_str = Chem.MolFromSmarts('[SX2][C]')
if m.HasSubstructMatch(sub_str):#substructure
    rxn = AllChem.ReactionFromSmarts('[SX2:1][C:2]>>[C:2][S:1]=[O:3]') #disulfide
    product = rxn.RunReactants((Chem.MolFromSmiles('CSC'),))
    print Chem.MolToSmiles(product[0][0])

else:
    print "Error! This is not the right substructure. No disulfide formation is possible. How is your chemistry?"

