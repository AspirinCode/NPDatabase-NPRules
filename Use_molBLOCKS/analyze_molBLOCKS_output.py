#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 09:52:34 2018

Script that processes molBLOCKS output, calculates and prints several
substructure results.

Use RDkit version 2018.09.2

Usage:
Command line: python3 analyze_molBLOCKS_output.py [output_file_molBLOCKS]

@author: stokm006
"""


from sys import argv
from rdkit import Chem 
from rdkit import RDLogger


def calculate_structure_results(molblocks_results):
    """
    This function calculates and prints the number of substructures generated
    per structure, the nr of input structures, the total number of output
    substructures and the number of unique substructures.
    
    molblocks_results: txt file with substructure smiles and structure_id
    """
        

    substructure_results = molblocks_results.split('\n')
    nr_of_input_structures = 0
    rdkit_smiles_list = [] 
    subs_per_structure_list = []
      

    print ('\nNr.\t structure_id\t\t Nr. of substructures')
    for i, subs in enumerate(substructure_results[:-1]):
        subs = subs.split('\t')
        structure_id = subs[1]
        subs = subs[0].split('.')
        
        if subs[0] == '<NA>':
            nr_of_substructures = '0'
        if subs[0] != '<NA>':
            nr_of_substructures = len(subs)
        print (i,'\t',structure_id,'\t\t', nr_of_substructures)
        # to prevent warning print statements from RDkit 
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        
        nr_of_input_structures += 1
        subs_per_structure_list += [subs]
        for sub in subs:
            mol = Chem.MolFromSmiles(sub)
            rdkit_smiles_list += [Chem.MolToSmiles(mol)]  # RDKit canonical SMILES for consistency

    print ('\n          ~~~~~TOTAL~~~~~')
    print ('Total nr. of input structures      ', nr_of_input_structures)
    
    # Generate the total number of (unique) substructures
    unique_sub_list = []
    total_nr_of_subs = 0
    for smiles in rdkit_smiles_list:
        if smiles != '<NA>' and smiles not in unique_sub_list:
            unique_sub_list += [smiles]
        if smiles != '<NA>':
            total_nr_of_subs += 1
    
    print ('Total nr. of output substructures  ', total_nr_of_subs)
    print ('Total nr. of unique substructures  ', len(unique_sub_list))
    
    # Generate the number of fragmented structures in order to calculate the % of fragmented structures
    nr_of_fragmented_structures = 0
    for subs in subs_per_structure_list:
        if len (subs) >= 1 and subs[0] != '<NA>':
            nr_of_fragmented_structures += 1
            

    print ('Fragmented input structures        ', nr_of_fragmented_structures/nr_of_input_structures*100, '%')



if __name__ == "__main__":
    with open(argv[1]) as file_object:
        molblocks_results = file_object.read()
        calculate_structure_results(molblocks_results)

