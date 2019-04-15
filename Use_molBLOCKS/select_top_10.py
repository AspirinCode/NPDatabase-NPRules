#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29-03-19

Script that selects the top 10 most occurring substructures (molBLOCKS output),
writes the results into a csv file and shows the substructures in a 
pop-up display

This scripts requires ImageMagick, see: https://anaconda.org/conda-forge/imagemagick
Use RDkit version 2018.09.2

Usage:
Command line: python3 select_top_10.py [output_file_molBLOCKS]

@author: stokm006
"""


from sys import argv
from rdkit import Chem 
from rdkit import RDLogger
from rdkit.Chem import Draw 
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG = True


def select_top10(molblocks_results, output_name):
    """
    This function first selects the top 10 most occurring substructures
    from the molBLOCKS output file. The SMILES are written into a csv 
    file together with the total number of that substructure (frequency 
    in molBLOCKS output file) and the percentage. Also the top 10 
    substructures are shown.
    
    
    molblocks_results: txt file with substructure smiles and structure_id
    output_name: str, name for the output csv file
    """
        
    substructure_results = molblocks_results.split('\n')
    smiles_list = []
  
    # Parse data and generate rdkit canonical smiles 
    for subs in substructure_results[0:-1]:
        subs = subs.split('\t')
        subs = subs[0].split('.')

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)

        for sub in subs:
            if sub != '<NA>':
                mol = Chem.MolFromSmiles(sub)
                rdkit_smiles_list += [Chem.MolToSmiles(mol)]
    total_nr_of_subs = len(smiles_list)
    
    # Determine the total number of each unique substructure
    hit_list = []
    top_dict = {}
    unique_sub_list = []
    for smiles in rdkit_smiles_list:
        count_smiles = smiles_list.count(smiles)
        percentage = count_smiles/total_nr_of_subs*100
        percentage = round(percentage, 4)
        if smiles not in hit_list:
            hit_list += [smiles, count_smiles, percentage]

    hit_list2 = []
    for i in range(0,(len(hit_list)), 3):
        hit_list2 += [[hit_list[i],hit_list[i+1],hit_list[i+2]]]
    
    # Sort list to get the top 10 ASC
    sorted_hit_list2 = sorted(hit_list2, key=lambda x: int(x[1]), reverse=True)
   
    
    # Writes data top 10 into a csv file and creates draw mol list
    mol_draw_list = []
    with open(output_name, 'w') as db_file:
        db_file.write('Nr,SubstructureSMILES,Total,Percentage\n')
        for i, hit in enumerate(sorted_hit_list2[0:10]):
             sub_smiles = hit[0]
            mol_draw_list += [Chem.MolFromSmiles(sub_smiles)]
            count = hit[1]
            percentage = hit[2]
            db_file.write(str(i+1)+ ','+str(sub_smiles)+','+str(count)+','+str(percentage)+'\n')
    
    
    # Shows the top 10 structures in a pop-up display
    im = Draw.MolsToGridImage(mol_draw_list, 
                                  useSVG=False,
                                  legends=['1','2','3','4','5','6','7', '8','9','10'],
                                  molsPerRow=5,
                                  subImgSize = (250,300))
    im.show()



if __name__ == "__main__":
    output_name = argv[1]
    output_name = output_name[0:-4]
    output_name = output_name + "_top_hits.csv"
    
    with open(argv[1]) as file_object:
        result_substructures = file_object.read()
        select_top10(molblocks_results, output_name)

