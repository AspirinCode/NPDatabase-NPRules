#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 16:32:37 2018

Script that extracts all structure data from a given superclass/ class/
subclass present in NPDatabase and writes it into a txt file which can 
be used as molBLOCKS input. NPDatabase_path should be specified in this
script.

Usage:
Command line: python3 create_molblocks_input.py [taxonomy] [taxonomy_name]
taxonomy: superclass, class or subclass
taxonomy_name: the exact name of superclass, class or subclass, including
capitals and punctuation marks.
Example:
Command line: python3 create_molblocks_input.py class Polypeptides
@author: stokm006
"""

import sys
import sqlite3

def store_structures(taxonomy, taxonomy_name, NPDatabase_path):
    """
    function that gets all structure data and writes it into a file suitable
    as input for molBLOCKS
    
    
    taxonomy: str, superclass, class or subclass
    taxonomy_name: str, name of taxonomy
    NPDatabasePath: path to NPDatabase
    
    """
    # Connect to NPDatabase
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Generate the nr of structures present in given superclass/ class/ subclass
    numrows = c.execute("SELECT count(*) FROM structure where cf_{0} = '{1}';".format(taxonomy, taxonomy_name))
    numrows = str(numrows.fetchone()).lstrip("(").rstrip(",)")
    numrows = int(numrows)   
    
    # Select all structure data from superclass/ class/ subclass
    c.execute("SELECT * FROM structure where cf_{0} = '{1}';".format(taxonomy, taxonomy_name))
    
    structure_id_list = []    
    canonical_smiles_list = []
    per_struct_list = []
    for x in range(0, numrows):
        row = c.fetchone()
        canonical_smiles = row[2]    
        # Get only unique canonical SMILES and structure id's
        if canonical_smiles not in canonical_smiles_list:
            structure_id_list += [row[0]]
            canonical_smiles_list += [row[2]]
    
    
    # Write structures into a txt file
    output_file_name = taxonomy+  '_' + taxonomy_name + '_' + 'molBLOCKS_input.txt'
    with open(output_file_name, 'w') as db_file:
        for i in range(len(canonical_smiles_list)):
            db_file.write(canonical_smiles_list[i] + '\t' + structure_id_list[i] + '\n')
    
    print ("----------------------------------------\n")
    print (len(canonical_smiles_list), taxonomy_name, 'structures written into', output_file_name, 'file\n')  

        
if __name__ == '__main__':
   
   NPDatabase_path = "/path/to/NPDatabase.sqlite"
   argv = sys.argv[1:]
   taxonomy = argv[0]
   taxonomy_name = argv[1]
    
   store_structures(taxonomy, taxonomy_name, NPDatabase_path)
