
"""
Input: parsed file with structure data from new_input.py

Checkes whether the new structure is already in NPDatabase, if not, then 
the structure is added to the output file. The external_source_id's and 
isomeric smiles (csv) are written into the file which can be used as input
for molconvert to create inchi keys, see https://github.com/NP-Plug-and-Play-Scripts/inchiKeyCreatorPipeline
    
Use RDkit version 2018.09.2

Usage:
Command line: python3 create_molconvert_input.py [parsed new input file] [path/to/NPDatabasefile.sqlite]

"""
  
  
from rdkit import Chem
from sys import argv
import sqlite3


def create_molconvert_input(input_file, input_file_name, NPDatabase_path):
    """
    Checkes whether the new structure is already in NPDatabase, if not,
    then the structure is added to the output file with the external source id 
    and the isomeric SMILES.
    
    input_file: new parsed input data with external source id and isomeric smiles
    input_file_name: name that will be used for the molconvert input file
    NPDatabase_path: path to NPDatabase
    """
    
    input_file_name = input_file_name[0:-13]
    input_file_name = input_file_name + "_molconvert_input.csv"
  
    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Extract all smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure')
    iso_smiles_sqlite_list = []
    for row in c:
        iso_smiles_sqlite_list += [row[1]]

    not_yet_in_db_count = 0
    in_db_count = 0
    
    all_lines = input_file.split('\n')
    for i, line in enumerate(all_lines[1:-1]):
        line = line.split('|||')
        iso_smiles_mibig = line[1]
        external_source_id = line[0]
        
        # if isomeric smiles not yet in db write it into a csv file for molconvert
        temp_structure_has_new_structure_id_list = []
        if iso_smiles_mibig not in iso_smiles_sqlite_list and "*" not in iso_smiles_mibig:  # molconvert cannot proccess SMILES having "*"
            not_yet_in_db_count += 1

            with open(input_file_name, 'a') as db_file:
                db_file.write(external_source_id + ',' + iso_smiles_mibig + '\n')      

        if iso_smiles_mibig in iso_smiles_sqlite_list:
            in_db_count += 1
               
    print (not_yet_in_db_count, "structures not present in NPDatabase") 
    print ("external source id's and isomeric SMILES written in a csv file suitable for molconvert pipeline")
    print (in_db_count, "structures already present in NPDatabase")
     
     
if __name__ == "__main__":
    
    input_file_name = argv[1]
    NPDatabase_path = argv[2]
    
    
    with open(argv[1]) as file_object:

        input_file = file_object.read()
        create_molconvert_input(input_file, input_file_name, NPDatabase_path)
   
   
