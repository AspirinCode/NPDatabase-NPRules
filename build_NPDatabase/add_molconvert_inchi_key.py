"""
Script that adds the molconvert inchi key generated with Rutger's
pipeline to NPDatabase.

Usage:
Command line: python3 add_molconvert_inchi_key.py [molconvert_output] [path/to/NPDatabasefile.sqlite]

"""

import sys
import sqlite3

def add_molconvert_inchi_keys(inchikey_data, NPDatabase_path):
    """
    This script add all inchi keys to the structure table in the inchi_key_molconvet
    table
   
    inchikey_data: txt file with structure id, isomeric SMILES and molconvert
    inchi key
    NPDatabase_path: path to NPDatabase    
    """
  
    with open(inchikey_data) as file_object:
        input_file = file_object.read()  

    all_lines = input_file.split('\n')
    inchikey_dict = {}
    for line in all_lines[0:-1]:
        line = line.split(' ')
        structure_id = line[0]
        molconvert_inchikey = line[2]
        molconvert_inchikey = molconvert_inchikey[9:]
        inchikey_dict[structure_id]=[molconvert_inchikey]

    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
        
    sql_list = [] 
    for identifier in inchikey_dict:
        inchi_key = str(inchikey_dict[identifier])
        inchi_key = inchi_key.lstrip("['").rstrip("']")
        sql_string = "update structure set inchi_key_molconvert = '%s' where structure_id = '%s';"% (inchi_key, identifier)
        sql_list += [sql_string] 
        
        
    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
  
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()



if __name__ == '__main__':

    
    inchi_key_input = sys.argv[1]
    NPDatabase_path = sys.argv[2]
    
    add_molconvert_inchi_keys(inchi_key_input, NPDatabase_path)

