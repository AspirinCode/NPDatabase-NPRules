"""
This script generates all structure data and stores it in NPDatabase.
Structure identifiers and all attributes generated with rdkit are included.

Usage: python3 create_NPDatabase.py [all_collected_structure_data] [path/to/NPDatabase.sqlite]
Note: A new sqlite file is generated automatically

Use RDkit version 2018.09.2

Example:
Command line: python3 create_NPDatabase.py all_structure_input.txt /mnt/nexenta/stokm006/NPDatabase.sqlite
"""


import sys, getopt
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdinchi
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)
import sqlite3

def create_canonical_data(input_file):
    """
    This function transforms all structures from the combined database file
    into the isomeric and canonical (NOTE: use the last RDkit version from 2019) 
    and returns them into a list togetheter with the external source id and 
    db_source name
    
    
    input_file: txt file with all structure SMILES, external source id's
    and database name (from collect_all_structures.py)
    
    """

    all_lines = input_file.split('\n')

    fails = 0
   
    canonical_SMILES_list = []
    identifier_list = []
    db_name_list = []
    
    all_info_list = []
    
    for line in all_lines[1:-1]:
        line = line.split('\t')
        SMILES = line[0]
        identifier = line[1]
        db_name = line[2]
        mol_SMILES = Chem.MolFromSmiles(SMILES)
        str_mol_SMILES = str(mol_SMILES)
        all_info_list_temp = []
        if str_mol_SMILES.startswith('<'):  
            canonical_SMILES = Chem.MolToSmiles(mol_SMILES)
            all_info_list_temp += [canonical_SMILES, identifier, db_name]
            all_info_list += [all_info_list_temp]

   
    # The same smiles are grouped together in the list
    canonical_SMILES_list_sorted = sorted(all_info_list, key = lambda x:(x[0])) # key tells which attribute, in this case the smiles, will be sorted
    
    return (canonical_SMILES_list_sorted)


def create_structure_data(canonical_smiles_list):
    """
   This function creates structure_id's, canonical smiles (non isomeric),
   molecular weight, molecular formula, InChI and InChIkey (RDkit)
   
   
   canonical_smiles_list: list of lists with isomeric canonical smiles, external source id's
   and source names
    """


    
    SMILES_list = []
    for line in canonical_smiles_list:
        SMILES = line[0]
        if SMILES not in SMILES_list:
            SMILES_list += [SMILES]
    
    all_RDkit_structure_data = []
    for i, SMILES in enumerate(SMILES_list):
       
        #Create unique identifier (structure_id)              
        structure_id = format(i+1, "07")
        structure_id = 'NPDB:{}'.format(structure_id)
        temp_all_RDkit_structure_data = [structure_id]        
        
        #Isomeric and canonical SMILES
        temp_all_RDkit_structure_data += [SMILES]
        mol = Chem.MolFromSmiles(SMILES)

        #Canonical SMILES list
        can_SMILES = Chem.MolToSmiles(mol, isomericSmiles = False)
        temp_all_RDkit_structure_data += [can_SMILES]

        #Monoisotopic mass 
        mol_weigth = Descriptors.ExactMolWt(mol)
        temp_all_RDkit_structure_data += [mol_weigth]
        #Molecular formula
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        temp_all_RDkit_structure_data += [mol_formula]
        #InChI         
        inchi = Chem.MolToInchi(mol)
        temp_all_RDkit_structure_data += [inchi]  
        #InChIKey
        inchi_key = Chem.MolToInchiKey(mol)
        temp_all_RDkit_structure_data += [inchi_key]
        all_RDkit_structure_data += [temp_all_RDkit_structure_data]

   
    return (all_RDkit_structure_data)

def create_structure_source_data(canonical_smiles_list, structure_table_data): 
    """
    This function generates all structure_has_source_table data

    canonical_smiles_list: list of lists with isomeric canonical smiles, external source id's
    and source names
    structure_table_data: list of lists from create_structure_data()
    """

    id_and_isoSMILES_list = []
    for line in structure_table_data:
        id_and_isoSMILES_list += [[line[0], line[1]]]    
    
    structure_id_list = []
    external_db_identifier_list = []
    db_name_list = []


    for line in canonical_smiles_list:
        SMILES = line[0]
        for id_SMILES in id_and_isoSMILES_list:
            if SMILES == id_SMILES[1]:
                stucture_id = id_SMILES[0]
                external_db_identifier = line[1]
                db_name = line[2]
                structure_id_list += [stucture_id]
                external_db_identifier_list += [external_db_identifier]
                db_name_list += [db_name]

    
    identifier_list = []
    for i, line in enumerate(external_db_identifier_list):
        structure_source_id = format(i+1, "08")             
        structure_source_id = 'NPSRC:{}'.format(structure_source_id)
        identifier_list += [structure_source_id]
        
    all_structure_source_data = []  
    for i in range(len(identifier_list)):
        all_structure_source_data += [[identifier_list[i], structure_id_list[i], external_db_identifier_list[i], db_name_list[i]]]
    
    return (all_structure_source_data)
 
def create_source_data(structure_source_data): 
    """
    Generates the data for the structure_source table
    
    structure_source_data: list of lists from create_structure_source_data()
    """

    db_name_dict = {}
    for line in structure_source_data:
        db_name = line[3]
        if db_name in db_name_dict:
             db_name_dict[db_name] += 1
        
        if db_name not in db_name_dict:
             db_name_dict[db_name] = 1
             
    data_source_list = []
    for line in db_name_dict.items():
        data_source_list += [[line[0], line[1]]]

    return (data_source_list)


def store_structure_table_data(structure_table_data, NPDatabase_path):
    """
    This function stores all structure data in the structure table in
    NPDatabase
    
    structure_table_data: list of lists from create_structure_data()
    NPDatabase_path: path to NPDatabase
    """


    sql_list = [] 
    for line in structure_table_data:
        sql_string = "INSERT INTO structure VALUES ('%s', '%s', '%s', '%s', '%s',\
        '%s', '%s');"% (line[0], line[1], line[2], line[3], line[4], \
        line[5], line[6])
        sql_list += [sql_string] 


    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
     
    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS structure")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE structure 
             (structure_id text PRIMARY KEY, isomeric_smiles text, canonical_smiles text, \
             monoisotopic_mass num, molecular_formula text, \
             inchi text, inchi_key_rdkit text)''')
    
    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
  
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()


def store_structure_has_data_source_table_data(structure_source_table_data, NPDatabase_path):
    """
    This function stores all structure source data in the structure_has_data_source
    table in NPDatabase
    
    structure_source_table_data: list of lists from create_structure_source_data
    NPDatabase_path: path to NPDatabase
    """



    sql_list = [] 
    for line in structure_source_table_data:
        sql_string = "INSERT INTO structure_has_data_source VALUES\
        ('%s','%s', '%s', '%s');"% (line[0], line[1], line[2], line[3])
        sql_list += [sql_string] 
  
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
     
    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS structure_has_data_source")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE structure_has_data_source
              (structure_source_id text PRIMARY KEY, structure_id text,\
              source_id text, source_name text,
              FOREIGN KEY(structure_id) REFERENCES structure(structure_id),
              FOREIGN KEY(source_id) REFERENCES data_source(source_id))''')

    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
 
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()


def store_source_table_data(source_table_data, NPDatabase_path):
    """
    This function stores all source data in the data_source table 
    in NPDatabase
    
    source_table_data: list of lists from create_data_source_data()
    NPDatabase_path: path to NPDatabase
    """

    sql_list = []
    for line in source_table_data:
        line = line
        sql_string = "INSERT INTO data_source VALUES ('%s', %s);"% (line[0], line[1])
        sql_list += [sql_string]
 

     # Connect
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
     
    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS data_source")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE data_source 
             (source_name text PRIMARY KEY, nr_of_structures num)''')


    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])

    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()
     
def main(argv):
    """
    Main function that calls all other functions 
    
    """
    
    input_file = argv[0]
    NPDatabase_path = argv[1]

   
    # generate all structure data
    with open(input_file) as file_object:
        input_file = file_object.read()  
    canonical_smiles_list = create_canonical_data(input_file)
    structure_table_data = create_structure_data(canonical_smiles_list)
    structure_source_table_data = create_structure_source_data(canonical_smiles_list, structure_table_data)
    source_table_data = create_source_data(structure_source_table_data)
    
    
    # store all structure data in NPDatabase 
    store_structure_table_data(structure_table_data, NPDatabase_path)
    store_structure_has_data_source_table_data(structure_source_table_data, NPDatabase_path)
    store_source_table_data(source_table_data, NPDatabase_path)

if __name__ == '__main__':
    
   main(sys.argv[1:])

