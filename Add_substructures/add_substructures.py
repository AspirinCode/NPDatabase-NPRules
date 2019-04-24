
"""
Created on 19-04-'19
Script that processes molBLOCKS output, generates new substructure data
and adds it to the substructure table. Also data for structure_has_substructure
is generated and stored.

Path to NPDatabase should be specified in this script.

Use RDkit version 2018.09.2

Output: NPDatabase with updated substructure and structure_has_substructure tables

Usage:
Command line: python3 add_substructures.py [input__file_molBLOCKS] [output_file_molBLOCKS]
@author: stokm006
"""

from sys import argv
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdinchi
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)
import sqlite3


def get_substructure_id(NPDatabase_path):
    """ 
    function that extracts and returns the substructure id with the highest
    number from NPDatabase
    
    NPDatabase_path: path to NPDatabase    
    """


    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Extract all substructure_id's from sqlite database
    c.execute('SELECT * from substructure')
    subtructure_id_sqlite_list = []

    for row in c:
        subtructure_id_sqlite_list += [row[0]]

    # For substructure_id's: delete "NPSB:" and only keep the nr
    substructure_id_list = []
    for structure_id in subtructure_id_sqlite_list:
        substructure_id_list += [structure_id[7:]]
    
    
    # Get the highest structure_source_id nr    
    substructure_id_nr = max(substructure_id_list)
    return substructure_id_nr

def parse_molblocks_data(molblocks_results, NPDatabase_path):
    """ 
    This function converts the substructures from molblocks into the 
    canonical RDkit SMILES. A dictionary (substructure_output_dict) will
    be returned that contains the structure id and all the canonical RDkit 
    substructures SMILES. Then the substuctures get checked whether
    they are already present in the substructure table (NPDatabase), if not, 
    then the substructure will be added to list for new substructures.
    
    
    molblocks_results: molblocks with substructures and structure id
    NPDatabase_path: path to NPDatabase    
    """
    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Check if substructure is already in substructure table
    # Extract all substructure_id's from sqlite database
    c.execute('SELECT * from substructure')
    subtructure_smiles_sqlite_list = []

    for row in c:
        subtructure_smiles_sqlite_list += [row[1]]

    rdkit_smiles_list = [] 
    substructure_output_dict = {}
    substructure_results = molblocks_results.split('\n')
    for subs in substructure_results[:-1]:
        subs = subs.split('\t')
        structure_id = subs[1]
        subs = subs[0].split('.')

        # To prevent warning print statements from RDkit 
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        
        for sub in subs: 
            if sub != '<NA>':
                mol = Chem.MolFromSmiles(sub)
                if mol is not None:
                    rdkit_smiles = Chem.MolToSmiles(mol)
                    if structure_id in substructure_output_dict.keys():
                        substructure_output_dict[structure_id] += [rdkit_smiles]
                    if structure_id not in substructure_output_dict.keys():
                        substructure_output_dict[structure_id] = [rdkit_smiles]
                    if rdkit_smiles not in subtructure_smiles_sqlite_list:
                        rdkit_smiles_list += [rdkit_smiles]
    # Just keep the unique substructures
    distinct_rdkit_smiles_list = set(rdkit_smiles_list)    

    return distinct_rdkit_smiles_list, substructure_output_dict


def create_substructure_data(parsed_data, substructure_id_nr):
    """ 
    This function generated and returns a list with all substructure data, including the subtructure id, mol formula,
    mol weigth, inchi and inchi key (all RDkit).
    
    
    parsed_data: list with all the distinct canonical rdkit substructure SMILES,
    which are not present in the substructure table 
    substructure_id_nr: the number of the highest substructure id  
    """
    all_RDkit_structure_data = []
    substructure_id_nr = int(substructure_id_nr)
    substructure_id_nr += 1

    
    for i, SMILES in enumerate(parsed_data):
       
        #Create unique identifier (substructure_id)              
        substructure_id = '000000{}'.format(substructure_id_nr)
        substructure_id = substructure_id[-8:]
        substructure_id = 'NPSB:{}'.format(substructure_id)
        temp_all_RDkit_structure_data = [substructure_id]        
        substructure_id_nr += 1
        #Canonical SMILES list
        temp_all_RDkit_structure_data += [SMILES]
        mol = Chem.MolFromSmiles(SMILES)

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
    
    
    print ('All RDkit attributes generated for ', len(all_RDkit_structure_data), ' new substructures' )
    return all_RDkit_structure_data
    
    
def store_substructure_table_data(substructure_table_data, NPDatabase_path):
    """
    This function stores all substructure data in the substructure table in
    NPDatabase
    
    substructure_table_data: list of lists from create_substructure_data()
    NPDatabase_path: path to NPDatabase
    """


    sql_list = [] 
    for line in substructure_table_data:
        sql_string = "INSERT INTO substructure VALUES ('%s', '%s', '%s', '%s',\
        '%s', '%s');"% (line[0], line[1], line[2], line[3], line[4], \
        line[5])
        sql_list += [sql_string] 
    

    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
     

    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
    
    
    print ('New substructures added to substructure table')
    
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()

def structure_has_substructure_table_data(molblocks_input, substructure_output_dict, NPDatabase_path):
    """ 
    This function creates the data for the structure_has_substructure table.
    
    molblocks_input: text file with structure SMILES and structure id's which
    is used as input for molBLOCKS
    substructure_output_dict: dict with canonical_smiles from the substructures and 
    the structure_id's they where generated from. Also the substructures which
    are already in the substructure table.
    NPDatabase_path: path to NPDatabase 
    """
  

    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Extract all smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure')
    can_smiles_sqlite_list = []
    structure_id_sqlite_list = []
    
    structure_dict = {}
    reverse_structure_dict = {}
    list_in_list = []
    
    for row in c:
        temp_list_in_list = []
        can_smiles = row[2]
        structure_id = row[0]
        temp_list_in_list = [can_smiles, structure_id]
        list_in_list += [temp_list_in_list]
                
    # Make a dict with the canonical smiles from NPDatabase with all structure ids (overlapping 
    # ones, since multiple isomeric smiles can have the same canonical smiles)
   
    smiles_all_structure_ids_dict = {}       
    structure_input = molblocks_input.split('\n')
    for line in structure_input[:-1]:
       line = line.split('\t')
       structure_id = line[1]
       smiles = line[0]
       for combo in list_in_list:
           NPID = combo[1]
           can_smiles = combo[0]
           if smiles == can_smiles:
               if smiles in smiles_all_structure_ids_dict.keys():
                  smiles_all_structure_ids_dict[smiles] += [NPID]
               if smiles not in smiles_all_structure_ids_dict.keys():
                  smiles_all_structure_ids_dict[smiles] = [NPID]

    
    # Extract all substructure smiles and substructure_id's from sqlite database
    c.execute('SELECT * from substructure')
    
    substructure_id_smiles_list = []   
    for row in c:
        temp_substructure_id_smiles_list = []
        structure_id_sqlite = row[0]
        can_smiles_sqlite = row[1]
        temp_substructure_id_smiles_list = [structure_id_sqlite, can_smiles_sqlite] 
        substructure_id_smiles_list += [temp_substructure_id_smiles_list]
        
    # Use output from molblocks, combine the structure ids and all substructure ids 
    # (replace the substructure SMILES with the corresponding substructure id from the sqlite db)
    sub_and_structure_id_dict = {}
    for key, values in substructure_output_dict.items():
        for value in values:
            smiles = value
            for combo in substructure_id_smiles_list:
                can_smiles = combo[1]
                substructure_id = combo[0]
                if smiles == can_smiles:
                    if key in sub_and_structure_id_dict.keys():
                       sub_and_structure_id_dict[key] += [substructure_id]
                    if key not in sub_and_structure_id_dict.keys():
                       sub_and_structure_id_dict[key] = [substructure_id]
    
    # Collect all structure id and substructure ids
    overall_list = []
    for values in smiles_all_structure_ids_dict.values():
        overall_list_temp = []
        for value in values:
            structure_id = value
            for key in sub_and_structure_id_dict.keys():
                if key == structure_id:
                    overall_list_temp = [values, sub_and_structure_id_dict[structure_id]]
                    overall_list += [overall_list_temp]
    
    # make structure id - substructure id pairs for all combinations
    structure_id_substructure_id_combo_list = []
    for line in overall_list:
        temp_structure_id_substructure_id_combo_list = []
        structure_id = line[0]
        substructure_id = line[1]
        for i in range(len(structure_id)):
            for j in range(len(substructure_id)):
                temp_structure_id_substructure_id_combo_list = [structure_id[i],substructure_id[j]]
                structure_id_substructure_id_combo_list += [temp_structure_id_substructure_id_combo_list]

    
    return structure_id_substructure_id_combo_list
    
def store_structure_has_substructure_table_data(combo_list, NPDatabase_path):
    """ 
    This function stores all structure id and substructure id combinations
    in the structure_has_substructure table in NPDatabase
    
    combo_list: list of lists with structure id-substructure id pairs
    NPDatabase_path: path to NPDatabase    
    """
    
    sql_list = [] 
    for line in combo_list:
        sql_string = "INSERT OR REPLACE INTO structure_has_substructure VALUES ('%s', '%s');"% (line[0], line[1])
        sql_list += [sql_string] 
    

    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
 

    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
    
    print ('Structure_has_substructure (junction) table updated') 
    
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()



if __name__ == "__main__":
    with open(argv[1]) as file_object1:
        molblocks_input = file_object1.read()

    with open(argv[2]) as file_object:
        molblocks_results = file_object.read()
        
        NPDatabase_path = "/path/to/NPDatabase.sqlite"
        
        substructure_id_nr = get_substructure_id(NPDatabase_path) 
        distinct_parsed_data, substructure_output_dict = parse_molblocks_data(molblocks_results, NPDatabase_path)
        
        # Create and store substructure table data
        substrucuture_table_data = create_substructure_data(distinct_parsed_data, substructure_id_nr)
        store_substructure_table_data(substrucuture_table_data, NPDatabase_path)
        
        # Create and store structure_has_substructure table data
        structure_id_substructure_id_combo_list = structure_has_substructure_table_data(molblocks_input, substructure_output_dict, NPDatabase_path) 
        store_structure_has_substructure_table_data(structure_id_substructure_id_combo_list, NPDatabase_path)
