
"""
Created on 19-04-'19
Script that processes molBLOCKS output, generates new substructure data,
creates the substructure table and adds the substructure data to it.

The new substructure table then is used for the generation of data for the
structure_has_substructure table. The latter is also created with this script
and the data will be stored.

Use RDkit version 2018.09.2

Usage:
Command line: python3 create_substructures.py [input__file_molBLOCKS] [output_file_molBLOCKS]
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

def create_substructure_data(molblocks_results):
    """ 
    This function converts the substructures into the substructures canonical
    rdkit SMILES. A dictionary (substructure_output_dict) will be returned that contains
    the structure id and all the canonical RDkit substructures SMILES. Also a list
    is returned with all substructure data, including the subtructure id, mol formula,
    mol weigth, inchi and inchi key (all RDkit).
    
    
    molblocks_results: output from molblocks with substructures and substructure id
    """
    
    
    rdkit_smiles_list = [] 
    substructure_output_dict = {}
    substructure_results = molblocks_results.split('\n')
    for subs in substructure_results[:-1]:
        subs = subs.split('\t')
        structure_id = subs[1]
        subs = subs[0].split('.')

        # to prevent warning print statements from RDkit 
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        

        for sub in subs:
             
            if sub != '<NA>':
                mol = Chem.MolFromSmiles(sub)
                rdkit_smiles = Chem.MolToSmiles(mol)
                rdkit_smiles_list += [rdkit_smiles]  # RDKit canonical SMILES for consistency
                if structure_id in substructure_output_dict.keys():
                    substructure_output_dict[structure_id] += [rdkit_smiles]
                if structure_id not in substructure_output_dict.keys():
                    substructure_output_dict[structure_id] = [rdkit_smiles]
    # Only keep the unique substructures
    distinct_rdkit_smiles_list = set(rdkit_smiles_list)
    
    all_RDkit_structure_data = []
    for i, SMILES in enumerate(distinct_rdkit_smiles_list):
       
        #Create unique identifier (substructure_id)              
        substructure_id = format(i+1, "07")
        substructure_id = 'NPSB:{}'.format(substructure_id)
        temp_all_RDkit_structure_data = [substructure_id]        

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

    return all_RDkit_structure_data, substructure_output_dict
    
    
def store_substructure_table_data(substructure_table_data, NPDatabase_path):
    """
    This function creates the substructure table and stores all substructure 
    data in the substructure table in NPDatabase
    
    structure_table_data: list of lists from create_substructure_data()
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
     
    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS substructure")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE substructure 
             (substructure_id text PRIMARY KEY, canonical_smiles text, \
             monoisotopic_mass num, molecular_formula text, \
             inchi text, inchi_key_rdkit text)''')
    
    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
    
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()

def structure_has_substructure_table_data(molblocks_input, substructure_output_dict, NPDatabase_path):
    """ 
    This function creates the data for the structure_has_substructure table.
     
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
                
    
    # make a dict with the canonical smiles from NPDatabase with all structure ids (overlapping 
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

    
    # Extract all substructure SMILES and substructure_id's from sqlite database
    c.execute('SELECT * from substructure')
    
    substructure_id_smiles_list = []
    
    for row in c:
        temp_substructure_id_smiles_list = []
        structure_id_sqlite = row[0]
        can_smiles_sqlite = row[1]
        temp_substructure_id_smiles_list = [structure_id_sqlite, can_smiles_sqlite] 
        substructure_id_smiles_list += [temp_substructure_id_smiles_list]
        
    # Use output from molblocks, combine the structure ids and all substructure ids 
    #(replace the substructure SMILES with the corresponding substructure id from the sqlite db)
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
    
    # collect all structure id and substructure ids
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
    This function creates the structure_has_substructure table and stores all 
    structure id and substructure id combinations in the structure_has_substructure 
    table in NPDatabase
    
    combo_list: list of lists with structure id-substructure id pairs
    NPDatabase_path: path to NPDatabase    
    """
    
    sql_list = [] 
    for line in combo_list:
        sql_string = "INSERT INTO structure_has_substructure VALUES ('%s', '%s');"% (line[0], line[1])
        sql_list += [sql_string] 
    

    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
 

    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS structure_has_substructure")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE structure_has_substructure 
             (structure_id text, substructure_id text, \
             PRIMARY KEY (structure_id, substructure_id),
             FOREIGN KEY(structure_id) REFERENCES structure(structure_id),
             FOREIGN KEY(substructure_id) REFERENCES substructure(substructure_id))''')

    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
    
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()



if __name__ == "__main__":
    with open(argv[1]) as file_object1:
        molblocks_input = file_object1.read()

    with open(argv[2]) as file_object:
        molblocks_results = file_object.read()
        substrucuture_table_data, substructure_output_dict = create_substructure_data(molblocks_results)
        
        NPDatabase_path = "/mnt/nexenta/stokm006/NPDatabase_05-4.sqlite"
        
      ##  store_substructure_table_data(substrucuture_table_data, NPDatabase_path)
        structure_id_substructure_id_combo_list = structure_has_substructure_table_data(molblocks_input, substructure_output_dict, NPDatabase_path)
         
        store_structure_has_substructure_table_data(structure_id_substructure_id_combo_list, NPDatabase_path)
