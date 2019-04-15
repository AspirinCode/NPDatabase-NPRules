"""

Adds the classified temp structure table to NPDatabase, also the
structure_has_data_source and the data_source table will be updated

Useage:
Command line: python3 update_NPDatabase.py [parsed new input file] [path/to/NPDatabase.sqlite] [path/to/temp_NPDatabase.sqlite]

Use RDkit version 2018.09.2

"""
  
  
  
from rdkit import Chem
from sys import argv
import sqlite3


def update_structure_table(NPDatabase_temp_path, NPDatabase_path):
    """
    Takes all data from the classified temp structure table and adds
    the data to NPDatabase (in total 23 attribiutes/ columns per structure) 
    
    NPDatabase_temp_path: path to temp NPDatabase
    NPDatabase_path: path to NPDatabase   
    """
  
   
    conn = sqlite3.connect(NPDatabase_temp_path)
    c = conn.cursor()
    
    # Retrieve all smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure_temp')
    
    temp_list = []
    
    row_count = 0
    sql_list = []
    for row in c:
        structure_id = row[0]
        isomeric_smiles = row[1]
        canonical_smiles = row[2]
        monoisotopic_mass = row[3]
        molecular_formula = row[4]
        inchi = row[5].replace("'","")
        inchi_key_rdkit = row[6]
        inchi_key_molconvert = row[7]
        cf_direct_parent = row[8].replace("'","")
        cf_kingdom = row[9].replace("'","")
        cf_superclass = row[10].replace("'","")
        cf_class = row[11].replace("'","")
        cf_subclass = row[12].replace("'","")
        cf_intermediate_0 = row[13].replace("'","")
        cf_intermediate_1 = row[14].replace("'","")
        cf_intermediate_2 = row[15].replace("'","")
        cf_intermediate_3 = row[16].replace("'","")
        cf_intermediate_4 = row[17].replace("'","")
        cf_intermediate_5 = row[18].replace("'","")
        cf_molecular_framework = row[19].replace("'","")
        cf_alternative_parents = row[20].replace("'","")
        cf_substituents = row[21].replace("'","")
        cf_description = row[22].replace("'","")
        

        row_count += 1
        sql_string = "INSERT INTO structure (structure_id, isomeric_smiles, \
        canonical_smiles, monoisotopic_mass, molecular_formula, inchi,\
        inchi_key_rdkit, inchi_key_molconvert, cf_direct_parent, cf_kingdom,\
        cf_superclass, cf_class, cf_subclass, cf_intermediate_0, \
        cf_intermediate_1, cf_intermediate_2, cf_intermediate_3, \
        cf_intermediate_4, cf_intermediate_5, cf_molecular_framework,  \
        cf_alternative_parents, cf_substituents, cf_description) VALUES ('%s', \
        '%s', '%s', '%s','%s','%s','%s','%s','%s', '%s', '%s', '%s','%s','%s',\
        '%s','%s', '%s','%s','%s','%s','%s','%s','%s');"%(structure_id, isomeric_smiles, \
        canonical_smiles, monoisotopic_mass, molecular_formula, inchi,\
        inchi_key_rdkit, inchi_key_molconvert, cf_direct_parent, cf_kingdom,\
        cf_superclass, cf_class, cf_subclass, cf_intermediate_0, \
        cf_intermediate_1, cf_intermediate_2, cf_intermediate_3, \
        cf_intermediate_4, cf_intermediate_5, cf_molecular_framework,  \
        cf_alternative_parents, cf_substituents, cf_description)

        sql_list += [sql_string]


    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
    print (len(sql_list), " new structures added to structure table")
    # Commit your changes
    conn.commit()


def update_structure_has_source_table(input_file, NPDatabase_path, new_source_name):
    """
    Updates the structure_has_source table. This function first checkes 
    whether the external source id is already in the table, if not, a new 
    structure_source_id is created and the new data is added.
    
    
    input_file: txt file, new parsed input data with external source id and isomeric smiles
    NPDatabase_path: path to NPDatabase   
    new_source_name: name of database where the new data came from 
    """

    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Retrieve all isomeric smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure')

    sqlite_combo_list = []
    for row in c:
        structure_id = row[0]
        iso_smiles = row[1]
        temp_combo_list = [structure_id, iso_smiles]
        sqlite_combo_list += [temp_combo_list]
        temp_combo_list = []

    # Adds the correct structure id to the isomeric smiles from the parsed input data
    structure_has_source_list = []
    all_lines = input_file.split('\n')
    for line in all_lines[1:-1]:
        line = line.split('|||')
        iso_smiles_mibig = line[1]
        external_id = line[0]
        for combo in sqlite_combo_list:
            sqlite_iso_smiles = combo[1]
            structure_id = combo[0]
            if iso_smiles_mibig == sqlite_iso_smiles:
                temp_structure_has_source_list = [structure_id, iso_smiles_mibig]
                structure_has_source_list += [temp_structure_has_source_list]
                temp_structure_has_source_list = []
       
    # Retrieve all isomeric smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure_has_data_source')
    structure_source_id_sqlite_list = []
    external_source_id_sqlite_list = []
    for row in c:
        structure_source_id_sqlite_list += [row[0]]
        external_source_id_sqlite_list += [row[3]]
    
    # For structure_id's: delete "NPSRC:" and only keep the nr
    structure_source_id_list = []
    for structure_id in structure_source_id_sqlite_list[1:]:
        structure_source_id_list += [structure_id[7:]]

    # Get the highest structure_source_id nr    
    structure_source_id_nr = max(structure_source_id_list)
  
    
    external_source_name = new_source_name
    structure_source_count = 0
    sql_list = []
  
    # Check whether external source id from parsed input is already in structure_has_data_source table and created a sql list
    for combo in structure_has_source_list:

        structure_id = combo[0]
        external_source_id = combo[1]

        if external_source_id in external_source_id_sqlite_list:
            print (external_source_id, '   already in database')
        if external_source_id not in external_source_id_sqlite_list:
            structure_source_count += 1
            structure_source_id_nr = int(structure_source_id_nr)+1
            new_structure_source_id = 'NPSRC:00{}'.format(structure_source_id_nr)
            sql_string = "INSERT INTO structure_has_data_source VALUES\
            ('%s','%s', '%s', '%s');"% (new_structure_source_id, structure_id, external_source_id, external_source_name)
            sql_list += [sql_string] 
            structure_source_id_nr = structure_source_id_nr+1
  
    source_name_count_list = [external_source_name, structure_source_count]

    # Add new structure source data to structure_has_data_source table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])

    # Commit changes
    conn.commit()

    # Close the connection
    conn.close()
    
    return source_name_count_list 
    
def update_source_table(source_name_count_list, NPDatabase_path):
    """ adds the new database name and nr of structures or updates the nr of 
    stuctures in the data_source table.
        
    source_name_count_list: source name and number of new structures
    NPDatabase_path: path to NPDatabase  
    """

    
    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Retrieve all smiles and structure_id's from sqlite database
    c.execute('SELECT * from data_source')

    source_name_sqlite_list = []
    structure_number_sqlite_list = []
    

    source_number_dict = {}
    for row in c:
        source_name = row[0]
        nr_of_structures = int(row[1])
        source_name_sqlite_list += [source_name]
        structure_number_sqlite_list += [nr_of_structures]
        source_number_dict[source_name] = nr_of_structures
    source_name = source_name_count_list[0]
    structures_in_source = source_name_count_list[1]
    
    
    # if new source name not in data_source table, it will be added together with the number of structures added
    sql_list = []
    not_in_sqlite_count = 0
    if source_name not in source_name_sqlite_list:
         not_in_sqlite_count += 1
         sql_string = "INSERT INTO data_source VALUES\
            ('%s','%s');"% (source_name, structures_in_source)
         sql_list += [sql_string] 
         print (source_name, "and number of structures (", structures_in_source, ") added to data_source table") 
    
    # if new source name already in data_source table, the number of structure will be updated
    delete_list = []
    if source_name in source_name_sqlite_list and structures_in_source != 0:
         current_number = int(source_number_dict[source_name])
         new_number = int(structures_in_source) + current_number

         delete_string = "DELETE from data_source where source_name = ('%s');"% (source_name)
         delete_list +=  [delete_string]
         sql_string = "INSERT INTO data_source VALUES ('%s', %s);"% (source_name, new_number)
         sql_list += [sql_string]
         print (source_name, "already in data_source tabel", structures_in_source, 'structures added')

    
   # Add new structure source data to structure source table
    for i in range(len(sql_list)):
        if len(delete_list) != 0:
             c.execute(delete_list[i])
        c.execute(sql_list[i])


    # Commit changes
    conn.commit()

    # Close the connection
    conn.close()
    


if __name__ == "__main__":

    new_source_name = argv[1]
    new_source_name = source_name[:-11]

    NPDatabase_path = argv[2]
    NPDatabase_temp_path = argv[3]
    
    
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        update_structure_table(NPDatabase_temp_path, NPDatabase_path)

        source_name_count_list = update_structure_has_source_table(input_file, NPDatabase_path, new_source_name)
        update_source_table(source_name_count_list, NPDatabase_path)
