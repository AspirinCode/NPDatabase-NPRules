
"""
Creates structure data; new structure_id, canonincal_smiles (RDkit), 
monoisotopic_mass (RDkit), mol_formula (RDkit), inchi (RDkit),
inchi_key (RDkit) and ichi_key (molconvert). This data will be stored
in a temp structure table in a specified temp NPDatabse file. This temp
Database file is suitable for classification with
with https://github.com/OscarHoekstra/ClassifyNPDB 

Use RDkit version 2018.09.2

Usage:
Command line: python3 create_temp_structure_table.py [molconvert output] [parsed new input file] [path/to/NPDatabase.sqlite] [path/to/temp_NPDatabase.sqlite]
[path/to/temp_NPDatabase.sqlite] will become a new sqlite file

Example:
Command line: python3 create_temp_structure_table.py molconvert_output.txt new_input_parsed.txt /mnt/nexenta/stokm006/NPDatabase.sqlite /mnt/nexenta/stokm006/NPDatabase_temp.sqlite
"""
 
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdinchi
from sys import argv
import sqlite3



def get_structure_id(NPDatabase_path):
    """ 
    function that extracts and returns the structure id with the highest
    number from NPDatabase
    
    NPDatabase_path: path to NPDatabase    
    """
    # Connect to sqlite database
    conn = sqlite3.connect(NPDatabase_path)
    c = conn.cursor()
    
    # Retrieve all smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure')
    iso_smiles_sqlite_list = []
    structure_id_sqlite_list = []
    for row in c:
        structure_id_sqlite_list += [row[0]]

   # For structure_id's: delete "NPDB:" and only keep the nr
    structure_id_list = []
    for structure_id in structure_id_sqlite_list:
       structure_id_list += [int(structure_id[5:])]
    
    # Get the highest structure_id nr    
    structure_id_nr = max(structure_id_list)
    structure_id_nr = structure_id_nr+1
 
    return structure_id_nr
    
def generate_structure_data(input_file_molconvert, input_file, structure_id_nr):
    """ takes all text from the input structure data file and returns a list of
    lists with all generated data needed for the sqlite database. The isomeric smiles 
    from the parsed data is required, since the smiles from the molconvert output can
    be altered (neutralized).
  
    input_file_molconvert: txt file, with output from molconvert pipeline, 
    containing structure id, isomeric smiles and molconvert inchi key  
    input_file: txt file, new parsed input data with external source id and isomeric smiles
    structure_id_nr: int, from get_structrure_id()
    
    """
    all_input_combo_list = []
    all_lines1 = input_file.split('\n')
    for line in all_lines1[1:-1]:
        line = line.split("|||")
        external_source_id = line[0]
        iso_smiles = line[1]
        temp_all_input_combo_list = [[external_source_id], [iso_smiles]]
        all_input_combo_list += [temp_all_input_combo_list]
        temp_all_input_combo_list = []

    external_source_id = []
    overall_new_structure_list1 = []
    all_lines = input_file_molconvert.split('\n')
    for line in all_lines[:-1]:

        overall_temp_list1 = []
        new_structure_list = []
        line = line.split(' ')
        external_source_id = line[0]
        inchi_key = line[2]
        for combo in all_input_combo_list:
            external_source_id_mibig = combo[0][0]
            isomeric_smiles = combo[1][0]

            if external_source_id_mibig == external_source_id:
                overall_temp_list1 = [external_source_id, isomeric_smiles, inchi_key]
                overall_new_structure_list1 += [overall_temp_list1]


    sorted_overall_list1 = sorted(overall_new_structure_list1)
   
    overall_new_structure_list = []
    for i, line in enumerate(sorted_overall_list1):

        new_structure_list = []        
        external_source_id = line[0]
        isomeric_smiles = line[1]
        inchi_key_molconvert = line[2]
        inchi_key_molconvert = inchi_key_molconvert[9:]
        
        #New structure id
        structure_id = 'NPDB:0{}'.format(structure_id_nr)
        structure_id_nr += 1
        new_structure_list += [structure_id]
        #Generate rdkit mol
        mol = Chem.MolFromSmiles(isomeric_smiles)
        #SMILES, isomeric
        new_structure_list += [isomeric_smiles] 
        #SMILES, canonical
        can_smiles = Chem.MolToSmiles(mol, isomericSmiles = False)
        new_structure_list += [can_smiles] 
        #Monoisotopic mass
        mol_weigth = Descriptors.ExactMolWt(mol)
        new_structure_list += [mol_weigth] 
        #Mol Forumula
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        new_structure_list += [mol_formula]
        #InChI         
        inchi = Chem.MolToInchi(mol)
        new_structure_list += [inchi]
        #InChIKey RDkit
        inchi_key = Chem.MolToInchiKey(mol)     
        new_structure_list += [inchi_key]
        #InChiKey molconvert
        new_structure_list += [inchi_key_molconvert]
        overall_new_structure_list += [new_structure_list]


    return (overall_new_structure_list)

def create_temp_structure_table(all_structure_data, NPDatabase_temp_path):
    """
    Creates a temporary structure table which can be classified easily
    with https://github.com/OscarHoekstra/ClassifyNPDB
    
    all_structure_data: list of lists with new structure data, generated
    with RDkit and the molconvert inchi key.
    NPDatabase_temp_path: path to temp NPDatabase
    """
         

    sql_list = [] 
    for line in all_structure_data:
        sql_string = "INSERT INTO structure_temp VALUES ('%s', '%s', '%s', '%s', '%s',\
        '%s', '%s', '%s');"% (line[0], line[1], line[2], line[3], line[4], \
        line[5], line[6], line[7])
        sql_list += [sql_string]
     

    
    conn = sqlite3.connect(NPDatabase_temp_path) # path where temp structure table will be stored
    c = conn.cursor()
     
     # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS structure_temp")
    
     # Add table with atrribute names
    c.execute('''CREATE TABLE structure_temp 
              (structure_id text PRIMARY KEY, isomeric_smiles text, canonical_smiles text, \
              monoisotopic_mass num, molecular_formula text, \
              inchi text, inchi_key_rdkit text, inchi_key_molconvert text)''')

     # Add data to table
    for i in range(len(sql_list)):
       c.execute(sql_list[i])
  
     # Commit your changes
    conn.commit()

        
if __name__ == "__main__":
    
    NPDatabase_path = argv[3]
    NPDatabase_temp_path = argv[4]

    with open(argv[1]) as file_object:
        input_file_molconvert = file_object.read()
    with open(argv[2]) as file_object2:
        new_input_structures = file_object2.read()
       
        structure_id_nr = get_structure_id(NPDatabase_path)
        all_structure_data = generate_structure_data(input_file_molconvert, new_input_structures, structure_id_nr)
        create_temp_structure_table(all_structure_data, NPDatabase_temp_path)
