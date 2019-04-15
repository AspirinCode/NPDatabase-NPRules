
"""

This script is just an example to parse new structure data from a new
data source. The minimum input for this script requires the structure's
SMILES and external source id.

The external source id has to be unique. 

Input: all structures from MIBiG: https://github.com/OscarHoekstra/ClassifyNPDB/blob/master/InFiles/All_MIBiG_compounds_with_SMILES_and_PMID_MAS.txt
Script that checks for structures with missing smiles and double accession numbers etc.
The external source id's for the MIBiG data were not unique and therefore, the 
compound name is added to it.

Result a text file with unique source_ids (accession number + compound name) and the 
isomeric/ canononical RDkit smiles.

Use RDkit version 2018.09.2

The name of the input file (excluding ".txt") will be used as Database name,
if this Database name is already in NPDatabas, copy that exact name.

Usage:
Command line: python3 parse_new_input.py [file with new structures]

"""

from rdkit import Chem
from sys import argv

def parse_data(input_file, database_name):
    """
    This scripts parses structure data and result a text file with 
    unique source_ids and the isomeric/ canononical RDkit smiles.
    
    input_file: txt file with new data
    database_name: name of database where the new data comes from
    """

    output_file_name = database_name + '_parsed.txt'
    all_lines = input_file.split('\n')

    with open(output_file_name, 'w') as db_file:
            db_file.write('Exteral_source_id_and_name'+ '\t' + 'Isomeric_RDkit_smiles' + '\n')
    
    iso_smiles_list = []
    for line in all_lines[1:-1]:
        line = line.split('\t')

        source_id = line[0]
        smiles = line[2]
        name = line[1]
        
        name = name.replace(",","-")   # , and ' and space characters must be excluded 
        name = name.replace("'","")
        name = name.replace(' ', '_')

    
        if smiles != '':   # the MIBiG data contains structures without SMILES strings
            mol = Chem.MolFromSmiles(smiles)
            iso_smiles = Chem.MolToSmiles(mol)   

            if iso_smiles not in iso_smiles_list:
               iso_smiles_list += [iso_smiles]

               with open(output_file_name, 'a') as db_file:
                  db_file.write(source_id + '-' + name + '|||' + iso_smiles + '\n')
    

if __name__ == "__main__":
    
    db_name = argv[1]
    db_name = db_name[:-4]
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        parse_data(input_file, db_name)
