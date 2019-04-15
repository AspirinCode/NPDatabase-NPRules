"""
Script 1 for creating NPDatabase
This script takes the data from the external database file and writes a
file with only the SMILES (directly copied from source, thus  not
RDkit isomerical/ canonical), external_db_id and the db name
Note that the input file name (minus .txt) will be used as external
database name. 

This script has to be used for each database file you want to add, in order
to create one large datafile with all database data from several databases.


Usage:
Command line: python3 collect_all_structures.py [input_external_database] [name_outputfile]

"""
import sys, getopt
import os.path


def create_all_input_data(input_file, output_file):
    """  
    This function takes the data from the external database file and writes a
    file with only the SMILES, external_db_id and the db name
    
    input_file: file with input data (CLASS format)
    output_file: name of output file, which collects all structure data 
    from all external databases
    
    """


    db_name = input_file[:-4]
    with open(input_file) as file_object:
        input_file = file_object.read()  
   
   
    SMILES_list = []
    identifier_list = []
    all_lines = input_file.split('\n')
    for line in all_lines[1:-1]:
        line = line.split('\t')
        if line[0] != '':  # for NuBBE db
            SMILES = line[2].strip('"')
            if not SMILES.startswith('no dat'): # for GNPS db
                SMILES_list += [SMILES]
                external_db_identifier = line[3].strip('"')
                identifier_list += [external_db_identifier]

    if not os.path.isfile(output_file):
        with open(output_file, 'w') as db_file:
            db_file.write('SMILES' + '\t' + 'External_db_identifier'+ '\t' + 'External_db_name' + '\n') 
        
    with open(output_file, 'a') as db_file:
        for i in range(len(SMILES_list)):
            db_file.write(SMILES_list[i] + '\t' + identifier_list[i] + '\t' + db_name+ '\n')
    


def main(argv):
    """
    Main function that calls all other functions 
    
    """
    
    input_file_name = argv[0]
    output_file_name = argv[1]
    create_all_input_data(input_file_name, output_file_name)


if __name__ == '__main__':
    
   main(sys.argv[1:])
