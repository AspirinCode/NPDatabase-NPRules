# NPDatabase-NPRules


Use RDkit version 2018.09.2 (install with Conda: https://anaconda.org/rdkit/rdkit)  
The command line to use is specified in each script

You can download the latest version of NPDatabase (the sqlite file) here: https://www.dropbox.com/s/hk5dnjq262zlks3/NPDatabase_05-4.sqlite?dl=0   


## Build NPDatabase
_all scripts in Build_NPDatabase directory_  

The instructions and scripts to build NPDatabase can be found in Build NPDatabase dir:  
https://github.com/SamStokman/NPDatabase-NPRules/tree/master/Build_NPDatabase

## Add new structure data to NPDatabase
_All scripts in Add_new_structure directory_  

__Step 1:__  
Input: New structure data  
Script: parse_new_input.py  
Output: text file with external_source_id and the isomeric/ canoninal SMILES  

Since each database with structure data uses different formats to store its data, first the data must be parsed. The minimal input requirements are the structures SMILES, the external source identifier and the name of the database where the new structure data comes from. The script, parse_new_input.py, can be used which has to be customized according to the input. The output of this script contains the external_source_id and the isomeric/ canoninal SMILES generated with RDKit.

__Step 2:__  
Input: text file with external_source_id and the isomeric/ canoninal SMILES  
Script: create_molconvert_input.py  
Output: csv file with external source id and SMILES  

The output from step 1, then can be used to check whether the structure is already present in NPDatabase. Structures which are not yet present in NPDatabase will be written into a csv file suitable for the molconvert pipeline.

__Step 3:__  
Input: csv file with external source id and SMILES
Ouput: several text files with external source id's, SMILES and inchi keys 

Create molconvert inchi keys. First, install the molconvert pipeline, created by Rutger Ozinga: https://github.com/NP-Plug-and-Play-Scripts/inchiKeyCreatorPipeline. You can find the instructions to install the pipeline in InchiKey pipeline manual.pdf.  
Command to use the pipeline:  
```
python3 [path/to/workplace/InchiKeyCreatorPipeline.py] [path/to/workplace/JChem/jchemsuite/bin/] [/path/to/workplace/input_SMILES_files] [inchi_key_input.csv]  
```
This command line contains the path to the script to run the pipeline, the path to molconvert program, the path to the folder that contains the files with the input SMILES and the name of the input file.

__Step 4:__  
Input: several text files with external source id's, SMILES and inchi keys  
Output: text file with external source id's, SMILES and inchi keys  

The pipeline outputs the inchi keys in separated files, the files that contain the inchi key data end with ‘_part_XX_dataFile.txt’ . To concatenate these files into one file use the cat command:
```
cat [file1.txt] [file2.txt] [file3.txt] > [all_output_collected.txt]
```

__Step 5:__  
Input: text file with external_source_id and the isomeric/ canoninal SMILES (from step 1)  
Input: text file with external source id's, SMILES and inchi keys  
Script: create_temp_structure_table.py  
Output: SQLite file with temp structure table  

Create a temporary structure table with all structure data from RDkit and the molconvert inchi keys. 
Attributes: structure_id, isomeric_smiles, canonical_smiles, monoisotopic_mass, molecular_formula, 
inchi, inchi_key_rdkit, inchi_key_molconvert

__Step 6:__  
Input: SQLite file with temp structure table  
Output: SQLite file with temp structure table with classified structures  

Classify the structures in the temporary structure table (in temporary NPDatabase) with the Classyfire pipeline, created by Oscar Hoekstra. First install the Classifire pipeline: https://github.com/OscarHoekstra/ClassifyNPDB.  

In the config file, adjust (also add a  “, “ at the end of line 83):  
-	Workbase  
-	Path to temporary NPDatabase   

Then, in ClassifyNPDB directory, use the command line:  
```
python3 ClassifyNPDB.py  
```
The temporary structure table is now fully classified and also contains attributes:   
cf_direct_parent, cf_kingdom, cf_superclass, cf_class, cf_subclass, cf_intermediate_0, cf_intermediate_1, cf_intermediate_2, cf_intermediate_3, cf_intermediate_4, cf_intermediate_5, cf_molecular_framework,  cf_alternative_parents, cf_substituents, cf_description

__Step 7:__  
Input: text file with external_source_id and the isomeric/ canoninal SMILES (from step 1)   
Input: SQLite file with temp structure table  
Input: NPDatabase SQLite file  
Output: Updated NPDatabase SQLite file  

Script to use: Update_NPDatabase.py
The data for the temporary structure table will be added to NPDatabase. The structure_has_data_source and data_source table will be updated accordingly


## Generate and analyse substructures
_All scripts in Use_molBLOCKS directory_  

__Install molBLOCKS__  

First download the modified version of molBLOCKS and upload it into the working directory of choice.
https://github.com/kheikamp/modified_molBLOCKS

molBLOCKS requires Openbabel, if you have root access, you just can install Openbabel as described in the original molBLOCKS user guide: http://compbio.cs.princeton.edu/molblocks/molblocksguide.pdf

If you do not have root access, you can install Openbabel as described in the link below:
https://openbabel.org/wiki/Install_(source_code)#Installing_locally_without_root_access

When Openbabel is installed successfully, the final step is to compilate the molBLOCKS suite. Enter the modified_molBLOCKS directory, and open Makefile. 

Adjust these two paths, usr should be replaced by the path to Openbabel directory:  
INCLUDES := -Iboost -I/usr/include/openbabel-2.0  
LIBS=-L/usr/local/lib -lm -lopenbabel  

Close and save the changes and then type and enter ‘make’. The installation should be finished now.

__Prepare input data__  
Input: NPDatabase SQLite file
Script: create_molBLOCKS_input.py
Output: text file with SMILES and structure id's 

First, you have to create input suitable for molBLOCKS. This script extracts structures from NPDatabase based on the taxonomy type (superclass/ class/ subclass) and name. The stucture_id and SMILES are written into a file suitable for molBLOCKS. The NPDatabase path should be specified in the script.

Store the output file in a new subdirectory in the modified_molBLOCKS directory.

__Generate substructures with molBLOCKS__  
Input: text file with SMILES and structure id's 
For further input, see parameter description
Output: text file with substructure SMILES and structure id's  

Rulesets: 
-	NPRules.txt
-	AmideRules.txt
-	LPRules.txt
-	BNRules.txt
-	OHRules.txt

Download the rulesets and upload them into the modified_molBLOCKS directory. Go into the directory where you stored the output file from create_molBLOCKS_input.py. From here you can use the command line:	
```
../fragment -r ../NPRules.txt -i input_structure_file.txt -n 2 -m 100 -k 1 -w 1000 -s 0.99 -o output_structure_file.txt
```
Parameter description:  
-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Input file with structures  
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;File with ruleset  
-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Maximum molecular weight  
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Minimum number of atoms  
-m&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Maximum number of atoms  
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fragment size relative to the parent structure  
-k&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of fragments that should be connected and considered as new fragment  
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Output file  

__Quantity results__  
Input: text file with substructure SMILES and structure id's (ouput from molBLOCKS) 
Script: analyze_molBLOCKS_output.py   

This script will calculate and output:  
-	Number of generated substructures per input structure
-	Total number of input structures
-	Total number of substructures
-	Total number of unique substructures
-	Percentage of fragmented input structures

__Quality results, select top 10__   
Input: output_from_molBLOCKS.txt  
Script: select_top_10.py   
Output: csv file with top 10 results

This script will select the top 10 most occurring substructures found in de molBLOCKS output file. The results are written into a csv file and the substructures will be shown in a pop-up display. This scripts requires ImageMagick, see: https://anaconda.org/conda-forge/imagemagick.

__Quality results, show substructures__   
Input: output_from_molBLOCKS.txt  
Script: Show_substructures.py  

This script is shows the generated substructures in a pop-up display, this script also requires ImageMagick. If you choose ‘y’ as method to show the substructure, keep in mind that the substructure shown is not necessarily has the correct position relative to the structure, it is just a substructure match. 
