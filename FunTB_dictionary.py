import pandas as pd
import numpy as np
import os
import sys

# Define directories
MTBSEQ_DIR = "MTBSeq_files"
METADATA_DIR = "Metadata_files"
SAVE_DIR = "Variations_dictionaries"

# Ensure directories exist
os.makedirs(SAVE_DIR, exist_ok=True)

# Get input files
mtbseq_filename = sys.argv[1]  # MTBseq file (.tab)
metadata_filename = sys.argv[2]  # Metadata file (.csv)

# Build full paths
mtbseq_filepath = os.path.join(MTBSEQ_DIR, mtbseq_filename)
metadata_filepath = os.path.join(METADATA_DIR, metadata_filename)

# Read input files
MTBseq_file = pd.read_csv(mtbseq_filepath, sep='\t', skiprows=1)
Metadata_file = pd.read_csv(metadata_filepath)

# Dictionary name based on MTBSeq file
dict_name = os.path.splitext(mtbseq_filename)[0]  # Remove file extension

# Reference sequences and gene names
ReferencSequence = MTBseq_file['Ref'].values
GeneNames = MTBseq_file['GeneName'].values

# Amino acid dictionary
amino_acid_dict = {
    'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',
    'glu': 'E', 'gln': 'Q', 'gly': 'G', 'his': 'H', 'ile': 'I',
    'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',
    'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V'
}

# Function to extract variation parameters
def GetVariationsParameters(variation):
    if variation[0] == variation[-1] == '_': 
        return None
    elif variation[0] == '_':  
        position = variation[1:-3]
        traduction = '_' + variation[1:-3] + amino_acid_dict.get(variation[-3:], '?')
        return position, traduction
    elif variation[-1] == '_':  
        position = variation[3:-1]
        traduction = amino_acid_dict.get(variation[:3], '?') + variation[3:-1] + '_'
        return position, traduction
    elif variation[:3] != variation[-3:]:  
        position = variation[3:-3]
        traduction = amino_acid_dict.get(variation[:3], '?') + variation[3:-3] + amino_acid_dict.get(variation[-3:], '?')
        return position, traduction
    return None

# Function to filter relevant data
def GetRelevantData(variations, genes=GeneNames):
    variations = np.array(variations)
    indxs = np.where((variations != '-') & (variations != ' ') & pd.notnull(variations) & (variations.astype(str) != 'nan'))
    return variations[indxs], GeneNames[indxs]

# Function to create the variation dictionary
def CreateDictionary():
    variations_dict = {}

    for sample in MTBseq_file.columns[-len(Metadata_file.index):]:
        sample_name = sample.split('.')[0]
        variations_dict[sample_name] = {}
        mutations, genes = GetRelevantData(MTBseq_file[sample].values)

        for gene, mutation in zip(genes, mutations):
            variation = mutation.split(' ')[0].lower()
            variation_params = GetVariationsParameters(variation)

            if gene != '-' and variation_params:
                if gene not in variations_dict[sample_name]:
                    variations_dict[sample_name][gene] = {
                        'total_variations': 1,
                        'variation_positions': {variation_params[0]: 1},
                        'symbolic_mutations': {variation_params[1]: 1}
                    }
                else:
                    variations_dict[sample_name][gene]['total_variations'] += 1
                    variations_dict[sample_name][gene]['variation_positions'][variation_params[0]] = \
                        variations_dict[sample_name][gene]['variation_positions'].get(variation_params[0], 0) + 1
                    variations_dict[sample_name][gene]['symbolic_mutations'][variation_params[1]] = \
                        variations_dict[sample_name][gene]['symbolic_mutations'].get(variation_params[1], 0) + 1

    return variations_dict

# Function to save the dictionary
def SaveDictionaryFile(variation_dict):
    file_path = os.path.join(SAVE_DIR, f"{dict_name}.txt")
    try:
        with open(file_path, 'w') as file:
            file.write(str(variation_dict))
        print(f"Dictionary saved as {dict_name}.txt in {SAVE_DIR}")
    except Exception as e:
        print(f"Error saving dictionary: {e}")

# Generate and save the dictionary
variation_dict = CreateDictionary()
SaveDictionaryFile(variation_dict)