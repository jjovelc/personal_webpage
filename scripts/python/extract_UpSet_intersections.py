#!/usr/bin/python

import pandas as pd
from itertools import combinations, chain
import glob

# Function to read DESeq2 file and return upregulated and downregulated transcripts
def read_deseq2_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    upregulated = df[df['log2FoldChange'] > 0]['transcript'].unique()
    downregulated = df[df['log2FoldChange'] < 0]['transcript'].unique()
    return upregulated, downregulated

# Function to find exclusive intersections
def find_exclusive_intersections(sets_list):
    all_intersections = {}
    exclusive_intersections = {}
    # First find all intersections
    for i in range(len(sets_list), 0, -1):  # Start from the largest intersection
        for subset in combinations(sets_list, i):
            key = '_AND_'.join(sorted([s[1] for s in subset]))
            intersection_set = set.intersection(*[s[0] for s in subset])
            all_intersections[key] = intersection_set

    # Then find exclusive sets by subtracting larger sets from smaller ones
    for smaller_key in all_intersections:
        exclusive_set = all_intersections[smaller_key].copy()
        for larger_key in all_intersections:
            if len(larger_key.split('_AND_')) > len(smaller_key.split('_AND_')):
                exclusive_set -= all_intersections[larger_key]
        exclusive_intersections[smaller_key] = exclusive_set

    return exclusive_intersections

# Directory where your DESeq2 files are stored
directory_path = '/Users/juanjovel/jj/data_analysis/coltonUnger/DE_analysis/mm10_plus_RNAspades_assembly/UpSet'

# Automatically find all files ending in *_q0.05.tsv
files = glob.glob(f"{directory_path}/*_q0.05.tsv")

# Read files and prepare lists for intersections
upregulated_sets = []
downregulated_sets = []

for file_path in files:
    up, down = read_deseq2_file(file_path)
    # Extract filename for labeling purposes
    filename = file_path.split('/')[-1].replace('_q0.05.tsv', '')
    upregulated_sets.append((set(up), filename))
    downregulated_sets.append((set(down), filename))

# Calculate exclusive intersections
upregulated_intersections = find_exclusive_intersections(upregulated_sets)
downregulated_intersections = find_exclusive_intersections(downregulated_sets)

# Function to save intersections to a file
def save_transcripts_to_files(intersections, prefix):
    for key, transcripts in intersections.items():
        # Format the filename based on the key
        filename = key + prefix + '.txt'
        # Write the transcripts to the file
        with open(filename, 'w') as f:
            for transcript in transcripts:
                f.write(f"{transcript}\n")

# Call the function for upregulated and downregulated intersections
save_transcripts_to_files(upregulated_intersections, '_up')
save_transcripts_to_files(downregulated_intersections, '_down')
