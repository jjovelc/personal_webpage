import pandas as pd
import sys
import re

"""
USAGE: python normalize_kraken_counts.py <input_file.tsv>

Output file will have the same prefix with suffix CPM, 
e.g. <input_file_CPM.tsv>
"""


# Read the input table with taxa names and count data
infile = sys.argv[1]
input_table = pd.read_csv(infile, sep='\t')

# Extract the count data columns (excluding the taxa names column)
count_data = input_table.iloc[:, 1:]

# Calculate the sum of counts for each column (sample)
sum_counts = count_data.sum()

# Normalize the count data to counts per million (CPM)
normalized_data = (count_data / sum_counts) * 1000000

# Add the taxa names column back to the normalized data
normalized_table = pd.concat([input_table.iloc[:, 0], normalized_data], axis=1)

# Export results to a tsv file
outfile = re.sub('.tsv','_CPM.tsv', infile)
normalized_table.to_csv(outfile, index=False, sep='\t')
