#!/usr/bin/python

'''
This scripts copies the ID and Sequence
from a fastq file and write them into a new
fasta-like file.

usage fastq_2_fasta.py <infile> <outfile>
'''

import sys

def main():
    out =    open (sys.argv[1], 'w')
    infile = open (sys.argv[2], 'r')

    with infile:
        while True:    
            id1  = infile.readline()
            seq  = infile.readline()
            id2  = infile.readline()
            qual = infile.readline()
                
            out.write (id1)
            out.write (seq)

    out.close()
    infile.close()

if __name__ == "__main__":
    main ()







