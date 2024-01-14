#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse
import pandas as pd
import sys

#help file of arguments
parser = argparse.ArgumentParser(description='Calculates GC content of every sequence in fasta.')
parser.add_argument('--input', required=True, help='file of fasta sequences')
parser.add_argument('--out', required=True, help='output file (csv)')
args=parser.parse_args()

#print warning if no arguments given
if not len(sys.argv) > 1: 
  parser.print_help()
  sys.exit(0)
return args

#pull out sequences from fasta
sequences = SeqIO.parse(sys.argv[1], "fasta")

#print warning if there are no sequences in the fasta file
if len(sequences) < 1: 
  print("\nError: no sequences were detected in your input fasta file.")
  sys.exit(0)

#parse through all the sequences and calculate GC content
ids = []
contents = []
for seqs in sequences:
  ids.append(seq.id)
  contents.append(GC(seq.seq))

#add to dataframe
d = {'id':ids,'GC_%':contents}
df = pd.DataFrame(data=d)

#extract as csv
df.to_csv(args.out,index=False)
