#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import sys

#help file that deliniates required arguments
parser = argparse.ArgumentParser(description='Enter a list of headers to remove and a fasta file to remove them from.')
parser.add_argument('--headers', required=True, help='file of headers to remove')
parser.add_argument('--input', required=True, help='file of fasta sequences')
parser.add_argument('--out', required=True, help='file of results')
args=parser.parse_args()

#print warning if no arguments given
if not len(sys.argv) > 1: 
  parser.print_help()
  sys.exit(0)
return args

#identify each of the headers
headers = set(line.strip() for line in open(args.headers))

#parse through all the sequences and remove any with a header matching the inputted headers to be removed
for seqs in SeqIO.parse(args.input, "fasta"):
  headers.remove(seqs.name)
  if seqs.id not in headers:
    SeqIO.write(seqs, args.out, "fasta")
  
#report any headers that were not found in the fasta file and therefore not removed
if len(headers) != 0:
  print(len(headers),'of the headers were not in the fasta file.', file=sys.stderr)
