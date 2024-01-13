#!/usr/bin/env python3

from Bio import SeqIO
import sys
import argparse

#help file that deliniates required arguments
parser = argparse.ArgumentParser(description='Enter a list of headers to remove and a fasta file to remove them from.')
parser.add_argument('--headers', required=True, help='file of headers to remove')
parser.add_argument('--input', required=True, help='file of fasta sequences')
parser.add_argument('--out', required=True, help='file of results')
args=parser.parse_args()

#print warning if not enough arguments
if len(sys.argv) < 3: 
  print("\nYou need to provide: \n(1) a list of headers to remove\n(2) a fasta file to remove them from\n(3) output file.")
  sys.exit(0)

# print warning if too many arguments.
if len(sys.argv) > 4: 
  print("\nExtra arguments input!")
  sys.exit(0)

#identify each of the headers

headers = set(line.strip() for line in open(sys.argv[2]))

#parse through all the sequences and remove any with a header matching the inputted headers to be removed
for seqs in SeqIO.parse(sys.argv[1], "fasta"):
  if seqs.id not in headers:
    SeqIO.write(seqs, args.out, "fasta")
  try:
      headers.remove(seqs.name)
  except KeyError:
      print(seqs.format("fasta"))
      continue

#report any headers that were not found in the fasta file and therefore not removed
if len(headers) != 0:
  print(len(headers),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)
