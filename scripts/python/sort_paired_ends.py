#!/usr/bin/env python

import os.path
import sys

def main():
  if len(sys.argv) != 4:
    sys.stderr.write("Usage: ./%s <read_ends_one> <read_ends_two> <output_prefix>\n" % sys.argv[0])
    sys.exit(1)

  try:
    file1 = open(sys.argv[1])
  except:
    sys.stderr.write("Error, couldn't open %s" % sys.argv[1])
    sys.exit(1)

  try:
    file2 = open(sys.argv[2])
  except:
    sys.stderr.write("Error, couldn't open %s" % sys.argv[2])
    sys.exit(1)

  try:
    out1 = open("%s_paired_1.fa" % sys.argv[3], "w")
    out2 = open("%s_paired_2.fa" % sys.argv[3], "w")
    out3 = open("%s_single.fa" % sys.argv[3], "w")
  except:
    sys.stderr.write("Error, couldn't open prefix files at %s" % sys.argv[3])
    sys.exit(1)

  reads = {}
  fasta1 = Fasta(file1)
  fasta2 = Fasta(file2)

  print "Reading first file"
  for entry in fasta1.entries():
    split_name = entry[0].split("/")
    read_name = split_name[0]
    read_end = split_name[1][0]
    if not read_name in reads:
      reads[read_name] = ["",""]
    if read_end == '1' and len(entry[1]) > len(reads[read_name][0]):
      reads[read_name][0] = entry[1]
    elif read_end == '2' and len(entry[1]) > len(reads[read_name][1]):
      reads[read_name][1] = entry[1]

  print "Reading second file"
  for entry in fasta2.entries():
    split_name = entry[0].split("/")
    read_name = split_name[0]
    read_end = split_name[1][0]
    if not read_name in reads:
      reads[read_name] = ["",""]
    if read_end == '1' and len(entry[1]) > len(reads[read_name][0]):
      reads[read_name][0] = entry[1]
    elif read_end == '2' and len(entry[1]) > len(reads[read_name][1]):
      reads[read_name][1] = entry[1]

  print "Writing new files"
  for read_name, read in reads.items():
    if len(read[0]) > 0 and len(read[1]) > 0:
      out1.write(">%s/1\n%s\n" % (read_name, read[0]))
      out2.write(">%s/2\n%s\n" % (read_name, read[1]))
    elif len(read[0]) > 0:
      out3.write(">%s/1\n%s\n" % (read_name, read[0]))
    elif len(read[1]) > 0:
      out3.write(">%s/2\n%s\n" % (read_name, read[1]))

class Fasta:
  def __init__(self, file):
    self.__fasta_f = file
    self.__last_line = self.__fasta_f.readline()
  def next(self):
    result = []
    seq = []
    if len(self.__last_line) == 0:
      return None
    while self.__last_line[0] != '>':
      self.__last_line = self.__fasta_f.readline()
      if len(self.__last_line) == 0:
        return None
    result.append(self.__last_line[1:].strip())

    self.__last_line = self.__fasta_f.readline()
    if len(self.__last_line) == 0:
      return None
    while self.__last_line[0] != '>':
      seq.append(self.__last_line.strip())
      self.__last_line = self.__fasta_f.readline()
      if len(self.__last_line) == 0:
        break
    result.append("".join(seq))
    return result
  def entries(self):
    entry = self.next()
    while entry != None:
      yield(entry)
      entry = self.next()

if __name__ == "__main__":
  main()
