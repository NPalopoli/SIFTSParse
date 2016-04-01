#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Parse pdb_chain_uniprot.tsv from SIFTS to map ELM entries to PDB ss and disorder information
  File name: parseSIFTS_PDBssSEQATOM.py
  Author: Nicolas Palopoli
  Date created: 2015/10/30
  Date last modified: 2015/11/02
  Python Version: 2.7
'''

import sys
import csv
from collections import OrderedDict
from Bio import SeqIO

# Read input files

try:
  inELMinstances = open(sys.argv[1])
  inSIFTSpdbuniprot = open(sys.argv[2])
except IndexError:
  print("Input file not specified. Format: ./parseSIFTS_PDBssSEQATOM.py <elm_instances[.date].tsv> <pdb_chain_uniprot.tsv>")
  exit()
except IOError:
  print("Input file not found. Format: ./parseSIFTS_PDBssSEQATOM.py <elm_instances[.date].tsv> <pdb_chain_uniprot.tsv>")
  exit()

# Functions

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def mapELMpositions(parsedELM):
  '''Make dict with {ELMAccession:[Primary_acc,Start,End,[PDB1,...,PDBN]]}'''
  ELMpos = {}
  for row in parsedELM:
    ELMpos[row['Accession']] = [row['Primary_Acc'],row['Start'],row['End'],row['PDB'].split()]
  return ELMpos

def readSIFTSpdbuniprot(infile):
  '''Store SIFTS information as list of dicts'''
  sifts = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return sifts

def mapSIFTSpdbpositions(parsedSIFTS):
  '''Make dict with {'SP_PRIMARY:PDB':[PDB,CHAIN,RES_BEG,RES_END,PDB_BEG,PDB_END,SP_PRIMARY,SP_BEG,SP_END]}'''
  SIFTSpos = {}
  for row in parsedSIFTS:
    key = '%s:%s' % (row['SP_PRIMARY'],row['PDB'].upper())  # index by SP_PRIMARY:PDB
    SIFTSpos[key] = [row['PDB'].upper(),row['CHAIN'].upper(),row['RES_BEG'],row['RES_END'],row['PDB_BEG'],row['PDB_END'],row['SP_PRIMARY'],row['SP_BEG'],row['SP_END']]
  return SIFTSpos

def makeUniprotGaps(infile,SIFTSentry):
  fastaseqs = SeqIO.parse(open(infile),'fasta')
# with open(output_file) as out_file:
  for fasta in fastaseqs:
#    print "%s\n%s" % (SIFTSentry[6],fasta.id)
    if SIFTSentry[6] in fasta.id:
#      SIFTSentry = ">%s:%s" % (line.id[0:6],'disorder')
      seqgaps = '-' * len(fasta.seq)
  return seqgaps
    
def getPDBpos(ELMpos,SIFTSpdbpos):
  ''' '''
  for keyELM,valELM in ELMpos.items():
    for PDB in valELM[3]:  # list of PDBs with the ELM
      target = '%s:%s' % (valELM[0],PDB)  # PrimaryAcc:PDB 
      if target in SIFTSpdbpos:  # subset of ELMs where PDB sequence is included in UniprotID -- need to extend this for cases where PDB is not in UniprotID
        if SIFTSpdbpos[target][0] == PDB and SIFTSpdbpos[target][6] == valELM[0] and int(SIFTSpdbpos[target][7]) <= int(valELM[1]) and int(SIFTSpdbpos[target][8]) >= int(valELM[2]):  # PDB and Chain match, within limits
#        if SIFTSpdbpos[target][0] == PDB and SIFTSpdbpos[target][6] == valELM[0] and SIFTSpdbpos[target][7] <= valELM[1] and SIFTSpdbpos[target][8] >= valELM[2]:  # PDB and Chain match, within limits
          print parsePDBss('../PDBss/PDB_ss_dis_SEQATOM_all.txt',SIFTSpdbpos[target][0],SIFTSpdbpos[target][1],SIFTSpdbpos[target][2],SIFTSpdbpos[target][6],SIFTSpdbpos[target][7],keyELM)
#  return PDBss

def parsePDBss(infile,pdb,chain,resbeg,uniprot,spbeg,keyELM):
  '''Read 2nd struct from PDB_ss_dis_SEQATOM_all.txt (unified from PDB files ss.txt and ss_dis.txt and with SEQRES and SEQATOM from SEQATOMs db) in FASTA format'''
  resbeg = int(resbeg)
  seqgaps = '-' * (int(spbeg)-1)
#  seqgaps = '-' * (int(spbeg)-1+resbeg-1)
#  seqgaps = '-' * (int(spbeg)+resbeg-1)
  seq = '-'
  secstr = '-'
  disorder = '-'
  seqres = '-'
  seqatom = '-'
  fastaseqs = SeqIO.parse(open(infile),'fasta')
  for fasta in fastaseqs:
    if fasta.id[0:4] == pdb and fasta.id[5] == chain:
      if 'sequence' in fasta.id:
        seq = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'secstr' in fasta.id:
        secstr = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'disorder' in fasta.id:
        disorder = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'SEQRES' in fasta.id:
        seqres = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'SEQATOM' in fasta.id:
        seqatom = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
  name = "%s:%s:%s:%s" % (keyELM,uniprot,pdb,chain)
  res = ">%s:sequence\n%s\n>%s:secstr\n%s\n>%s:disorder\n%s\n>%s:SEQRES\n%s\n>%s:SEQATOM\n%s\n" % (name,seq,name,secstr,name,disorder,name,seqres,name,seqatom)
#  res = ">%s:%s\n%s\n>%s:%s\n%s\n>%s:%s\n%s\n" % (uniprot,pdb,seq,uniprot,pdb,secstr,uniprot,pdb,disorder)
  return res

# START

# Make list of one dict per ELM instance in input file
parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()

# Get ELM start and end positions on UniProt Primary Accession, and get all PDBs annotated for the ELM
ELMpos = mapELMpositions(parsedELM)

# Make list of one dict per SIFTS entry in input file
parsedSIFTS = readSIFTSpdbuniprot(inSIFTSpdbuniprot)
inSIFTSpdbuniprot.close()

# Make list of dicts with positions from UniProt Primary Accession to PDB Chain SEQRES
# e.g.: { (...), '3DNL:P04578': ['3DNL', 'E', '1', '67', '330', '396', 'P04578', '330', '396'], (...) }
SIFTSpdbpos = mapSIFTSpdbpositions(parsedSIFTS)

# Get PDB, chain, start and end positions for entries in ELM
PDBss = getPDBpos(ELMpos,SIFTSpdbpos)
#print PDBss

exit()

#PDBss = parsePDBss('../PDBss/PDB_ss_dis.txt',SIFTSpdbpos['Q60795'][0].upper(),SIFTSpdbpos['Q60795'][1].upper())
#PDBss = parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[primaryAcc][0].upper(),SIFTSpdbpos[primaryAcc][1].upper())
PDBss = {}
for key in SIFTSpdbpos:
  print parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[key][0],SIFTSpdbpos[key][1])
  PDBss[key] = parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[key][0],SIFTSpdbpos[key][1])
print PDBss
# {'disorder': '-----------------------------------', 'seq': 'MDLIDILWRQDIDLGVSREVFDFSQRQKDYELEKQ', 'pdb': '3WN7', 'chain': 'M', 'secstr': '--HHHHHHTTTGGGT--GGGG--------------'}


