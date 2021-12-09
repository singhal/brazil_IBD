import argparse
import glob
import os
import pandas as pd
import re
import subprocess

# where the alignments are
locfiles = glob.glob('/scratch/drabosky_root/drabosky1/sosi/brazil/mtDNA_loci/*aln')

file1 = '/scratch/drabosky_root/drabosky1/sosi/brazil/mtDNA_loci/concatenated.fasta'

loci = {}
seq = {}

for f in locfiles:
        locus = re.sub('^.*/', '', f)
        locus = re.sub('.fasta', '', locus)

        f = open(f, 'r')
        id = ''
        s = {}
        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l).group(1)
                        # get rid of reverse
                        id = re.sub('^_R_', '', id)
                        id = re.sub('/scratch/drabosky_root/drabosky1/sosi/brazil/mt_genomes_anno/', '', id)
                        s[id] = ''
                else:
                        s[id] += l.rstrip()
        f.close()

        for sp, tmpseq in s.items():
                tmpseq = re.sub('\s+', '', tmpseq)
                if sp not in seq:
                        seq[sp] = {}
                seq[sp][locus] = tmpseq


        loci[locus] = len(s[s.keys()[0]])


ordloc = sorted(list(loci.keys()))
f = open(file1, 'w')
for sp in seq:
        totseq = ''
        for loc in ordloc:
                if loc in seq[sp]:
                        totseq += seq[sp][loc]
                else:
                        totseq += '-' * loci[loc]
        print(len(totseq))
        f.write('>%s\n%s\n' % (sp, totseq))
f.close()
