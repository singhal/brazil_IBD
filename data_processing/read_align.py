import pandas as pd
import subprocess
import argparse
import re

d = pd.read_csv("/scratch/drabosky_root/drabosky1/sosi/brazil/brazil_samples_v7.csv")
lins = d['lineage'].value_counts()
tlins = lins[lins > 3].index.tolist()
d1 = d.ix[d.lineage.isin(tlins), ]

def calc_cover(bam):
    p = subprocess.Popen("samtools flagstat %s" % bam, stdout = subprocess.PIPE, shell=True)
    x = [l.rstrip() for l in p.stdout]
    tot = re.search('([\d+|\.]+)\%', x[4]).group(1)
    return tot

out = '/home/sosi/brazil/read_align.csv'
o = open(out, 'w')

for ix, row in d1.iterrows():
    ind = row['sample']

    bam1 = "/scratch/drabosky_root/drabosky1/sosi/brazil/alignments/%s.dup.rg.mateFixed.sorted.bam" % ind
    bam2 = "/scratch/drabosky_root/drabosky1/sosi/brazil/alignments/%s.dup.rg.mateFixed.sorted.recal.bam" % ind
    lineage = row['lineage']

    tot1 = calc_cover(bam1)
    tot2 = calc_cover(bam2)

    o.write('%s,%s,%s,%s\n' % (ind, lineage, tot1, tot2))

o.close()
