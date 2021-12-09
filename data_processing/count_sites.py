import pandas as pd
import subprocess
import argparse
import re
import gzip

parser = argparse.ArgumentParser(description="run for the lineage",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                )
parser.add_argument('--lin', type=str, default=None, help='lineage for which to run.')
args = parser.parse_args()
lin = args.lin

vcf = '/scratch/drabosky_root/drabosky1/sosi/brazil/variants/%s.qual_filtered20.cov_filtered10.vcf.gz' % lin
f = gzip.open(vcf, 'r')

d = {}

for l in f:
    l = l.decode('utf-8').rstrip()
    if re.search('^#CHROM', l):
        inds = re.split('\t', l)[9:]
        for ind in inds:
            d[ind] = {'loci': {}, 'num_sites': 0, 'tot_cov': 0}
    elif not re.search('^#', l):
        x = re.split('\t', l)
        dpix = re.split(':', x[8]).index('DP')

        for ind, g in zip(inds, x[9:]):
            dp = int(re.split(':', g)[dpix])

            if dp >= 10:
                d[ind]['loci'][x[0]] = 1
                d[ind]['num_sites'] += 1
                d[ind]['tot_cov'] += dp
f.close()

for ind in d:
    print("%s,%s,%s,%s" % (ind, len(d[ind]['loci']), d[ind]['num_sites'], d[ind]['tot_cov']))
