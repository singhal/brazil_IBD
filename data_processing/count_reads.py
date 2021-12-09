import pandas as pd
import subprocess
import argparse
import re

parser = argparse.ArgumentParser(description="run for the ind",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                )
parser.add_argument('--ind', type=str, default=None, help='individual for which to run.')
args = parser.parse_args()
ind = args.ind

d = pd.read_csv("/scratch/drabosky_root/drabosky1/sosi/brazil/brazil_samples_v7.csv")

row = d.ix[d['sample'] == ind, ].to_dict('list')
orig_reads = [row['read1'][0], row['read2'][0]]
outdir = '/scratch/drabosky_root/drabosky1/sosi/brazil/trim_reads/'
new_reads = ['%s%s_R1.final.fq.gz' % (outdir, ind),
             '%s%s_R2.final.fq.gz' % (outdir, ind),
             '%s%s_unpaired.final.fq.gz' % (outdir, ind)]

def run_reads(reads):
    tot = 0
    bp = 0
    for r in reads:
        p = subprocess.Popen("~/bin/count_reads/kseq_test %s" % r, stdout = subprocess.PIPE, shell=True)
        x = [l.rstrip() for l in p.stdout]
        tot += int(re.search('(\d+)', x[0]).group(1))
        bp += int(re.search('(\d+)', x[1]).group(1))
    return tot, bp

tot1, bp1 = run_reads(orig_reads)
tot2, bp2 = run_reads(new_reads)

print '%s,%s,%s,%s,%s' % (ind, tot1, bp1, tot2, bp2)
