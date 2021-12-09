import argparse
import gzip
import os
import numpy as np
import pandas as pd
import re
import subprocess

'''
Sonal Singhal
created on 29 June 2016
'''


def get_args():
	parser = argparse.ArgumentParser(
		description="Calculate pi for a sample.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# ample
	parser.add_argument(
		'--sample',
		type=str,
		default=None,
		help='samplee for which to make calculations.'
		)

	# sample file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with sample info.'
		)
		
	# base dir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help="Base directory as necessary"
			 " when used with pipeline"
		)

	# output dir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Directory to output pop gen stats, '
			 'only necessary if not running '
			 'in context of pipeline'
		)

	# vcfdir
	parser.add_argument(
		'--vcfdir',
		type=str,
		default=None,
		help='Directory with VCFs, '
			 'only necessary if not running '
					 'in context of pipeline'
		)

	return parser.parse_args()


def get_diversity(ind, vcf, outdir):
	# keep track of it all
	all  = {'sum': 0, 'sites': 0}

	allowed = ['0/0', '0/1', '1/1']

	f = gzip.open(vcf, 'r')
	for l in f:
		l = l.decode('utf-8')
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				geno = re.search('^(\S\/\S)', d[9]).group(1)
				if geno in allowed:
					all['sites'] += 1
					if geno == '0/1':
						all['sum'] += 1
	f.close()
		
	out = os.path.join(outdir, '%s_diversity.csv' % ind)
	o = open(out, 'w')

	o.write('ind,het,het_denom\n')
	if all['sites'] > 0:
		het = all['sum'] / float( all['sites'] )
	o.write('%s,%.5f,%s\n' % (ind, het, all['sites']))
	o.close()
			

def get_data(args):
	ind = args.sample

	if args.dir:
		outdir = os.path.join(args.dir, 'pop_gen')
		vcf = os.path.join(args.dir, 'ind_variants', 
								   '%s.qual_filtered20.cov_filtered10.vcf.gz' % ind)
	else:
		outdir = args.outdir
		vcf = os.path.join(args.vcfdir, 
								   '%s.qual_filtered20.cov_filtered10.vcf.gz' % ind)

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return ind, vcf, outdir


def main():
	args = get_args()
	ind, vcf, outdir = get_data(args)
	get_diversity(ind, vcf, outdir)


if __name__ == "__main__":
	main()
