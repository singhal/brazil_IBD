import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 22 June 2016
Written assuming:
	* bcftools 1.3.1
	* samtools 1.3.1
	* GATK 3.6
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Align reads to lineage, step 2. "
					" Assumes bcftools 1.3.1 "
					" and GATK 4",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
	parser.add_argument(
		'--lineage',
		type=str,
		default=None,
		help='Lineage for which to run script.'
		)

	# file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with sample info.'
		)
				
	# basedir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help="Full path to base dir with reads & assemblies "
			 "everything else."
		)

	# bcftools
	parser.add_argument(
		'--bcftools',
		type=str,
		default=None,
		help='bcftools executable, full path.'
		)

	# GATK
	parser.add_argument(
		'--gatk',
		type=str,
		default=None,
		help='GATK executable, full path.'
		)
	
	# memory
	parser.add_argument(
		'--mem',
		type=int,
		default=1,
		help='Memory available, as an int, in terms of Gb.'
	   )

	# qual
	parser.add_argument(
		'--qual',
		type=int,
		default=20,
		help='Minimum quality to retain variant for '
					 'creating validated call set.'
		)

	# depth
	parser.add_argument(
		'--dp',
		type=int,
		default=10,
		help='Minimum depth required per individual to retain '
			 'variant for creating validated call set.'
		)

	# depth
	parser.add_argument(
		'--CPU',
		type=int,
		default=1,
		help='# of CPUs that can be used in script'
		)
		
	
	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Output directory for alignments, only needed '
			 'if not running in context of pipeline.'
		)
				
	# bamfiles
	parser.add_argument(
		'--bamfile',
		type=str,
		default=None,
		help="Full path to file with BAM files, listed one "
			 "per line if running not in context of pipeline. "
		)

	# PRG
	parser.add_argument(
		'--prg',
		type=str,
		default=None,
		help="Full path to pseudoref genome if "
			 "you aren't running in context of pipeline."
		)

	return parser.parse_args()


def get_files(args):
	# gets the bam files
	if args.bamfile:
		f = open(args.bamfile, 'r')
		files = []
		for l in f:
			files.append(l.rstrip())
		f.close()
	else:
		d = pd.read_csv(args.file)
		samps = d.ix[d['lineage'] == args.lineage, 'sample'].tolist()
		files = []
		for samp in samps:
			file = os.path.join(args.dir, 'alignments', 
											'%s.dup.rg.mateFixed.sorted.bam' % samp)
			files.append(file)
	# makes sure the order stays consistent
	files = sorted(files)

	# gets the genome
	if args.prg:
		genome = args.prg
	else:
		genome = os.path.join(args.dir, 'PRG', '%s.fasta' % args.lineage)

	# gets the outdir
	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'alignments')
	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return files, genome, outdir


def get_qual(args, files, genome, dir):
	raw_vcf = os.path.join(dir, '%s.raw.vcf' % args.lineage)
	filt_vcf = os.path.join(dir, '%s.filt.vcf' % args.lineage)

	bam = ' '.join(files)

	# makes the raw VCFs, only outputting variant SNPs
	subprocess.call("%s mpileup -A -f %s -a DP -Ou %s | %s call -vmO v -o %s" % (args.bcftools, genome, bam, args.bcftools, raw_vcf), shell = True)

	f = open(raw_vcf, 'r')
	o = open(filt_vcf, 'w')

	for l in f:
		if re.match('^#', l):
			o.write(l)
		else:
			d = re.split('\t', l.rstrip())
			
			# check if indel
			alleles = [d[3]] + re.split(',', d[4])
			indel = False
			for a in alleles:
				if len(a) > 1:
					indel = True
			
			# depth cover
			n_inds = len(d[9:])
			dp = args.dp * n_inds
			snp_depth = int(re.search('DP=(\d+)', d[7]).group(1))

			if not indel and float(d[5]) >= args.qual and snp_depth >= dp:
				o.write(l)

	f.close()
	o.close()

	os.remove(raw_vcf)
	# os.remove(raw_vcf + '.idx')
	return filt_vcf


def recalibrate(args, files, genome, vcf, dir):
	subprocess.call("%s CreateSequenceDictionary --REFERENCE %s" % (args.gatk, genome), shell = True)
	subprocess.call("%s IndexFeatureFile --feature-file %s" % (args.gatk, vcf), shell = True)

	for file in files:
		stem = file.replace('.bam', '')
		out = stem + '.recal.bam'
		recal = '%s.recal.table' % stem
		# generate recal table
		subprocess.call("%s BaseRecalibrator --reference %s --known-sites %s --input %s --output %s" %
					(args.gatk, genome, vcf, file, recal), shell=True)
		# print the recal reads
		subprocess.call("%s ApplyBQSR --input %s --bqsr-recal-file %s --output %s" %
				(args.gatk, file, recal, out), shell=True)

		os.remove(recal)
		# os.remove(file)
		# os.remove(file.replace('bam', 'bai'))
	os.remove(vcf)
	os.remove(vcf + '.idx')


def main():
	args = get_args()
	files, genome, outdir = get_files(args)	
	vcf = get_qual(args, files, genome, outdir)
	recalibrate(args, files, genome, vcf, outdir)


if __name__ == "__main__":
	main()
