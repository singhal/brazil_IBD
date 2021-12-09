import re
import glob
import os
import pandas as pd
import subprocess

genes = {'ND4': 'c', 'COI': 'c', 'CYTB': 
		'c', 'ND1': 'c', 'ND2': 'c', 'rRNA_12S': 
		'nc', 'rRNA_16S': 'nc', 'ATP6': 'c', 
		'ATP8': 'c', 'COII': 'c', 'COIII': 'c',
		'ND3': 'c', 'ND4': 'c', 'ND4L': 'c',
		'ND5': 'c', 'ND6': 'c'}

indir = '/scratch/drabosky_root/drabosky1/sosi/brazil/mt_genomes_anno/'
outdir = '/scratch/drabosky_root/drabosky1/sosi/brazil/mtDNA_loci/'

def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)

	return bases

for gene in genes:
	outfiles = glob.glob('%s*%s' % (indir, gene))
	matches = {}

	print(gene)
	for ix, outfile in enumerate(outfiles):
		sample = re.sub(outdir, '', outfile)
		sample = re.sub('_mitogenome..*', '', sample)

		seqfile = re.sub('.fasta.*$', '.fasta', outfile)
		seqfile = re.sub('mt_genomes_anno', 'mt_genomes', seqfile)


		o = open(outfile, 'r')
		for l in o:
			if re.search('similarity', l):
				d = re.split('\t', l.rstrip())
				if sample in matches:
					length = matches[sample]['end'] - matches[sample]['start']
					new_len_diff = (int(d[4]) - int(d[3])) / float(length)
					if new_len_diff >= 0.9 and int(d[5]) > matches[sample]['score']:
						matches[sample] = {'start': int(d[3]), 'end': int(d[4]),
										  'seqid': d[0], 'score': int(d[5]),
										  'seq': seqfile, 'ind': sample, 'orr': d[6]}
				else:
					matches[sample] = {'start': int(d[3]), 'end': int(d[4]), 
										'seqid': d[0], 'score': int(d[5]), 
										'seq': seqfile, 'ind': sample,
										'orr': d[6]}
		o.close()

	out = os.path.join(outdir, '%s.fasta' % gene)
	o = open(out, 'w')
	for sample in matches:
		call = subprocess.Popen("samtools faidx %s %s:%s-%s" % (matches[sample]['seq'], 
			matches[sample]['seqid'], matches[sample]['start'], matches[sample]['end']), 
			shell=True, stdout=subprocess.PIPE)
		lines = [line.rstrip() for line in call.stdout]


		o.write('>%s\n' % (sample))
		seq = ''
		for l in lines[1:]:
			seq += l
		if matches[sample]['orr'] == '-':
			seq = reverse_complement(seq)
		o.write('%s\n' % seq)
	o.close()