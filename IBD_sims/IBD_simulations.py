import re
import glob
import math
import numpy as np
import os

# bdir = '/Users/singhal/Dropbox (Personal)/brazil_gene_flow/data/IBD_simulations/'
bdir = '/Users/singhal/Desktop/sims/'

def fst_reich(counts, sample_sizes):
	counts1 = np.array(counts)[0]
	counts2 = np.array(counts)[1]

	sample_sizes1 = np.array(sample_sizes).astype('float')[0]
	sample_sizes2 = np.array(sample_sizes).astype('float')[1]

	h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
	h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))
	
	N = []
	D = []

	for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
		n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
		N.append(n)
		d = n + _h1 + _h2
		D.append(d)

	F = np.sum(N) / np.sum(D)

	return F

def get_divergence(sigma, rep, outdir, inds, inds2, vcf):
	diff = { '0|0': {'0|1': 0.5, '1|1': 1, '0|0': 0, '1|0': 0.5},
	        '0|1': {'0|1': 0, '1|1': 0.5, '0|0': 0.5, '1|0': 0},
	        '1|0': {'0|1': 0, '1|1': 0.5, '0|0': 0.5, '1|0': 0},
	        '1|1': {'0|1': 0.5, '1|1': 0, '0|0': 1, '1|0': 0.5} }
	count = {'0|0': 0, '1|1': 2, '0|1': 1, '.|.': np.nan, '1|0': 1}

	div = {}
	for ix, ind1 in enumerate(inds2):
		div[ind1] = {}
		for ind2 in inds2[(ix + 1):]:
			div[ind1][ind2] = {'diff': 0, 'denom': 0}
		
	# for calculating fst
	counts = dict([(ind, []) for ind in inds2])

	f = open(vcf, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = d[9:]
				genos = [re.search('^(\S\|\S)', x).group(1) for x in genos]

				# variable site to be used in fst
				if d[4] in ['A', 'T', 'C', 'G']:
					for ind, geno in zip(inds2, genos):
						counts[ind].append(count[geno])

				# get divergence data
				genos = dict(zip(inds2, genos))
				for ind1 in div:
					for ind2 in div[ind1]:
						if genos[ind1] != './.' and genos[ind2] != './.':
							div[ind1][ind2]['denom'] += 1
							div[ind1][ind2]['diff'] += diff[genos[ind1]][genos[ind2]]
	f.close()


	out = os.path.join(outdir, '%s_%s_divergence.csv' % (sigma, rep))
	o = open(out, 'w')
	o.write('sigma,rep,ind1,ind2,distance,nuc_dxy,nuc_denom,fst,fst_denom\n')
	for ind1 in div:
		for ind2 in div[ind1]:
			dxy_denom = 5000 * 3000
			if dxy_denom > 0:
				dxy = div[ind1][ind2]['diff'] / float(dxy_denom)
			else:
				dxy = np.nan

			alleles = np.array([counts[ind1], counts[ind2]])
			to_mask = np.any(np.isnan(alleles), axis=0)
			alleles = alleles[:, ~to_mask]
			if len(alleles[0]) > 0:
				sizes = [[2] * len(alleles[0]), [2] * len(alleles[0])]
				fst = fst_reich(alleles, sizes)
			else:
				fst = np.nan

			pt1 = inds[ind1]
			pt2 = inds[ind2]

			dist = math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)

			o.write('%s,%s,%s,%s,%.6f,%.6f,%s,%.6f,%s\n' % 
				(sigma, rep, ind1, ind2, dist, dxy, dxy_denom, fst, len(alleles[0])))

	o.close()

dirs = glob.glob(bdir + "out*")
for dir in dirs:
	info = re.sub(".*\/", "", dir)
	(sigma, pop, rep) = re.findall('([\d+|\.]+)', info)

	pfile = dir + '/positions1.txt'
	inds = []
	p = open(pfile, 'r')
	for l in p:
		if not re.search("index", l):
			d = re.split('\s+', l.rstrip())
			inds.append([float(d[1]), float(d[2])])


	inds2 = list(range(0, 100))
	vcf = dir + '/sample1.vcf'
	outdir = bdir + 'summary/'
	get_divergence(sigma, rep, outdir, inds, inds2, vcf)
