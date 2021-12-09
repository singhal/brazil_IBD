
sigmas = [2.0]
nrep = 5
pops = [20]
basedir = '/Users/singhal/Desktop/sims/'
sfile = '/Users/singhal/Desktop/sims/spaceness.slim'

for i in range(0, nrep):
    for sigma in sigmas:
        for pop in pops:
            outdir = basedir + 'out%s_n%s_%s/' % (sigma, pop, i)
            print("mkdir %s" % outdir)
            print("slim -d sigma=%s -d pop=%s -d fpath=\'\"%s\"\' %s" % (sigma, pop, outdir, sfile))
