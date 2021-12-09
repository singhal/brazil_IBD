import re
import pandas
import os

dir = '/scratch/drabosky_root/drabosky1/sosi/brazil/'
file = '/scratch/drabosky_root/drabosky1/sosi/brazil/brazil_samples_v6.csv'
trimjar = '/home/sosi/bin/Trimmomatic-0.39/trimmomatic-0.39.jar'

name = 'mitoz'
time = '20:00:00'
cpu = 16
d = pandas.read_csv(file)
mem_cpu = 8

for ix, row in d.iterrows():
    ind = row['sample']

    sh = "%s%s.sh" % (name, ix)
    o = open(sh, 'w')
    o.write("#!/bin/bash\n")
    o.write("#SBATCH --job-name %s%s\n" % (name, ix))
    o.write("#SBATCH --nodes=1\n")
    o.write("#SBATCH --cpus-per-task=%s\n" % cpu)
    o.write("#SBATCH --mem-per-cpu=%sg\n" % mem_cpu)
    o.write("#SBATCH --time=%s\n" % time)
    o.write("#SBATCH --account=drabosky1\n")
    o.write("#SBATCH --partition=standard\n")
    o.write("#SBATCH --mail-type=NONE\n")

    o.write("module unload python2.7-anaconda\n")
    o.write("module load Bioinformatics\n")
    o.write("module load samtools/1.9\n")
    # o.write("python ~/brazil/remove_adaptor.py --trimjar %s --dir %s --file %s --CPU %s --sample %s\n" % (trimjar, dir, file, cpu, ind))
    o.write("mkdir %s%s\n" % (dir, ind))
    o.write("module load python3.6-anaconda\n")
    o.write("source activate mitozEnv\n")
    o.write("cd %s%s\n" % (dir, ind))

    o.write("python /home/sosi/bin/release_MitoZ_v2.3/MitoZ.py assemble --clade Chordata \
    --outprefix %s \
    --thread_number 16 \
    --fastq1 /scratch/drabosky_root/drabosky1/sosi/brazil/trim_reads/%s_R1.final.fq.gz \
    --fastq2 /scratch/drabosky_root/drabosky1/sosi/brazil/trim_reads/%s_R2.final.fq.gz \
    --fastq_read_length 101 --insert_size 150 --min_abundance 0.9 --filter_taxa_method 3\n" % (ind, ind, ind))
    o.close()

    # outfile = '%s/trim_reads/%s_R1.final.fq.gz' % (dir, ind)
    outfile = '%s/mt_genomes/%s_mitogenome.fasta' % (dir, ind)
    if os.path.isfile(outfile):
        os.remove(sh)
    # else:
    #    print(ind)
