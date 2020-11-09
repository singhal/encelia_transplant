import argparse
import os
import subprocess
import pandas as pd

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
parser.add_argument('--dir', help="directory to use")
args = parser.parse_args()

ind = args.ind
dir = args.dir

samtools = 'samtools'
bwa = 'bwa'
picard = '/Volumes/heloderma4/sonal/bin/picard.jar'
gatk = '/Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar'

c_file = os.path.join(dir, 'encelia_samples_v4.csv')
seq_dir = os.path.join(dir, 'ref_bias')
read_dir = os.path.join(dir, 'trim_reads')
out_dir = os.path.join(dir, 'ref_alignments')

def get_cluster(c_file, ind):
	d = pd.read_csv(c_file)
	cl = d[d['sample'] == ind].lineage.tolist()[0]
	return cl


def prepare_seq(cl, seq_dir):
	seq = os.path.join(seq_dir, '%s.fasta' % cl)
	if not os.path.isfile(seq + '.bwt'):
		subprocess.call("%s index %s" % (bwa, seq), shell=True)
	if not os.path.isfile(seq + '.fai'):
		subprocess.call("%s faidx %s" % (samtools, seq), shell=True)
	if not os.path.isfile(seq.replace('.fasta', '.dict')):
		subprocess.call("java -jar %s CreateSequenceDictionary R=%s O=%s" % (picard, seq, seq.replace('.fasta', '.dict')), shell=True)


def align_seq(ind, cl, out_dir, read_dir, seq_dir):
	r1 = '%s/%s_R1.final.fq.gz' % (read_dir, ind)
	r2 = '%s/%s_R2.final.fq.gz' % (read_dir, ind)
	rU = '%s/%s_unpaired.final.fq.gz' % (read_dir, ind)
	seq = '%s/%s.fasta' % (seq_dir, cl)

	out1 = '%s/%s.sam' % (out_dir, ind)
	out1b = '%s/%s_u.sam' % (out_dir, ind)
	out2 = '%s/%s.mateFixed.bam' % (out_dir, ind)
	out3a = '%s/%s.mateFixed.sorted1.bam' % (out_dir, ind)
	out3b = '%s/%s.mateFixed.sorted2.bam' % (out_dir, ind)
	out3 = '%s/%s.mateFixed.sorted.bam' % (out_dir, ind)
	out4 = '%s/%s.rg.mateFixed.sorted.bam' % (out_dir, ind)

	tmpdir = '%s/%s/' % (out_dir, ind)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call("%s mem -t 4 %s %s %s > %s" % (bwa, seq, r1, r2, out1), shell=True)
	subprocess.call("%s mem -t 4 %s %s > %s" % (bwa, seq, rU, out1b), shell=True)
	# fixmate
	subprocess.call("java -jar %s FixMateInformation I=%s O=%s" % (picard, out1, out2), shell=True)
	# sorted
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (samtools, out3a, tmpdir, out2), shell=True)
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (samtools, out3b, tmpdir, out1b), shell=True)
	# merge
	subprocess.call("%s merge %s %s %s" % (samtools, out3, out3a, out3b), shell=True)
	# readgroup
	subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % (picard, out3, out4, ind, ind, ind), shell=True)
	subprocess.call("%s index %s" % (samtools, out4), shell=True)
	# indeltarget
	# subprocess.call("java -Xmx10g -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt 4" % (gatk, seq, out4, intervals), shell=True)
	# indelrealigner
	# subprocess.call("java -Xmx10g -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (gatk, seq, out4, intervals, out5), shell=True)

	call = [os.remove(x) for x in [out1, out1b, out2, out3a, out3b, out3]]
	os.rmdir(tmpdir)

# get cluster
cl = get_cluster(c_file, ind)
# prepare seqfiles
prepare_seq(cl, seq_dir)
# align it all the way until time to call SNPs
align_seq(ind, cl, out_dir, read_dir, seq_dir)
