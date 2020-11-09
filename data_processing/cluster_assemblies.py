import re
import pandas as pd
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description="Run simple homology search.")
parser.add_argument('--cl', help="Cluster for which to run homology search.")
args = parser.parse_args()
cl = args.cl

# clustering file
c_file = '/Volumes/heloderma4/sonal/encelia/encelia_samples_v2.csv'
# directory
dir = '/Volumes/heloderma4/sonal/encelia/'
vsearch = 'vsearch'

WCLUST = 0.95

def get_clusters(c_file, cl):
        '''
        get the inds in a given cluster
        '''
        d = pd.read_csv(c_file)
        
        inds = d.ix[d.species == cl, "sample"].tolist()
       
        return inds


def get_homolog(cluster, inds, dir):
	subdir = dir + 'ind_assemblies/'
	outdir = dir + 'cluster_assemblies/'
	tmp1 = '%s%s.tmp1.fa' % (outdir, cluster)
	tmp2 = '%s%s.tmp2.fa' % (outdir, cluster)
	tmp3 = '%s%s.tmp3.fa' % (outdir, cluster)
	tmp4 = '%s%s.tmp4.fa' % (outdir, cluster)
	final = '%s%s.fasta' % (outdir, cluster)

	inds = ['%s%s.fasta' % (subdir, ind) for ind in inds]
	inds2 = []
	for ind in inds:
		if os.path.isfile(ind):
			inds2.append(ind)
	new_inds = []

	for ind in inds2:
		indtmp1 = ind + '_1'
		indtmp2 = ind + '_2'
		indtmp3 = ind + '_3'
		subprocess.call("%s --derep_fulllength %s --output %s --fasta_width 0 --strand both" % (vsearch, ind, indtmp1), shell=True)
		subprocess.call("%s --sortbylength %s --output %s" % (vsearch, indtmp1, indtmp2), shell=True)
        	subprocess.call("%s --cluster_smallmem %s --centroids %s --id %s --usersort --fasta_width 0 --strand both --minsl 0.5 --query_cov 0.7" % (vsearch, indtmp2, indtmp3, WCLUST), shell=True)
		new_inds.append(indtmp3)
		os.remove(indtmp1)
		os.remove(indtmp2)

	subprocess.call("cat %s > %s" % (' '.join(new_inds), tmp1), shell=True)
	call = [os.remove(x) for x in new_inds]
	subprocess.call("%s --derep_fulllength %s --output %s --sizeout --fasta_width 0 --strand both" % (vsearch, tmp1, tmp2), shell=True)
	subprocess.call("%s --sortbylength %s --output %s" % (vsearch, tmp2, tmp3), shell=True)
	subprocess.call("%s --cluster_smallmem %s --centroids %s --sizein --sizeout --id %s --usersort --fasta_width 0 --strand both --minsl 0.5 --query_cov 0.7" % (vsearch, tmp3, tmp4, WCLUST), shell=True)

	# the number of times a contig should be sampled
	ninds = len(new_inds)
	if ninds < 3:
		cl_depth = 1
	elif ninds in [3, 4, 5]:
		cl_depth = 2
	else:
		cl_depth = round(ninds / 3.0)


	f = open(tmp4, 'r')
	o = open(final, 'w')
	ix = 0
	for l in f:
	        if re.search('>', l):
			seq = f.next().rstrip()
                	size = int(re.search('size=(\d+)', l).group(1))
                	if size >= cl_depth:
				o.write('>%s_%s\n%s\n'% (cluster, ix, seq))
				ix += 1
	f.close()
	o.close()

	# os.remove(tmp1)
	# os.remove(tmp2)
	# os.remove(tmp3)
				
inds = get_clusters(c_file, cl)
get_homolog(cl, inds, dir)
