#!/usr/bin/python
# This is the master python script for the bootstrap tool

import numpy as np
# import scipy as sp
import sys
from StringIO import StringIO
from datetime import datetime

tstart = datetime.now();

def gen_samples(c1a, ac2, c1c2, K):
    "Generate bootstrap samples of PSI"
    n = c1a.size; #print n
    idx1 = np.random.randint(0,n,[n,K]);
    idx2 = np.random.randint(0,n,[n,K]);
    idx3 = np.random.randint(0,n,[n,K]);
    ninc = np.sum(c1a[idx1],axis=0) + np.sum(ac2[idx2],axis=0);
    nexc = np.sum(c1c2[idx3],axis=0);
    incl = np.zeros(K); excl = np.zeros(K);
    for i in np.arange(K):
        incl[i] = np.random.gamma(ninc[i]/2+1);
        excl[i] = np.random.gamma(nexc[i]+1);
    
    all_l = incl + excl;    
    psit = incl / all_l;
    return psit;

print 'Argument List:', str(sys.argv)

K = 10000; # number of bootstrap samples
suffix = "psi_bootstrap";
file_result = suffix + ".result";
file_sample = suffix + ".sample";
header = "#ID\tPSI_standard\tPSI_bootstrap\tPSI_bootstrap_std\n";
#print file_result, header
f_r = open(file_result, 'w');
f_r.write(header); #f_r.close();

file_junc_rd = 'triplet_junc_rd.tsv';
data_rd = np.genfromtxt(file_junc_rd, delimiter="\t", dtype="O,f,f,f,O,O,O,O", 
                        names=True);
ID = data_rd['ID']; N = ID.size;
n_C1A = data_rd['n_C1A']; n_AC2 = data_rd['n_AC2']; n_C1C2 = data_rd['n_C1C2'];
n_inc = (n_C1A + n_AC2)/2;
psi_standard = (n_inc + 1)/(n_inc + n_C1C2 + 2);
psi_bootstrap = np.zeros(N); psi_b_std = np.zeros(N);
#psi_b_std = psi_bootstrap; #this will make both the same all the time!
rd_C1A = data_rd['RD_C1A']; rd_AC2 = data_rd['RD_AC2']; 
rd_C1C2 = data_rd['RD_C1C2'];

for i in np.arange(N):
    c1a = np.genfromtxt(StringIO(rd_C1A[i].strip(",")), delimiter=",");
    # print i, c1a;
    ac2 = np.genfromtxt(StringIO(rd_AC2[i].strip(",")), delimiter=",");
    c1c2 = np.genfromtxt(StringIO(rd_C1C2[i].strip(",")), delimiter=",");
    # generate bootstrap samples here!
    psit = gen_samples(c1a, ac2, c1c2, K);
    # process the samples
    psi_bootstrap[i] = np.mean(psit);
    psi_b_std[i] = np.std(psit);
    str_out = ("%s\t%f\t%f\t%f\n" % (ID[i], psi_standard[i], 
                                     psi_bootstrap[i], psi_b_std[i]));
    f_r.write(str_out);                            
    
    # print progress
    if ((i+1) % 100 == 0): 
        print "%d out of %d triplets have been processed..." % (i+1, N+1)

f_r.close()

tend = datetime.now();
tdiff = tend - tstart;
runtime = tdiff.seconds/60.0;
print "Finished processing all %d triplets in %g minutes!" % (N+1, runtime)
