#!/usr/bin/python
# This is the python script to compute bootstrap PSI

import numpy as np
# import scipy as sp
import sys
from StringIO import StringIO
from datetime import datetime

tstart = datetime.now();
    
def gen_pdf(inc, exc, K,G):
    "Generate bootstrap PDF of PSI"
    grid=np.arange(1./(2*G),1.,1./G);
    a=1.; b=1.;r=0.; 
    pinc = inc.size; #print n
    pexc = exc.size; #print n
    i_inc = np.random.randint(0,pinc,[pinc,K]);
    i_exc = np.random.randint(0,pexc,[pexc,K]);
    
    ninc = np.sum(inc[i_inc],axis=0);
    nexc = np.sum(exc[i_exc],axis=0);
    pdf = np.zeros(G);
    #pdf.dtype
    for i in np.arange(K):
        logpdf=(ninc[i]+a-1)*np.log(grid)+(nexc[i]+b-1)*np.log(1-grid)-(ninc[i]+nexc[i])*np.log(grid*pinc+(1-grid)*pexc+r);
        maxpdf=np.max(logpdf);
        logpdf=logpdf-maxpdf;
        pdf_current=np.exp(logpdf);
        pdf_current=np.divide(pdf_current,np.sum(pdf_current));
        pdf=np.add(pdf,pdf_current);
    pdf=np.divide(pdf,K);    
    return (pdf,grid);
#
#print 'Argument List:', str(sys.argv)

K = 1000; # number of bootstrap samples
G = 100; #number of points in PSI pdf

# path =r"";
# file_result = path + "psi_bootstrap.result";
# file_sample = path + "psi_bootstrap.sample";
suffix = "psi_bootstrap";
file_result = suffix + ".result";
file_sample = suffix + ".sample";
header = "#ID\tn_inc\tn_exc\tp_inc\tp_exc\tPSI_standard\tPSI_bootstrap\tPSI_bootstrap_std\n";
#print file_result, header
f_r = open(file_result, 'w');
f_r.write(header); #f_r.close();

# file_junc_rd = path+'event_inc_exc.rd';
file_junc_rd = 'event_inc_exc.rd';
data_rd = np.genfromtxt(file_junc_rd, delimiter="\t", dtype="O,O,O", 
                        names=True);
ID = data_rd['ID'];
n_inc= data_rd['RD_inc']; n_exc = data_rd['RD_exc'];
N=np.size(n_inc);

psi_bootstrap = np.zeros(N); psi_b_std = np.zeros(N);psi_standard = np.zeros(N);
#psi_b_std = psi_bootstrap; #this will make both the same all the time!

for i in np.arange(N):
    inc = np.genfromtxt(StringIO(n_inc[i].strip(",")), delimiter=",");
    exc = np.genfromtxt(StringIO(n_exc[i].strip(",")), delimiter=",");
    p_inc=np.size(inc);
    p_exc=np.size(exc);
    min_p=min(p_inc,p_exc);
    scaled_inc=np.sum(inc)/p_inc*min_p;
    scaled_exc=np.sum(exc)/p_exc*min_p;
    psi_standard[i]=(scaled_inc+1)/(scaled_inc+scaled_exc+2);
    
    # generate bootstrap PDF here!
    pdf,grid = gen_pdf(inc,exc,K,G);
    # process the PDF to get PSI
    psi_bootstrap[i] = sum(pdf*grid);
    psi_b_std[i] = np.sqrt(sum(pdf*np.square(grid-psi_bootstrap[i])));
    str_out = ("%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n" % (ID[i], np.sum(inc), np.sum(exc), p_inc,p_exc, psi_standard[i], psi_bootstrap[i], psi_b_std[i]));
    f_r.write(str_out);                            
    
    # print progress
    if ((i+1) % 1000 == 0): 
        print "%d out of %d AS events have been processed..." % (i+1, N)

f_r.close()

tend = datetime.now();
tdiff = tend - tstart;
runtime = tdiff.seconds/60.0;
print "Finished processing all %d AS events in %g minutes!" % (N, runtime)
