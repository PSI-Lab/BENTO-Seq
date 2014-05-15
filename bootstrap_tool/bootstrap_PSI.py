#!/usr/bin/env python
# This is the python script to compute bootstrap PSI

import numpy as np
from  numpy import newaxis as na
# import scipy as sp
import sys, os
from datetime import datetime
from itertools import izip

tstart = datetime.now()
    
def gen_pdf(inc, exc, n_bootstrap_samples=1000, n_grid_points=100, a=1., b=1., r=0.):
    """Generate bootstrap PDF of PSI"""

    grid = np.arange(1. / (2 * n_grid_points), 1., 1. / n_grid_points)
    pinc = inc.size
    pexc = exc.size
    i_inc = np.random.choice(pinc, (pinc, n_bootstrap_samples))
    i_exc = np.random.choice(pexc, (pexc, n_bootstrap_samples))
    
    ninc = np.sum(inc[i_inc], axis=0)
    nexc = np.sum(exc[i_exc], axis=0)

    logpdf = (ninc + a - 1)[:, na] * np.log(grid) + (nexc + b - 1)[:, na] * np.log(1 - grid) - \
             (ninc + nexc)[:, na] * np.log(grid * pinc + (1 - grid) * pexc + r)
    logpdf -= logpdf.max(1, keepdims=True)
    pdf = np.exp(logpdf)
    pdf /= pdf.sum(1, keepdims=True)
    pdf = pdf.sum(0) / n_bootstrap_samples

    return pdf, grid

def main(n_bootstrap_samples=1000, n_grid_points=100):
    suffix = "psi_bootstrap"
    file_result = suffix + ".result"
    file_sample = suffix + ".sample"
    header = "#ID\tn_inc\tn_exc\tp_inc\tp_exc\tPSI_standard\tPSI_bootstrap\tPSI_bootstrap_std\n"

    f_r = open(file_result, 'w')
    f_r.write(header) 

    file_junc_rd = 'event_inc_exc.rd'
    data_rd = np.genfromtxt(file_junc_rd, delimiter="\t", dtype="O,O,O", 
                            names=True)
    ID = data_rd['ID']
    n_inc = [np.array(map(int, row.strip(',').split(','))) for row in data_rd['RD_inc']]
    n_exc = [np.array(map(int, row.strip(',').split(','))) for row in data_rd['RD_exc']]    
    N = np.size(n_inc)

    psi_bootstrap = np.zeros(N)
    psi_b_std = np.zeros(N)
    psi_standard = np.zeros(N)

    for i, id, inc, exc in izip(np.arange(N), ID, n_inc, n_exc):
        p_inc = np.size(inc)
        p_exc = np.size(exc)
        min_p = min(p_inc,p_exc)
        scaled_inc = np.sum(inc) / p_inc * min_p
        scaled_exc = np.sum(exc) / p_exc * min_p
        psi_standard[i] = (scaled_inc + 1) / (scaled_inc + scaled_exc + 2)

        # generate bootstrap PDF here!
        pdf, grid = gen_pdf(inc, exc, n_bootstrap_samples, n_grid_points)

        # process the PDF to get PSI
        psi_bootstrap[i] = sum(pdf * grid)
        psi_b_std[i] = np.sqrt(sum(pdf * np.square(grid - psi_bootstrap[i])))
        str_out = ("%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n" % (id, np.sum(inc),
                                                         np.sum(exc),
                                                         p_inc, p_exc,
                                                         psi_standard[i],
                                                         psi_bootstrap[i],
                                                         psi_b_std[i]))
        f_r.write(str_out)                            

        # print progress
        if not (i + 1) % 1000:
            print "%d out of %d AS events have been processed..." % (i + 1, N)

    f_r.close()

    tend = datetime.now()
    tdiff = tend - tstart
    runtime = tdiff.seconds / 60.0
    print "Finished processing all %d AS events in %g minutes!" % (N, runtime)

if __name__ == '__main__':
    main()
