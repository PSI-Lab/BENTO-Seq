import numpy as np
from  numpy import newaxis as na
from datetime import datetime
from itertools import izip

tstart = datetime.now()
    
def gen_pdf(inc, exc, n_bootstrap_samples=1000, n_grid_points=100, a=1., b=1., r=0.):
    """Generate bootstrap PDF of PSI"""

    grid = np.arange(1. / (2 * n_grid_points), 1., 1. / n_grid_points)
    pinc = inc.size
    pexc = exc.size
    i_inc = np.random.randint(0, pinc, (pinc, n_bootstrap_samples))
    i_exc = np.random.randint(0, pexc, (pexc, n_bootstrap_samples))
    
    ninc = np.sum(inc[i_inc], axis=0)
    nexc = np.sum(exc[i_exc], axis=0)

    logpdf = (ninc + a - 1)[:, na] * np.log(grid) + (nexc + b - 1)[:, na] * np.log(1 - grid) - \
             (ninc + nexc)[:, na] * np.log(grid * pinc + (1 - grid) * pexc + r)
    logpdf -= logpdf.max(1, keepdims=True)
    pdf = np.exp(logpdf)
    pdf /= pdf.sum(1, keepdims=True)
    pdf = pdf.sum(0) / n_bootstrap_samples

    return pdf, grid
