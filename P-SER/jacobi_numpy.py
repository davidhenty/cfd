# 
# Jacobi routine for CFD calculation
#
import numpy as np

import sys

def jacobi(niter, psi):

    (m, n) = psi.shape
    m = m - 2
    n = n - 2

    tmp = np.zeros((m+2, n+2))
    for iter in range(1,niter+1):
        # Use index notation and offsets to compute the stream function
        tmp[1:m+1,1:n+1] = 0.25 * (psi[2:m+2,1:n+1]+psi[0:m,1:n+1]+psi[1:m+1,2:n+2]+psi[1:m+1,0:n])

        # Update psi
        np.copyto(psi[1:m+1,1:n+1], tmp[1:m+1,1:n+1])

        if iter%1000 == 0:
            sys.stdout.write("completed iteration {0}\n".format(iter))
