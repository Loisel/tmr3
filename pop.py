"""
Module to calculate the populations from the rates.
Different algorithms were implemented, but only
the asymptotic method works reliably.
"""

import scipy as sp
import cfg

from scipy import linalg,matrix


def asymptotic_ssp(ME):
    """
    A asymptotic evolution of the density matrix to calculate
    the steady state populations.
    This proved to work fast and reliably.
    The timescale is chosen as a multiple of the coupling over
    hbar.

    Args:
      ME: the density matrix of size N,N
    Returns:
      The population vector of size N
    """

    # the time scale. This has to be a large number > 10**5
    tscale = 10.*cfg.conf["G_scale"]/cfg.conf["HBAR"]

    # initial guess for the populations: equally populated states
    P0 = sp.ones(ME.shape[0])*1./ME.shape[0]

    # the calculation itself is just a exponential matrix multiplication
    # provided by scipy
    result =  sp.dot(linalg.expm(ME*tscale),P0)

    # uncomment to check the conservation of the total probability
    # print "Probability sum: {}".format(sp.sum(result))

    return result

    


