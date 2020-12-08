import numpy as np


def init_cond(params):
    """
    Initial conditions.

    Parameters
    ----------
    params : dict
        Dictionary of values from parameters file

    Returns
    -------
    ic : dict
        Dictionary of initial conditions
    """

    # Parameters
    # ------------------------------------------------------------------------
    din_bio = params['din_bio']
    mf_bio = params['mfin_bio']
    rho_bio = params['rho_bio']

    # Determine initial conditions
    # ------------------------------------------------------------------------

    # cross-section area of biomass inlet [m²]
    ain_bio = (np.pi / 4) * (din_bio**2)

    # initial solid fuel particle velocity [m/s]
    v0 = mf_bio / (rho_bio * ain_bio)

    # Initial conditions
    # ------------------------------------------------------------------------

    ic = {
        'v0': v0
    }

    return ic
