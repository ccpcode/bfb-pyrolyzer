import numpy as np


def solver(params, initcond, dz, nz):
    """
    Solve equations for 1-D model of BFB biomass pyrolyzer.

    Parameters
    ----------
    params : dict
        Dictionary of values from parameters file
    initcond : dict
        Dictionary of initial conditions
    dz : array
        Distance between grid points, ∆z [m]
    nz : float
        Total number of grid points in z-direction [-]

    Returns
    -------
    results : dict
        Dictionary of results
    """

    # Parameters
    # ------------------------------------------------------------------------
    n0 = params['n0']
    n1 = params['n1']
    nt = params['nt']
    dt = params['dt']

    # Constants
    # ------------------------------------------------------------------------
    g = 9.81    # acceleration from gravity [m/s²]

    # Initial conditions
    # ------------------------------------------------------------------------

    # TODO: calculate rhos
    rhos = np.ones(nz - 1)
    rhos[0:n0 + n1] = 540
    rhos[n0 + n1:] = 1e-12

    # TODO: calculate rhog
    rhog = np.ones(nz - 1)
    rhog[0:n0 + n1] = 0.422
    rhog[n0 + n1:] = 1e-12

    # Solid phase momentum
    # ------------------------------------------------------------------------

    # Initialize array for solid fuel particle velocity [m/s]
    # rows = time steps ∆t
    # columns = grid points along reactor height in z-direction
    v = np.zeros((nt, nz))

    # Assign initial conditions
    # v[0, 1:n0 + n1] = initcond['v0']
    v[0, n0] = initcond['v0']

    # Calculate v for each time step
    for n in range(1, nt):
        vn = v[n - 1]

        v[n, 1:] = (
            vn[1:]
            - vn[1:] * (vn[1:] - vn[:-1]) * dt / dz
            + g * (rhos - rhog) * dt / rhos
        )

    # Results
    # ------------------------------------------------------------------------

    results = {
        'v': v
    }

    return results
