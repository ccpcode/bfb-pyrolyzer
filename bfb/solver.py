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
    din_bio = params['din_bio']
    mfin_bio = params['mfin_bio']
    rhob = params['rho_bio']
    rhoc = params['rho_char']

    # Constants
    # ------------------------------------------------------------------------
    g = 9.81        # acceleration from gravity [m/s²]
    r = 0.008314    # universal gas constant [kJ/(mol K)]

    # TODO: calculate rhog, this value is for nitrogen gas
    rhog = 0.422

    # TODO: calculate ts
    ts = 773.15     # solid fuel temperature [K]

    a1 = 1.3e8      # biomass -> volatiles frequency factor [1/s]
    a2 = 2.0e8      # biomass -> tar frequency factor [1/s]
    a3 = 1.08e7     # biomass -> char frequency factor [1/s]

    e1 = 140        # biomass -> volatiles activation energy [kJ/mol]
    e2 = 133        # biomass -> tar activation energy [kJ/mol]
    e3 = 121        # biomass -> char activation energy [kJ/mol]

    # Setup
    # ------------------------------------------------------------------------

    # cross-section area of biomass feed inlet [m²]
    ain_bio = (np.pi / 4) * (din_bio**2)

    # number of grid points to top of bed
    nbed = n0 + n1

    # Initialize arrays
    # ------------------------------------------------------------------------
    # rows = time steps ∆t
    # columns = grid points along reactor height in z-direction

    # arrays for biomass and char mass concentrations [kg/m³]
    # to prevent division by zero use 1e-12 instead of 0
    rhobb = np.ones((nt, nz)) * 1e-12
    rhocc = np.ones((nt, nz)) * 1e-12

    # array for solid fuel velocity [m/s]
    v = np.zeros((nt, nz))

    # Initial conditions
    # ------------------------------------------------------------------------

    # initial velocity of solid fuel [m/s]
    # v[0, 1:n0 + n1] = initcond['v0']
    v[0, n0] = initcond['v0']

    # Solve
    # ------------------------------------------------------------------------

    for n in range(1, nt):

        # state of the system at time step n
        rhobbn = rhobb[n - 1]
        rhoccn = rhocc[n - 1]
        vn = v[n - 1]

        # kinetic rate constants
        k1 = a1 * np.exp(-e1 / (r * ts))
        k2 = a2 * np.exp(-e2 / (r * ts))
        k3 = a3 * np.exp(-e3 / (r * ts))

        # mass generation terms
        sb0 = -(k1 + k2 + k3) * rhobbn[1:n0] * dt
        sb1 = -(k1 + k2 + k3) * rhobbn[n0 + 1:nbed] * dt
        sb = -(k1 + k2 + k3) * rhobbn[n0] * dt
        sc = k3 * rhobbn[1:n0 + n1] * dt

        # biomass mass concentration below inlet (from 0 to n0)
        rhobb[n, 1:n0] = (
            rhobbn[1:n0]
            - (vn[1:n0] * rhobbn[1:n0] - vn[0:n0 - 1] * rhobbn[0:n0 - 1]) * dt / dz[1:n0]
            + sb0
        )

        # biomass mass concentration at inlet (n0)
        rhobb[n, n0] = (
            rhobbn[n0]
            - (vn[n0] * rhobbn[n0] - vn[n0 - 1] * rhobbn[n0 - 1]) * dt / dz[n0]
            + mfin_bio / (ain_bio * dz[n0])
            + sb
        )

        # biomass mass concentration above inlet (from n0 to n1)
        rhobb[n, n0 + 1:nbed] = (
            rhobbn[n0 + 1:nbed]
            - (vn[n0 + 1:nbed] * rhobbn[n0 + 1:nbed] - vn[n0:nbed - 1] * rhobbn[n0:nbed - 1]) * dt / dz[n0 + 1:nbed]
            + sb1
        )

        # char mass concentration
        rhocc[n, 1:nbed] = (
            rhoccn[1:nbed]
            - (vn[1:nbed] * rhoccn[1:nbed] - vn[0:nbed - 1] * rhoccn[0:nbed - 1]) * dt / dz[1:nbed]
            + sc
        )

        # solid fuel velocity
        yc = rhoccn[1:] / (rhoccn[1:] + rhobbn[1:])
        rhos = ((yc / rhoc) + (1 - yc) / rhob)**(-1)

        v[n, 1:] = (
            vn[1:]
            - vn[1:] * (vn[1:] - vn[:-1]) * dt / dz
            + g * (rhos - rhog) * dt / rhos
        )

    # Results
    # ------------------------------------------------------------------------

    results = {
        'rhobb': rhobb,
        'rhocc': rhocc,
        'v': v
    }

    return results
