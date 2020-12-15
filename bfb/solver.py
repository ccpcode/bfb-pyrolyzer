import logging
import numpy as np
from scipy.integrate import solve_ivp
from dydt_system import dydt


def solver(params):
    """
    Solve one-dimensional biomass pyrolysis model.
    """

    # Parameters
    # ------------------------------------------------------------------------

    tf = params['tf']       # final time [s]

    ls = params['ls']       # distance to solid fuel inlet [m]
    lst = params['lst']     # static bed height [m]
    lt = params['lt']       # reactor height [m]

    n1 = params['n1']       # number of grid points below fuel inlet
    n2 = params['n2']       # number of grid points above fuel inlet
    n3 = params['n3']       # number of grid points in freeboard

    d_si = params['d_si']       # inner diameter of solid fuel inlet tube [m]
    mf = params['mf']           # biomass inlet feed flow rate [kg/s]
    rho_b = params['rho_b']     # biomass density [kg/m³]

    # One-dimensional grid
    # ------------------------------------------------------------------------

    ni = n1 + n2            # grid point at top of bed
    n = n1 + n2 + n3        # grid point at top of reactor

    lp = lst - ls           # distance from solid fuel inlet to bed top [m]
    dz1 = ls / n1           # dz in bed below fuel inlet [m]
    dz2 = lp / n2           # dz in bed above fuel inlet [m]
    dz3 = (lt - lst) / n3   # dz for freeboard [m]

    dz = np.zeros(n)        # dz vector for reactor [m]
    dz[0:n1] = dz1
    dz[n1:ni] = dz2
    dz[ni:n] = dz3

    z1 = np.linspace(0, ls, n1)         # z below fuel inlet [m]
    z2 = np.linspace(ls, lp, n2)        # z above fuel inlet [m]
    z3 = np.linspace(lp, lt, n3)        # z for freeboard [m]
    z = np.concatenate((z1, z2, z3))    # z vector for reactor [m]

    # Initial conditions
    # ------------------------------------------------------------------------

    ai = (np.pi / 4) * (d_si**2)    # inner cross-section area of solid fuel inlet tube [m²]
    vi = mf / (rho_b * ai)           # biomass velocity at fuel inlet [m/s]

    rhobb0 = np.zeros(n)    # initial biomass mass concentration [kg/m³]
    rhobb0[:] = 1e-12

    rhocc0 = np.zeros(n)    # initial char mass concentration [kg/m³]
    rhocc0[:] = 1e-12

    v0 = np.zeros(n)        # initial solid fuel velocity [m/s]
    v0[n1] = vi

    # Solve ODEs
    # ------------------------------------------------------------------------

    tspan = (0, tf)                             # time span [s]
    y0 = np.concatenate((rhobb0, rhocc0, v0))   # initial conditions vector

    sol = solve_ivp(dydt, tspan, y0, method='Radau', args=(dz, z, n, ni, params))

    # Log
    # ------------------------------------------------------------------------

    logging.basicConfig(format='%(message)s', level=logging.INFO)

    log = (
        f't0      {sol.t[0]}\n'
        f'tf      {sol.t[-1]}\n'
        f'n       {n}\n'
        f'len t   {len(sol.t)}\n'
        f'y shape {sol.y.shape}'
    )

    logging.info(log)

    # Results
    # ------------------------------------------------------------------------

    results = {
        'n': n,
        'z': z,
        't': sol.t,
        'rhobb': sol.y[0:n],
        'rhocc': sol.y[n:2 * n],
        'v': sol.y[2 * n:]
    }

    return results
