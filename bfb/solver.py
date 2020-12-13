import chemics as cm
import logging
import numpy as np
from scipy.integrate import solve_ivp
from dydt_system import dydt


def solver(params):
    """
    Solve one-dimensional biomass pyrolysis model.
    """

    # One-dimensional grid
    # ------------------------------------------------------------------------

    ls = params['ls']       # distance to solid fuel inlet [m]
    lst = params['lst']     # static bed height [m]
    lt = params['lt']       # reactor height [m]
    lp = lst - ls           # distance from solid fuel inlet to bed top [m]

    n1 = params['n1']       # number of grid points below fuel inlet
    n2 = params['n2']       # number of grid points above fuel inlet
    n3 = params['n3']       # number of grid points in freeboard
    ni = n1 + n2            # grid point at top of bed
    n = n1 + n2 + n3        # grid point at top of reactor

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

    # Bubble phase
    # ------------------------------------------------------------------------

    di = params['di']       # inner diameter of the reactor [m]
    dpp = params['dpp']     # mean diameter of bed particle [m]
    emf = params['emf']     # void fraction of the bed [-]
    phip = params['phip']   # sphericity of bed particle [-]
    rhop = params['rhop']   # bed particle density [kg/m³]
    ug = params['ug']       # velocity of inlet gas [SLM]

    ac = (np.pi / 4) * (di**2)  # inner cross-section area of the reactor [m²]

    # TODO: calculate rhog, this value is for nitrogen gas
    rhog = 0.422

    # TODO: define pressure and temp in parameters file
    # superficial gas velocity at gas inlet [m/s]
    q_lpm = cm.slm_to_lpm(ug, 101.3, 773.15)
    q_m3s = q_lpm / 60_000
    us = q_m3s / ac

    # TODO: calculate gas properties
    # minimum fluidization velocity [m/s]
    umf = cm.umf_ergun(dpp, emf, 3.6e-5, phip, rhog, rhop)

    db = 0.00853 * (1 + 27.2 * (us - umf))**(1 / 3) * (1 + 6.84 * z)**(1.21)
    vb = 1.285 * ((db / di)**1.52) * di
    ub = 12.51 * ((us - umf)**0.362) * ((db / di)**0.52) * di

    # Initial conditions
    # ------------------------------------------------------------------------

    dis = params['dis']     # inner diameter of biomass inlet tube [m]
    mf = params['mf']       # biomass inlet feed flow rate [kg/s]
    rhob = params['rhob']   # biomass density [kg/m³]

    ai = (np.pi / 4) * (dis**2)     # inner cross-section area of inlet [m²]
    vi = mf / (rhob * ai)           # biomass velocity at fuel inlet [m/s]

    rhobb0 = np.zeros(n)    # initial biomass mass concentration [kg/m³]
    rhobb0[:] = 1e-12

    rhocc0 = np.zeros(n)    # initial char mass concentration [kg/m³]

    v0 = np.zeros(n)        # initial solid fuel velocity [m/s]
    v0[n1] = vi

    # Solve ODEs
    # ------------------------------------------------------------------------

    tf = params['tf']   # final time [s]
    tspan = (0, tf)     # time span [s]

    y0 = np.concatenate((rhobb0, rhocc0, v0))   # initial conditions vector

    sol = solve_ivp(dydt, tspan, y0, method='Radau', args=(ac, dz, n, ni, ub, vb, params))

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
