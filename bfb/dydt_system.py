import numpy as np


def dydt(t, y, ac, dz, n, ni, ub, vb, params):
    """
    System of ODEs where
    y[0] = rhobb
    y[1] = rhocc
    y[2] = v
    """

    # Variables
    # ------------------------------------------------------------------------

    rhobb = y[0:n]
    rhocc = y[n:2 * n]
    v = y[2 * n:]

    # Constants and parameters
    # ------------------------------------------------------------------------

    g = 9.81        # acceleration from gravity [m/s²]
    r = 0.008314    # universal gas constant [kJ/(mol K)]

    rhog = 0.422    # nitrogen gas density [kg/m³]
    ts = 773.15     # solid fuel temperature [K]

    emf = params['emf']
    fw = params['fw']
    mf = params['mf']
    n1 = params['n1']
    rhob = params['rhob']
    rhoc = params['rhoc']
    rhop = params['rhop']

    # Biomass pyrolysis kinetics
    # ------------------------------------------------------------------------

    a1 = 1.3e8      # biomass -> volatiles frequency factor [1/s]
    a2 = 2.0e8      # biomass -> tar frequency factor [1/s]
    a3 = 1.08e7     # biomass -> char frequency factor [1/s]
    e1 = 140        # biomass -> volatiles activation energy [kJ/mol]
    e2 = 133        # biomass -> tar activation energy [kJ/mol]
    e3 = 121        # biomass -> char activation energy [kJ/mol]

    k1 = a1 * np.exp(-e1 / (r * ts))
    k2 = a2 * np.exp(-e2 / (r * ts))
    k3 = a3 * np.exp(-e3 / (r * ts))

    sb = -(k1 + k2 + k3) * rhobb
    sc = k3 * rhobb

    # Biomass mass concentration [kg/m³]
    # ------------------------------------------------------------------------

    # below fuel inlet
    drhobbdt = np.zeros(n)
    drhobbdt[0:n1] = - (-rhobb[1:n1 + 1] * v[1:n1 + 1] + rhobb[0:n1] * v[0:n1]) / dz[0:n1] + sb[0:n1]

    # at fuel inlet
    drhobbdt[n1] = - (rhobb[n1] * v[n1] - rhobb[n1 - 1] * v[n1 - 1]) / dz[n1] + mf / (ac * dz[n1]) + sb[n1]

    # above fuel inlet
    drhobbdt[n1 + 1:ni] = - (rhobb[n1 + 1:ni] * v[n1 + 1:ni] - rhobb[n1:ni - 1] * v[n1:ni - 1]) / dz[n1 + 1:ni] + sb[n1 + 1:ni]

    # Char mass concentration [kg/m³]
    # ------------------------------------------------------------------------

    # below bed top
    drhoccdt = np.zeros(n)
    drhoccdt[0:ni] = - (-rhocc[1:ni + 1] * v[1:ni + 1] + rhocc[0:ni] * v[0:ni]) / dz[0:ni] + sc[0:ni]

    # Solid fuel velocity [m/s]
    # ------------------------------------------------------------------------

    # char mass fraction and density of solid fuel particle
    yc = rhocc / (rhocc + rhobb)
    rhos = ((yc / rhoc) + (1 - yc) / rhob)**(-1)

    # below bed top
    dvdt = np.zeros(n)
    dvdt[0:ni] = (
        - v[0:ni] * (v[0:ni] - v[1:ni + 1]) / dz[0:ni]
        + g * (rhos[0:ni] - rhog) / rhos[0:ni]
        - (1 - emf) * rhop * fw * vb[0:ni] * (-ub[1:ni + 1] + ub[0:ni]) / (rhos[0:ni] * dz[0:ni])
        + v[0:ni] * (sb[0:ni] + sc[0:ni]) / rhos[0:ni]
    )

    dydt = np.concatenate((drhobbdt, drhoccdt, dvdt))
    return dydt
