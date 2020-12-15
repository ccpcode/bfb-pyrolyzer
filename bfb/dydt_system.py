import chemics as cm
import numpy as np


def dydt(t, y, dz, z, n, ni, params):
    """
    System of ODEs where
    y[0] = rho_bb (ρ︦𝚋)
    y[1] = rho_cc (ρ︦𝚌)
    y[2] = v

    Parameters
    ----------
    t : float
        Time [s] in func(t, y)
    y : array
        Variables in func(t, y)
    dz : array
        Step in z-direction, Δz [m]
    n : float
        Grid point at top of reactor, total grid points
    ni : float
        Grid point at top of bed
    params : dict
        Parameters from JSON file

    Returns
    -------
    dydt : array
        System of ODEs
    """

    # Variables
    # ------------------------------------------------------------------------
    # rho_bb    biomass mass concentration, ρ︦𝚋 [kg/m³]
    # rho_cc    char mass concentration, ρ︦𝚌 [kg/m³]
    # v         solid fuel velocity, v [m/s]

    rho_bb = y[0:n]
    rho_cc = y[n:2 * n]
    v = y[2 * n:]

    # Constants and parameters
    # ------------------------------------------------------------------------
    # g    acceleration from gravity [m/s²]
    # r    universal gas constant [kJ/(mol K)]
    # ts   solid fuel temperature [K]

    g = 9.81
    r = 0.008314

    # TODO: define temperature in parameters file
    ts = 773.15

    d = params['d']
    d_b = params['d_b']
    d_p = params['d_p']
    e = params['e']
    emf = params['emf']
    theta_w = params['theta_w']
    mf = params['mf']
    nf1 = params["nf1"]
    phi_p = params['phi_p']
    psi_s = params['psi_s']
    mu_c = params['mu_c']
    n1 = params['n1']
    rho_b = params['rho_b']
    rho_c = params['rho_c']
    rho_p = params['rho_p']
    ug = params['ug']
    wa = params['wa']
    wc = params['wc']
    xc = params['xc']

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

    sb = -(k1 + k2 + k3) * rho_bb
    sc = k3 * rho_bb

    # Biomass mass concentration in bed ∂ρ︦𝚋/∂t
    # ------------------------------------------------------------------------
    # ac         inner cross-section area of reactor bed [m²]
    # drhobb_dt  biomass mass concentration ∂ρ︦𝚋/∂t where ρ︦𝚋 is [kg/m³]

    ac = (np.pi / 4) * (d**2)

    # below fuel inlet
    drhobb_dt = np.zeros(n)
    drhobb_dt[0:n1] = - (-rho_bb[1:n1 + 1] * v[1:n1 + 1] + rho_bb[0:n1] * v[0:n1]) / dz[0:n1] + sb[0:n1]

    # at fuel inlet
    drhobb_dt[n1] = - (rho_bb[n1] * v[n1] - rho_bb[n1 - 1] * v[n1 - 1]) / dz[n1] + mf / (ac * dz[n1]) + sb[n1]

    # above fuel inlet in bed
    drhobb_dt[n1 + 1:ni] = - (rho_bb[n1 + 1:ni] * v[n1 + 1:ni] - rho_bb[n1:ni - 1] * v[n1:ni - 1]) / dz[n1 + 1:ni] + sb[n1 + 1:ni]

    # Char mass concentration in bed ∂ρ︦𝚌/∂t
    # ------------------------------------------------------------------------
    # drhocc_dt  char mass concentration ∂ρ︦𝚌/∂t where ρ︦𝚌 is [kg/m³]

    drhocc_dt = np.zeros(n)
    drhocc_dt[0:ni] = - (-rho_cc[1:ni + 1] * v[1:ni + 1] + rho_cc[0:ni] * v[0:ni]) / dz[0:ni] + sc[0:ni]

    # Momentum exchange with inert bed particles Fʙ'
    # ------------------------------------------------------------------------
    # d_bub    bubble diameter, dʙ [m]
    # rho_b    biomass density, ρ𝖻 [kg/m³]
    # rho_c    char density, ρ𝖼 [kg/m³]
    # rho_g    gas density, ρɡ [kg/m³]
    # us       superficial gas velocity at inlet, u𝗌𝚏 [m/s]
    # umf      minimum fluidization velocity, u𝗆𝚏 [m/s]
    # u_bub    bubble velocity, uʙ [m/s]
    # v_bub    bubble volumetric flux, Vʙ [m/s]
    # f_bub    force on fuel particles by inert bed material due to bubble flow, Fʙ' [N/m³]

    # TODO: calculate rho_g, this value is for nitrogen gas
    rho_g = 0.422

    # TODO: define pressure and temp in parameters file
    q_lpm = cm.slm_to_lpm(ug, 101.3, 773.15)
    q_m3s = q_lpm / 60_000
    us = q_m3s / ac

    # TODO: calculate gas properties
    umf = cm.umf_ergun(d_p, emf, 3.6e-5, phi_p, rho_g, rho_p)

    d_bub = 0.00853 * (1 + 27.2 * (us - umf))**(1 / 3) * (1 + 6.84 * z)**(1.21)
    u_bub = 12.51 * ((us - umf)**0.362) * ((d_bub / d)**0.52) * d
    v_bub = 1.285 * ((d_bub / d)**1.52) * d

    f_bub = -(1 - emf) * rho_p * theta_w * v_bub[0:ni] * (-u_bub[1:ni + 1] + u_bub[0:ni]) / dz[0:ni]

    # Momentum exchange with opposite flow βg,𝗌
    # ------------------------------------------------------------------------
    # beta_gs   momentum transfer coefficient due to gas drag, βg,𝗌 [N⋅s/m⁴]
    # cd        here
    # d_s       average diameter of solid fuel particle, d𝗌 [m]
    # dr        here
    # psi       biomass shrinkage factor, Ѱ [-]
    # re_s      here
    # rho_ss    solid fuel mass concentration, ︦ρ𝗌 [kg/m³]
    # ysc       mass fraction of char in solid fuel [-]

    # TODO: calculate gas velocity and gas viscosity
    u = us
    mu_g = 3.6e-5

    rho_ss = rho_cc + rho_bb
    ysc = rho_cc / rho_ss

    psi = rho_c / (rho_b * (wc + wa))
    d_s = d_b / (1 + (1.25 * (nf1 * psi * (1 - xc))**(1 / 3) - 1) * ysc)

    re_s = (rho_g * d_s / mu_g) * np.abs(u + v)
    cd1 = (24 / re_s) * (1 + 8.1716 * np.exp(-4.0655 * psi_s)) * (re_s**(0.0964 + 0.5565 * psi_s))
    cd2 = (73.69 * re_s * np.exp(-5.0748 * psi_s)) / (re_s + 5.378 * np.exp(6.2122 * psi_s))
    cd = cd1 + cd2

    dr = (1 / 8) * np.pi * (d_s**2) * rho_g * cd * np.abs(u + v)

    beta_gs = (6 * dr) / (np.pi * d_s**3)

    # Momentum exchange with inert particles β𝗉,𝗌
    # ------------------------------------------------------------------------
    # alpha_p   solids volume fraction of inert bed material, ⍺𝗉 [-]
    # alpha_s   solids volume fraction of fuel particles, ⍺𝗌 [-]
    # beta_ps   momentum transfer coefficient due to collision with inert bed particles, β𝗉,𝗌 [N⋅s/m⁴]
    # d_bub_avg bubble diameter averaged over the bed height, ḏʙ [m]
    # de        bed expansion, Δe [-]
    # ef        bed voidage at fluidized state, e𝚏 [-]
    # go        radial distribution function [-]
    # rho_s    solid fuel density, ρ𝗌 [kg/m³]

    d_bub_avg = np.mean(d_bub[0:ni])
    de = (1 - 0.103 * ((us - umf)**-0.362) * (d_bub_avg / d))**(-1) - 1
    ef = 1 - (1 - emf) / (1 + de)

    rho_s = ((ysc / rho_c) * ((1 - ysc) / rho_b))**-1
    alpha_s = rho_ss / rho_s
    alpha_p = 1 - ef - alpha_s

    go = (1 / ef) + ((3 * d_s * d_p) / (ef**2 * (d_s + d_p))) * ((alpha_s / d_s) + (alpha_p / d_p))

    beta_ps1 = 3 * np.pi * (1 + e) * (0.5 + (mu_c * np.pi / 8)) * (d_s + d_p)**2
    beta_ps2 = (rho_p * d_p**3) + (rho_s * d_s**3)
    beta_ps = (beta_ps1 / beta_ps2) * alpha_p * rho_p * rho_ss * go * np.abs(v)

    # Solid fuel velocity ∂v/∂t
    # ------------------------------------------------------------------------
    # dv_dt     solid fuel velocity ∂v/∂t where v is [m/s]

    dv_dt = np.zeros(n)

    dv_dt[0:ni] = (
        - v[0:ni] * (v[0:ni] - v[1:ni + 1]) / dz[0:ni]
        + g * (rho_s[0:ni] - rho_g) / rho_s[0:ni]
        + f_bub / rho_s[0:ni]
        # + beta_gs[0:ni] * (-u - v[0:ni]) / rho_s[0:ni]
        + beta_ps[0:ni] * (-v[0:ni]) / rho_s[0:ni]
        + v[0:ni] * (sb[0:ni] + sc[0:ni]) / rho_s[0:ni]
    )

    dydt = np.concatenate((drhobb_dt, drhocc_dt, dv_dt))
    return dydt
