import numpy as np


def beta_gs():
    """
    Calculate momentum transfer coefficient due to drag by gas. See Equations
    8, 9, 10, 11 in Agu 2019 paper.

    Parameters
    ----------

    Returns
    -------

    References
    ----------
    Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre
    Tokheim, and Britt M.E. Moldestad. Detailed One-Dimensional Model for
    Steam-Biomass Gasification in a Bubbling Fluidized Bed. Energy and Fuels,
    vol. 33, pp. 7385-7397, 2019.
    """
    beta = (6 * dr) / (np.pi * ds**3)
    return beta


def beta_ps():
    """
    here
    """
    return 1


def drag_coeff(ds, mug, phi, rhog, u, v):
    """
    Calculate the gas-solid drag coefficient. See Equations 9, 10, 11 in Agu
    2019 paper.

    Parameters
    ----------
    ds : float or array
        Average diameter of a solid fuel particle [m]
    mug : float or array
        Gas dynamic viscosity [kg/m·s]
    phi : float
        Solid fuel particle sphericity [-]
    rhog : float or array
        Gas density [kg/m³]
    u : float or array
        Gas velocity [m/s]
    v : float or array
        Solid fuel particle velocity [m/s]

    Returns
    -------
    cd : float or array
        Gas-solid drag coefficient [-]

    References
    ----------
    Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre
    Tokheim, and Britt M.E. Moldestad. Detailed One-Dimensional Model for
    Steam-Biomass Gasification in a Bubbling Fluidized Bed. Energy and Fuels,
    vol. 33, pp. 7385-7397, 2019.
    """
    re = (rhog * ds / mug) * abs(u + v)
    cd1 = (24 / re) * (1 + 8.1716 * np.exp(-4.0655 * phi)) * (re**(0.0964 + 0.5565 * phi))
    cd2 = (73.69 * re * np.exp(-5.0748 * phi)) / (re + 5.378 * np.exp(6.2122 * phi))
    cd = cd1 + cd2
    return cd


def drag_resist(cd, ds, rhog, u, v):
    """
    Calculate drag resistance of the gas-solid. See Equation 9 in Agu 2019
    paper.

    Parameters
    ----------
    cd : float or array
        Gas-solid drag coefficient [-]
    ds : float or array
        Average diameter of a solid fuel particle [m]
    rhog : float or array
        Gas density [kg/m³]
    u : float or array
        Gas velocity [m/s]
    v : float or array
        Solid fuel particle velocity [m/s]

    Returns
    -------
    dr : float or array
        Drag resistance [N·s/m]

    References
    ----------
    Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre
    Tokheim, and Britt M.E. Moldestad. Detailed One-Dimensional Model for
    Steam-Biomass Gasification in a Bubbling Fluidized Bed. Energy and Fuels,
    vol. 33, pp. 7385-7397, 2019.
    """
    dr = 1 / 8 * np.pi * (ds**2) * rhog * cd * abs(u + v)
    return dr


def ds_average(db, n1, psi, xc, ysc):
    """
    Calculate the average diameter of the solid fuel particle. See Equation 12
    in Agu 2019 paper. Use xc=0 for biomass pyrolysis where there is no char
    gasification.

    Parameters
    ----------
    db : float
        Sauter mean diameter of the biomass particle [m]
    n1 : float
        Fragmentation factor of biomass particles [-]
    psi : float or array
        Shrinkage factor of biomass particle [-]
    xc : float or array
        Char conversion factor [-]
    ysc : float or array
        Mass fraction of char [-]

    Returns
    -------
    ds : float or array
        Average diameter of the solid fuel particle [m]

    References
    ----------
    Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre
    Tokheim, and Britt M.E. Moldestad. Detailed One-Dimensional Model for
    Steam-Biomass Gasification in a Bubbling Fluidized Bed. Energy and Fuels,
    vol. 33, pp. 7385-7397, 2019.
    """
    ds = db / (1 + (1.25 * (n1 * psi * (1 - xc))**(1 / 3) - 1) * ysc)
    return ds


def psi_biomass(rhob, rhoc, wa, wc):
    """
    Calculate biomass particle shrinkage factor.

    Parameters
    ----------
    rhob : float or array
        Biomass density [kg/m³]
    rhoc : float or array
        Char density [kg/m³]
    wa : float
        Weight fraction of ash in biomass particle [-]
    wc : float
        Weight fraction of char in biomass particle [-]

    Returns
    -------
    psi : float or array
        Biomass particle shrinkage factor [-]
    """
    psi = rhoc / (rhob * (wc + wa))
    return psi


def ysc_mass_frac(rhosb, rhosc):
    """
    Calculate mass fraction of char particles in solid fuel mixture. See
    Equation 14 in Agu 2019 paper.

    Parameters
    ----------
    rhosb : float or array
        Mass concentration of biomass in solid fuel [kg/m³]
    rhosc : float or array
        Mass concentration of char in solid fuel [kg/m³]

    Returns
    -------
    ysc : float or array
        Mass fraction of char [-]

    References
    ----------
    Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre
    Tokheim, and Britt M.E. Moldestad. Detailed One-Dimensional Model for
    Steam-Biomass Gasification in a Bubbling Fluidized Bed. Energy and Fuels,
    vol. 33, pp. 7385-7397, 2019.
    """
    ysc = rhosc / (rhosb + rhosc)
    return ysc
