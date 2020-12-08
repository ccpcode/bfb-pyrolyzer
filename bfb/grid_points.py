import matplotlib.pyplot as plt
import numpy as np
from read_json import read_json


def grid_points(params):
    """
    Create 1-D grid and determine distance between grid points.

    Parameters
    ----------
    params : dict
        Dictionary of values from parameters file

    Returns
    -------
    z : array
        Grid points along z-axis [m]
    dz : array
        Distance between grid points, ∆z [m]
    """

    # Parameters (n0, n1, n2, zin, zst, zt)
    n0 = params['n0']
    n1 = params['n1']
    n2 = params['n2']
    zin = params['zin']
    zst = params['zst']
    zt = params['zt']

    # Grid points
    z0 = np.linspace(0, zin, n0)
    z1 = np.linspace(zin, zst, n1)
    z2 = np.linspace(zst, zt, n2)
    z = np.concatenate((z0, z1[1:], z2[1:]))

    # Distance ∆z between all grid points
    dz = np.diff(z)

    return z, dz


if __name__ == '__main__':
    """
    Examples of creating 1-D grid and determining the distance between all the
    grid points.
    """

    # Example 1
    # -------------------------------------------

    # Parameters
    n0 = 7
    n1 = 5
    n2 = 8
    zin = 2
    zst = 5
    zt = 14

    # Distance ∆z between grid points for each section
    # Check with dz output from function
    dz0 = zin / (n0 - 1)
    dz1 = (zst - zin) / (n1 - 1)
    dz2 = (zt - zst) / (n2 - 1)

    # Grid points and distance between grid points ∆z
    params = {'n0': n0, 'n1': n1, 'n2': n2, 'zin': zin, 'zst': zst, 'zt': zt}
    z, dz = grid_points(params)

    # Print
    print('dz0 \t', dz0)
    print('dz1 \t', dz1)
    print('dz2 \t', dz2)
    print('nt \t', n0 + n1 + n2)
    print('len z \t', len(z))
    print('len dz \t', len(dz))
    print('\nz\n', z)
    print('\ndz\n', dz)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.8), tight_layout=True)
    ax1.plot(dz, marker='.')
    ax1.set_xlabel('Grid point delta [-]')
    ax1.set_ylabel('Distance between grid points, ∆z [m]')
    ax2.plot(z, marker='.')
    ax2.axvline(n0 - 1, color='C1', label='n0')
    ax2.axvline(n0 + n1 - 2, color='C2', label='n0 + n1')
    ax2.set_xlabel('Grid point, n [-]')
    ax2.set_ylabel('Reactor height, z [m]')
    ax2.legend()

    # Example 2
    # -------------------------------------------

    params = read_json('params.json')
    z, dz = grid_points(params)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.8), tight_layout=True)
    ax1.plot(dz, marker='.')
    ax1.set_xlabel('Grid point delta [-]')
    ax1.set_ylabel('Distance between grid points, ∆z [m]')
    ax2.plot(z, marker='.')
    ax2.set_xlabel('Grid point, n [-]')
    ax2.set_ylabel('Reactor height, z [m]')

    plt.show()
