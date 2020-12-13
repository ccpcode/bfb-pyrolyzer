import matplotlib.pyplot as plt


def _style_axis(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


def plot_rhobb(params, results):
    """
    Plot biomass mass concentration.
    """
    n1 = params['n1']

    n = results['n']
    rhobb = results['rhobb']
    t = results['t']
    z = results['z']

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 4.8), nrows=1, ncols=2, tight_layout=True)
    ax1.plot(t, rhobb[0], label='z @ 0')
    ax1.plot(t, rhobb[n1 - 1], label='z @ n1')
    ax1.plot(t, rhobb[n - 1], label='z @ n')
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('rhobb [kg/m³]')
    ax1.legend()
    _style_axis(ax1)

    ax2.plot(z, rhobb[:, 0], label='t0')
    ax2.plot(z, rhobb[:, -1], label='tf')
    ax2.set_xlabel('reactor height [m]')
    ax2.set_ylabel('rhobb [kg/m³]')
    ax2.legend()
    _style_axis(ax2)


def plot_rhocc(params, results):
    """
    Plot char mass concentration.
    """
    n1 = params['n1']

    n = results['n']
    rhocc = results['rhocc']
    t = results['t']
    z = results['z']

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 4.8), nrows=1, ncols=2, tight_layout=True)
    ax1.plot(t, rhocc[0], label='z @ 0')
    ax1.plot(t, rhocc[n1 - 1], label='z @ n1')
    ax1.plot(t, rhocc[n - 1], label='z @ n')
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('rhocc [kg/m³]')
    ax1.legend()
    _style_axis(ax1)

    ax2.plot(z, rhocc[:, 0], label='t0')
    ax2.plot(z, rhocc[:, -1], label='tf')
    ax2.set_xlabel('reactor height [m]')
    ax2.set_ylabel('rhocc [kg/m³]')
    ax2.legend()
    _style_axis(ax2)


def plot_v(params, results):
    """
    Plot solid fuel velocity.
    """
    n1 = params['n1']

    n = results['n']
    t = results['t']
    v = results['v']
    z = results['z']

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 4.8), nrows=1, ncols=2, tight_layout=True)
    ax1.plot(t, v[0], label='z @ 0')
    ax1.plot(t, v[n1 - 1], label='z @ n1')
    ax1.plot(t, v[n - 1], label='z @ n')
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('v [m/s]')
    ax1.legend()
    _style_axis(ax1)

    ax2.plot(z, v[:, 0], label='t0')
    ax2.plot(z, v[:, -1], label='tf')
    ax2.set_xlabel('reactor height [m]')
    ax2.set_ylabel('v [m/s]')
    ax2.legend()
    _style_axis(ax2)


def show_plots():
    plt.show()
