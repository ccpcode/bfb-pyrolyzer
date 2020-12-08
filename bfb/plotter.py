import matplotlib.pyplot as plt


def plot_rhobb(rhobb, z):
    """
    Plot biomass mass concentration.
    """
    fig, ax = plt.subplots(tight_layout=True)
    ax.plot(z, rhobb[0])
    ax.plot(z, rhobb[-1])
    ax.set_xlabel('Reactor height, z [m]')
    ax.set_ylabel('Biomass mass concentration, rhobb [kg/m³]')
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


def plot_rhocc(rhocc, z):
    """
    Plot char mass concentration.
    """
    fig, ax = plt.subplots(tight_layout=True)
    ax.plot(z, rhocc[0])
    ax.plot(z, rhocc[-1])
    ax.set_xlabel('Reactor height, z [m]')
    ax.set_ylabel('Char mass concentration, rhocc [kg/m³]')
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


def plot_v(v, z):
    """
    Plot solid fuel velocity.
    """
    fig, ax = plt.subplots(tight_layout=True)
    ax.plot(z, v[0])
    ax.plot(z, v[-1])
    ax.set_xlabel('Reactor height, z [m]')
    ax.set_ylabel('Solid fuel velocity, v [m/s]')
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


def show_plots():
    plt.show()
