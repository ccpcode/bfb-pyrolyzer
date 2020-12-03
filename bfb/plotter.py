import matplotlib.pyplot as plt


def plot_solid_velocity(v, z):
    """
    here
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
