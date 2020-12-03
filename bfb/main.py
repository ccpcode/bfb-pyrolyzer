import plotter
from grid_points import grid_points
from initial_conditions import init_cond
from read_json import read_json
from solver import solver


def main():
    """
    Run BFB biomass pyrolyzer 1-D model.
    """

    # Get parameters
    params = read_json('params.json')

    # Configuration
    t = params['nt'] * params['dt']
    print(f'Total time = {t} s')

    # Grid points and distance between all grid points
    z, dz = grid_points(params)
    nz = len(z)

    # Initial conditions
    ic = init_cond(params)

    # Run solver
    results = solver(params, ic, dz, nz)

    # Plot results
    plotter.plot_solid_velocity(results['v'], z)
    plotter.show_plots()


if __name__ == '__main__':
    main()
