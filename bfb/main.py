import plotter
from read_json import read_json
from solver import solver


def main():
    """
    Run one-dimensional BFB biomass pyrolysis model.
    """

    # Get parameters
    params = read_json('params.json')

    # Run solver
    results = solver(params)

    # Plot results
    plotter.plot_rhobb(params, results)
    plotter.plot_rhocc(params, results)
    plotter.plot_v(params, results)
    plotter.show_plots()


if __name__ == '__main__':
    main()
