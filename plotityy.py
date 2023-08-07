import argparse
from bokeh.plotting import figure, show, save
from bokeh.io import output_file
from bokeh.io import export_png


def energy_bokeh_plotter(filename, xlabel, ylabel, legend_names, colors):
    # Read data from the file
    with open('energy_data.dat', 'r') as f:
        data = f.readlines()

    time = []
    kinetic_energy = []
    potential_energy = []
    total_energy = []

    for line in data:
        if line.strip() and not line.startswith('#'):
            t, ke, pe, te = line.split()
            time.append(float(t))
            kinetic_energy.append(float(ke))
            potential_energy.append(float(pe))
            total_energy.append(float(te))

    # Create Bokeh plot
    output_file(filename + '.html')
    p = figure(title="Energy Plot", x_axis_label=xlabel, y_axis_label=ylabel,width=1200,height=600)

    p.line(time, kinetic_energy, legend_label=legend_names[0], line_color=colors[0])
    p.line(time, potential_energy, legend_label=legend_names[1], line_color=colors[1])
    p.line(time, total_energy, legend_label=legend_names[2], line_color=colors[2])

    p.legend.location = "top_right"
    p.legend.click_policy = "hide"

    # Show the plot
    show(p)

    # Save the plot as PNG
    export_png(p, filename=filename + '.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot energy data and save as PNG and HTML.')
    parser.add_argument('f1', type=str, help='Filename for saving the plot')
    parser.add_argument('f2', type=str, help='Filename for saving the plot')
    parser.add_argument('f3', type=str, help='Filename for saving the plot')
    parser.add_argument('f4', type=str, help='Filename for saving the plot')
    parser.add_argument('f5', type=str, help='Filename for saving the plot')
    parser.add_argument('f6', type=str, help='Filename for saving the plot')
    parser.add_argument('f7', type=str, help='Filename for saving the plot')
    args = parser.parse_args()

    filename=f"Final-V{args.f1}.{args.f2}_nP-{args.f3}_Bs-{args.f4}_T-{args.f5}_ls-{args.f6}_ns-{args.f7}"
    xlabel = "Time"
    ylabel = "Energy per Particle"
    legend_names = ["Kinetic Energy", "Potential Energy", "Total Energy"]
    colors = ['blue', 'green', 'red']
    energy_bokeh_plotter(filename, xlabel, ylabel, legend_names, colors)
