import argparse
import os
import plotly.graph_objects as go
import plotly.io as pio
from time import sleep

def energy_plotly_plotter(filenames, xlabel, ylabel, legend_names, colors, flg=0):
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

    # Create Plotly figure
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=time, y=kinetic_energy, mode='lines', name=legend_names[0], line=dict(color=colors[0])))
    fig.add_trace(go.Scatter(x=time, y=potential_energy, mode='lines', name=legend_names[1], line=dict(color=colors[1])))
    fig.add_trace(go.Scatter(x=time, y=total_energy, mode='lines', name=legend_names[2], line=dict(color=colors[2])))

    fig.update_layout(title="Energy Plot", xaxis_title=xlabel, yaxis_title=ylabel, width=1200, height=600)

    # Save or show the plot
    if flg == 0:
        fig.show()
        pio.write_image(fig, filenames + '.png')
        pio.write_html(fig, filenames + '.html')
    else:
        pio.write_html(fig, filenames + '.html')
        if flg <2:
            os.system(f"start {os.getcwd()}/tempi.html")

    # Save the plot as PNG
    

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
    # print("main")

    filename=f"M6-V{args.f1}.{args.f2}_nP-{args.f3}_Bs-{args.f4}_T-{args.f5}_ls-{args.f6}_ns-{args.f7}"
    xlabel = "Time"
    ylabel = "Energy per Particle"
    legend_names = ["Kinetic Energy", "Potential Energy", "Total Energy"]
    colors = ['blue', 'green', 'red']
    energy_plotly_plotter(filename, xlabel, ylabel, legend_names, colors)

else :
    co=1
    filename="tempi"
    xlabel = "Time"
    ylabel = "Energy per Particle"
    legend_names = ["Kinetic Energy", "Potential Energy", "Total Energy"]
    colors = ['blue', 'green', 'red']
    # print("not main")
    while True:
        energy_plotly_plotter(filename, xlabel, ylabel, legend_names, colors,co)
        co+=1
        sleep(10) #change to increase delay
        