import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import math
import telegram
import asyncio
import argparse
import json
import csv
import os
import re
import sys
import subprocess
from telegram import InputMediaPhoto, InputMediaVideo, InputFile
from matplotlib.animation import FuncAnimation

bot = telegram.Bot(token='6563925351:AAFEgghTh_Z5KKWANcKk76QI_2P-eCgUL-Q')
plot_file_path = ''  # Replace with the actual path
chat_id = '900396766'  # Replace with your chat ID


    
# Load the JSON file
with open('particle_positions.json', 'r') as file:
    data = json.load(file)

# Convert the data to a NumPy array
particle_positions = np.array(data)

def save_param(cpp_file_path="M6.cpp"):
    # Define the regex pattern to match C++ variable declarations and their values
    variable_pattern = r'\b(?:int|double|string)\s+(\w+)\s*=\s*([^;]+);'
    # Initialize a flag to check if we are within the desired block
    inside_block = False
    # Initialize an empty dictionary to store variable names and values
    variable_dict = {}
            
    print(cpp_file_path)
    # Read the file line by line and process it
    with open(cpp_file_path, 'r') as file:
        for line in file:
            # Check if we are entering the block
            if "//params-begin" in line:
                inside_block = True
                continue

            # Check if we are exiting the block
            if "//params-end" in line:
                inside_block = False
                break

            # If we are inside the block, try to match variable declarations and values
            if inside_block:
                match = re.match(variable_pattern, line.strip())
                if match:
                    variable_name = match.group(1)
                    variable_value = match.group(2)
                    variable_dict[variable_name] = variable_value

    # Print the extracted variable names and values
    with open("params.csv","w")as f:
        for variable_name, variable_value in variable_dict.items():
            f.write(f"{variable_name},{variable_value}\n")
    return variable_dict


def energy_plotly_plotter(filenames):
    # Read data from the file
    time = []
    kinetic_energy = []
    potential_energy = []
    total_energy = []
    
    xlabel = "Time"
    ylabel = "Energy per Particle"
    legend_names = ["Kinetic Energy", "Potential Energy", "Total Energy"]
    colors = ['blue', 'green', 'red']
    

    for t in range(len(particle_positions)):
        time.append(float(t)*dataCompression*timestep)
        kinetic_energy_total=0.0
        i=0
        for x, y ,vx ,vy , phi, PE in particle_positions[t]:
            kinetic_energy_total+=(float(0.5 * (vx**2 + vy**2)))
            if i==len(particle_positions[0])-1:
                potential_energy_per_particle=float(PE)
                kinetic_energy_per_particle=kinetic_energy_total/float(len(particle_positions[0]))
            i+=1
        potential_energy.append(potential_energy_per_particle)
        kinetic_energy.append(kinetic_energy_per_particle)
        total_energy.append(potential_energy_per_particle+kinetic_energy_per_particle)

    # Create Plotly figure
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=time, y=kinetic_energy, mode='lines', name=legend_names[0], line=dict(color=colors[0])))
    fig.add_trace(go.Scatter(x=time, y=potential_energy, mode='lines', name=legend_names[1], line=dict(color=colors[1])))
    fig.add_trace(go.Scatter(x=time, y=total_energy, mode='lines', name=legend_names[2], line=dict(color=colors[2])))

    fig.update_layout(title="Energy Plot", xaxis_title=xlabel, yaxis_title=ylabel, width=1200, height=600)

    # Save or show the plot
    pio.write_image(fig, filenames + '-Energy.png')
    pio.write_html(fig, filenames + '-Energy.html')
    return file

def persistence():

    # Create an empty list to store the persistence list of each particle at each gap to take average later
    persistencePerParticle = [[0.0] * int(math.floor(len(data) - 1.0)) for i in range(numParticles)]

    # Create a for loop running through each particle
    for i in range(numParticles):
        # Create an empty list to store the persistence values of each gap
        persistencePerGap = [0.0] * math.floor(len(data) -1.0)

        # Create a for loop running through half the size of the number of timesteps; here j will be called "gap"
        for j in range(math.floor(len(data) -1.0)):
            for k in range(len(data) - j - 1):
                # For each k calculate (cos(data[k][i][4])*cos(data[k+j+1][i][4]))/(cos(data[k][i][4]*cos(data[k][i][4]) and add it to the list
                # cos function is being called with an argument of type float which is a double in C++
                persistencePerGap[j] += (math.cos(float(data[k][i][4])) * math.cos(float(data[k+j+1][i][4])))

            # Take the average of the list and store it in the list
            persistencePerGap[j] /= (len(data) - j - 1)

        # Add the list to the list of lists
        persistencePerParticle[i] = persistencePerGap

    # Create a list which stores the average persistence of each particle at a given gap
    averagePersistencePerGap = [0.0] * math.floor(len(data) -1.0)
    # Create a list to store the standard deviation of persistence values at each gap
    stdPersistencePerGap = [0.0] * math.floor(len(data) -1.0)

    # Create a for loop running through half the size of the number of timesteps; here j will be called "gap"
    for j in range(math.floor(len(data)-1.0)):
        # Create a list to store the persistence values of each particle at a given gap
        persistenceValues = []

        # Create a for loop running through each particle
        for i in range(numParticles):
            # Add the persistence of each particle at a given gap to the list
            persistenceValues.append(persistencePerParticle[i][j])

        # Calculate the mean and standard deviation of the persistence values
        meanPersistence = np.mean(persistenceValues)
        stdPersistence = np.std(persistenceValues)

        # Store the mean and standard deviation in the respective lists
        averagePersistencePerGap[j] = meanPersistence
        stdPersistencePerGap[j] = stdPersistence


# Plot persistence.dat using matplotlib with x axis as i*dataCompression*timestep and y axis as averagePersistencePerGap[i]
    plt.errorbar([i*dataCompression*timestep for i in range(math.floor(len(data)-1.0))], averagePersistencePerGap, yerr=stdPersistencePerGap, ecolor = 'lightblue', capsize=0)
    #reduce font size of x axis and y axis
    
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    
    plt.grid()
    plt.xlabel("Time Gap")
    plt.ylabel("Persistence")
    plt.title("VC=" + theta)
    plt.savefig("persistence.png")

# Function to update the animation
def update(frame):
    ax.clear()
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)

    # Plot the particles at the given frame as circles
    for x, y , vx, vy, phi, PE in particle_positions[frame]:
        circle = plt.Circle((x, y), radius=r, linewidth=0)
        ax.add_patch(circle)
        
    # Add an arrow to represent the direction of movement
        arrow = plt.Arrow(x, y, vx, vy, width=0.1, color='red')
        ax.add_patch(arrow)

    # Add an arrow to represent the direction of movement
        arrow = plt.Arrow(x, y, np.cos(phi), np.sin(phi), width=0.1, color='green')
        ax.add_patch(arrow)

async def send_file(plot_file_path, param_file_path, gif_file_path, persistence_file_path, chat_id):

    try:
        with open(param_file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            # send plot_file_path with the message
            message = "\n".join([', '.join(row) for row in reader])  # Combine CSV rows into a single message with line breaks
            await bot.send_message(chat_id=chat_id, text=message)
    except Exception as e:
        print(f"Error: {e}")
        
    await bot.send_photo(chat_id=chat_id, photo=open(plot_file_path, 'rb'))
    await bot.send_photo(chat_id=chat_id, photo=open(persistence_file_path, 'rb'))
    # await bot.send_animation(chat_id=chat_id, animation=open(gif_file_path, 'rb'))

    

async def main(plot_file_path, param_file_path, gif_file_path, persistence_file_path, chat_id):
    await send_file(plot_file_path, param_file_path, gif_file_path, persistence_file_path, chat_id)


if __name__ == "__main__":
    variable_dict=save_param()
    
    r=float(variable_dict["sigma"])/2.0
    box_size=float(variable_dict["boxSize"])
    dataCompression=int(variable_dict["dataCompression"])
    timestep=float(variable_dict["timestep"])
    numParticles=int(variable_dict["numParticles"])
    theta=variable_dict["theta"]
    
    parser = argparse.ArgumentParser(description='Plot energy data and save as PNG and HTML.')
    parser.add_argument('f1', type=str, help='Code Version')
    parser.add_argument('f2', type=str, help='Plot Version')
    args = parser.parse_args()
    plot_file_path = f"M6-V{args.f1}.{args.f2}"
    

    
    # Create a figure and axis for the animation
    fig, ax = plt.subplots()
    # Create the animation
    animation = FuncAnimation(fig, update, frames=len(particle_positions), interval=100)

    # Save the animation as a GIF
    animation.save('particle_animation.gif', writer='pillow')  # For GIF
    # Save the animation as a video
    # animation.save('particle_animation.mp4', writer='ffmpeg')  # For MP4
    
    energy_plotly_plotter(plot_file_path)
    persistence()
    asyncio.get_event_loop().run_until_complete(main(plot_file_path + "-Energy.png", "params.csv", "particle_animation.gif", "persistence.png", chat_id))
