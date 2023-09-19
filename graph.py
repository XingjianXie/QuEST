import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import numpy.ma as ma
import re

num_qubits = 0

# Specify the layout
num_rows = 11
num_cols = 2

# Define colors for different prefixes
prefix_colors = {
    'Prompt': 'Red',
    'QuEST/Prompt': 'Green',
    'Prompt-512': 'Pink',
    'QuEST': 'Blue',
    'HpQC': 'Purple'
}

special = ["BV", "GS"]

# Read lines from stdin
lines = sys.stdin.readlines()

# Remove newline characters and filter out empty lines
lines = [line.strip() for line in lines if line.strip()]

# Create empty dictionaries to store the data
other_data = {}
controlled_data = {}

prefix = ""

i = 0

found = False

# Process each pair of lines
while i < len(lines):
    if "====" in lines[i]:
        i += 1
        found = False
        continue
    if not found and "----" not in lines[i]:
        i += 1
        continue
    if "----" in lines[i]:
        found = True
        num_qubits = int(re.findall(r'\d+', lines[i+1])[0])
        i += 2
        continue
    if "----" in lines[i + 1]:
        prefix = re.findall(r"\[([^\]]*)\]", lines[i])[0]
        i += 2
        continue
    # Get the key from the first line and remove the trailing colon
    key = lines[i].replace(":", "")
    if "controlled" in key or "multiControlled" in key:
        # values should be a matrix of num_qubits * num_qubits
        values = [[float(v) if 'N/A' not in v else np.nan for v in lines[i+1+j].split(',') if v] for j in range(num_qubits)]
        values = np.array(values)
        values = ma.masked_array(values, np.isnan(values))
        i += num_qubits + 1
        # Store the values in the controlled_data dictionary
        if key not in controlled_data:
            controlled_data[key] = {}
        controlled_data[key][prefix] = values
    else:
        # Convert the string of values on the next line into a list of floats
        values = [float(v) for v in lines[i+1].split(',') if v]
        values = np.array(values)
        i += 2
        # Store the values in the other_data dictionary
        if key not in other_data:
            other_data[key] = {}
        other_data[key][prefix] = values

# Create a figure
fig = plt.figure(figsize=(50, num_rows*32))

# Adjust the vertical spacing between subplots
plt.subplots_adjust(hspace=0.5)

# For each key, create a 2D line plot
for i, (key, key_data) in enumerate(other_data.items()):
    ax = fig.add_subplot(num_rows, num_cols, i + 1)

    # Process each prefix-value pair in the data dictionary
    for prefix, values in key_data.items():
        # Plot the values with the prefix as the label
        ax.plot(values, label=prefix, color=prefix_colors[prefix])

    ax.set_xticks(range(num_qubits))

    if key in special:
        ax.set_xticklabels(range(1, num_qubits + 1))

    # Add a legend
    ax.legend()

    # Set the title of the subplot to the key
    ax.set_title(key)

    if not key in special:
        ax.set_xlabel('targetQubit')
        ax.set_ylabel('time')
    else:
        ax.set_xlabel('number of qubits')
        ax.set_ylabel('time')
        ax.set_yscale("log")
    
    key_data["QuEST/Prompt"] = key_data["QuEST"] / key_data["Prompt"]
    ax2 = ax.twinx() 
    ax2.set_ylabel('QuEST time / Prompt time', color = prefix_colors["QuEST/Prompt"])
    ax2.plot(key_data["QuEST/Prompt"], color = prefix_colors["QuEST/Prompt"])
    ax2.tick_params(axis ='y', labelcolor = prefix_colors["QuEST/Prompt"]) 

# For each key in controlled_data, create a 2D subplot for the heatmaps
subplot_index = len(other_data)  # Create a separate subplot index

if subplot_index % num_cols != 0:
    subplot_index += num_cols - (subplot_index % num_cols)

for i, (key, key_data) in enumerate(controlled_data.items()):
    min = np.min([np.min(v) for v in key_data.values()])
    min = min - (min % 50)
    max = np.max([np.max(v) for v in key_data.values()])
    max = max + (50 - (max % 50))
    # Process each prefix-value pair in the data dictionary
    for j, (prefix, values) in enumerate(key_data.items()):
        ax = fig.add_subplot(num_rows, num_cols, subplot_index + j + 1)  # Use the separate subplot index

        # Add the heatmap to the subplot. The parameter 'cmap' is set to the color associated with the prefix.
        im = ax.imshow(values, cmap="plasma", alpha=0.5, interpolation='nearest', vmin=min, vmax=max)

        # Add a colorbar to the plot
        plt.colorbar(im, ax=ax)

        # Set xticks and yticks to every integer position
        ax.set_xticks(range(num_qubits))
        ax.set_yticks(range(num_qubits))

        # Adding numbers
        for (k, l), z in np.ndenumerate(values):
            if not values.mask[k, l]:  # Only add text if the value is not masked
                ax.text(l, k, '{:0.0f}'.format(z), ha='center', va='center', color='black')

        # Set the title of the subplot to the key
        ax.set_title(f"{key} - {prefix}")

        if "controlled" in key:
            # Add x and y axis labels
            ax.set_xlabel('targetQubit')
            ax.set_ylabel('controlQubit')
        
        if "multiControlled" in key:
            # Add x and y axis labels
            ax.set_xlabel('maskEndingQubit')
            ax.set_ylabel('maskStartingQubit')

        # Move the x-axis to the top of the plot
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 

    subplot_index += len(key_data)  # Increase the subplot index by the number of plots for each key

# Save the plot to a file
plt.savefig("output.png")
