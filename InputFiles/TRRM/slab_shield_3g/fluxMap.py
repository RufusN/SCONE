import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the JSON data from a file
file_name = '/Users/ru/SCONE/InputFiles/TRRM/slab_shield_3g/slab_shield_3g_planes.json'  # Replace with your actual file path

try:
    with open(file_name, 'r') as f:
        data = json.load(f)

    # Access the 'responseRate' field
    response_rates = data['fluxMap']['fluxMap']

    # response_rates = response_rates[:][:][1]

    # Extract the first elements from the nested lists and organize into a 36x20 grid
    heatmap_data = np.array([[pair[0] for pair in outer_list] for outer_list in response_rates])

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, cmap="viridis", cbar=True)

    # Modify the plot: flip y-axis and remove ticks
    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.yticks([])

    plt.title('')
    plt.xlabel('')
    plt.ylabel('')
    plt.show()

except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found. Please check the file path.")
except KeyError as e:
    print(f"Error: Missing expected key in the JSON data - {e}")
except Exception as e:
    print(f"An error occurred: {e}")
