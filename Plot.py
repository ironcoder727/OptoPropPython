# -*- coding: utf-8 -*-
from pathlib import Path


import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the text file (assuming it has headers)
file_path = Path.cwd() / 'output' / 'Baseline_Propeller_Sweep.txt' # Replace this with your actual file path

# Read the CSV data, pandas will automatically use the first row as headers
df = pd.read_csv(file_path, sep=';')

# Extract the 'T' (Thrust) column
rpm =  df['rpm']

thrust = df['T']

# Plot the thrust data
plt.plot(rpm,thrust, marker='o', linestyle='-', color='b')
plt.title('Thrust vs rpm')
plt.xlabel('rpm')
plt.ylabel('Thrust (N)')
plt.grid(True)

# Show the plot
plt.show()