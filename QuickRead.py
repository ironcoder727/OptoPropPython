# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 19:55:28 2024

@author: ENABLEH2_2
"""

import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the text file (assuming it has headers)
file_path = 'propeller_EXAMPLE_ACP_Sweep.txt'  # Replace this with your actual file path

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