# Author: Paige Pistono
# Email: ppistono@berkeley.edu

'''
This program uses the MDAnalysis Python packages and calculates MS2 viral capsid interdimer contact angles over 
the course of a trajectory. The centers of mass of the neighboring residues 36 and 98 were calculated, and a 
vector was drawn to represent the axis between the residues in each neighboring dimer.Then, the angle between 
the two vectors over the trajectory is calculated and plotted.
'''

import numpy as np
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt


# Load your trajectory and topology files (replace 'your_system.pdb' and 'your_trajectory.dcd' with your actual filenames)
u1 = mda.Universe("your_system.pdb", 
                 "your_trajectory.dcd")

# Select specific residues and chains for the interdimer interface angle calculation (residues 36 and 98)
residue_dimer1 = u1.select_atoms('chainID M and resid 36') # Replace 'M' with your actual chain ID
residue_dimer2 = u1.select_atoms('chainID L and resid 98') # Replace 'L' with your actual chain ID

# Initialize empty lists to store interdimer angles and frame numbers
interdimer_angles = []
frame_numbers = []

# Loop through frames and calculate interdimer angles
for ts in u1.trajectory:
    if len(residue_dimer1) > 0 and len(residue_dimer2) > 0:
        # Calculate vector between the two selected residues
        vector_dimer1 = residue_dimer1.positions.mean(axis=0)
        vector_dimer2 = residue_dimer2.positions.mean(axis=0)
        vector = vector_dimer2 - vector_dimer1

        # Calculate the interdimer angle using dot product
        interdimer_angle = np.arccos(np.dot(vector_dimer1, vector_dimer2) / (np.linalg.norm(vector_dimer1) * np.linalg.norm(vector_dimer2)))

        # Append angle and frame number to the lists
        interdimer_angles.append(np.degrees(interdimer_angle))
        frame_numbers.append(ts.frame)

# Calculate average, standard deviation, and standard error of the mean
average_angle = np.mean(interdimer_angles)
std_dev = np.std(interdimer_angles)
sem = std_dev / np.sqrt(len(interdimer_angles))

# Save interdimer angles to CSV
df = pd.DataFrame({'Frame': frame_numbers, 'Interdimer Angle (degrees)': interdimer_angles})
df.to_csv('009_transpro_ML_interdimer_angles_36_98.csv', index=False)

# Plotting
plt.plot(frame_numbers, interdimer_angles, label='Interdimer Angle')
plt.xlabel('Frame Number')
plt.ylabel('Interdimer Angle (degrees)')
plt.title('Interdimer Angle Between Residue 36 and Residue 98 Over Time')
plt.legend()
plt.ylim(0, 180)
plt.show()

# Print calculated statistics
print("Average Angle:", average_angle)
print("Standard Deviation:", std_dev)
print("Standard Error of the Mean:", sem)




