# Author: Paige Pistono
# Email: ppistono@berkeley.edu

'''
This program uses the MDAnalysis Python packages and calculates MS2 viral capsid 
interdimer contact angles over the course of a trajectory. The center of mass of the
αA helices (residues 98-112) are calculated and a vector is drawn to represent the axis
between the two sets of helices in each dimer. Then, the angle between the two vectors
over the trajectory is calculated and plotted.
'''

import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')


# Load your trajectory and topology files (replace 'your_system.pdb' and 'your_trajectory.dcd' with your actual filenames)
u1 = mda.Universe("your_system.pdb", 
                 "your_trajectory.dcd")

# Select specific residues for the first set of helices on the first chain for the first dimer
residues_dimer1_set1 = u1.select_atoms('chainID A and resid 98:112')  # Replace 'A' with your actual chain ID and adjust residue range

# Select specific residues for the second set of helices on the first chain for the first dimer
residues_dimer1_set2 = u1.select_atoms('chainID B and resid 98:112')  # Replace 'A' with your actual chain ID and adjust residue range

# Select specific residues for the first set of helices on the second chain for the second dimer
residues_dimer2_set1 = u1.select_atoms('chainID C and resid 98:112')  # Replace 'C' with your actual chain ID and adjust residue range

# Select specific residues for the second set of helices on the second chain for the second dimer
residues_dimer2_set2 = u1.select_atoms('chainID D and resid 98:112')  # Replace 'C' with your actual chain ID and adjust residue range

# Initialize empty lists to store interdimer angles and frame numbers
interdimer_angles1 = []
frame_numbers = []

# Loop through frames and calculate interdimer angles
for ts in u1.trajectory:
	# Checking to make sure input files contain a structure
    if len(residues_dimer1_set1) > 0 and len(residues_dimer1_set2) > 0 and len(residues_dimer2_set1) > 0 and len(residues_dimer2_set2) > 0:
        # Calculate center of mass for the first set of helices in the first dimer
        center_of_mass_dimer1_set1 = residues_dimer1_set1.positions.mean(axis=0)

        # Calculate center of mass for the second set of helices in the first dimer
        center_of_mass_dimer1_set2 = residues_dimer1_set2.positions.mean(axis=0)

        # Calculate vector representing the first axis between the two sets of helices in the first dimer
        axis_vector_dimer1 = center_of_mass_dimer1_set2 - center_of_mass_dimer1_set1

        # Normalize the first axis vector
        axis_vector_dimer1_normalized = axis_vector_dimer1 / np.linalg.norm(axis_vector_dimer1)

        # Calculate center of mass for the first set of helices in the second dimer
        center_of_mass_dimer2_set1 = residues_dimer2_set1.positions.mean(axis=0)

        # Calculate center of mass for the second set of helices in the second dimer
        center_of_mass_dimer2_set2 = residues_dimer2_set2.positions.mean(axis=0)

        # Calculate vector representing the second axis between the two sets of helices in the second dimer
        axis_vector_dimer2 = center_of_mass_dimer2_set2 - center_of_mass_dimer2_set1

        # Normalize the second axis vector
        axis_vector_dimer2_normalized = axis_vector_dimer2 / np.linalg.norm(axis_vector_dimer2)

        # Calculate the angle between the two normalized vectors using the dot product
        interdimer_angle = np.arccos(np.dot(axis_vector_dimer1_normalized, axis_vector_dimer2_normalized))

        # Append angle and frame number to the lists
        interdimer_angles1.append(np.degrees(interdimer_angle))
        frame_numbers.append(ts.frame)

# Save interdimer angles to CSV
df_u1 = pd.DataFrame({'Frame': frame_numbers, 'Interdimer Angle (degrees)': interdimer_angles1})
df_u1.to_csv('filename.csv', index=False)

# Plotting
plt.plot(frame_numbers, interdimer_angles1, label='Interdimer Angle')
plt.xlabel('Frame Number')
plt.ylabel('Interdimer Angle (degrees)')
plt.title('Interdimer Angle Between Two Sets of Helices Over Time')
plt.legend()
plt.ylim(0,180)
plt.show()