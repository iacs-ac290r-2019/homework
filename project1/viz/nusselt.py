"""
Michael S. Emanuel
Tue Mar 12 20:43:18 2019
"""

import numpy as np
import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set plot style
mpl.rcParams.update({'font.size': 20})

# open the Exodus file
path = '../data/case1.exo' # subject to change on Odyssey
ray = netCDF4.Dataset(path)

# grid size
ny = 1025
nx = 2049

# parameters of the simulation
H = 2.0 # height of the channel
k = 1.0 # thermal diffusivity kappa
dT = 1.0 # temperature difference at steady state 

# mapping between NetCDF variable names and actual meaning
# p = ray.variables['vals_nod_var1']
temp = ray.variables['vals_nod_var2']
# ux = ray.variables['vals_nod_var3']
uy = ray.variables['vals_nod_var4']
# x = ray.variables['coordx']
# y = ray.variables['coordy']
time = ray.variables['time_whole']

def calc_nusselt():
    # Number of time steps
    tMax = len(temp)
    # Initialize Nu to an array of zeros
    Nu = np.zeros(tMax)
    # Initialize tt to an array of zeros; the time
    tt = np.zeros(tMax)
    # Populate the Nusselt number at each time step
    for t in range(tMax):
        # The velocity in the y direction at this time step
        uy_t = uy[t].reshape(ny,nx)
        # The temperature field at this time step
        temp_t = temp[t].reshape(ny, nx)
        # The product of the y-velocity and temperature
        uT = uy_t * temp_t
        # the change in temperature dT from the bottom to the top
        # set dT to a minimum of 1E-6 to avoid divide by zero at simulation start
        # dT = max(temp_t[0, 0] - temp_t[ny-1, 0], 1E-6)
        # the nusselt number is Nu = 1 + H / k dT * <vT>
        Nu[t] = 1.0 + H / (k * dT) + np.mean(uT)
        # Save the time
        tt[t] = time[t] 
    # Return the nusselt number and time
    return Nu, tt

# Calculate the nusselt number
Nu, tt = calc_nusselt()

# Plot a chart of the Nusselt number
fig, ax = plt.subplots(figsize=[12,12])
ax.set_title('Nusselt Number for RBC Case 1')
ax.set_xlabel('Time t')
ax.set_ylabel('Nusselt Number')
ax.set_xlim(0.0, 0.25)
ax.set_ylim(2.0, 4.0)
ax.plot(tt, Nu, color='r')
ax.grid()
fig.savefig('frames_misc/nusselt.png', bbox_inches='tight')
