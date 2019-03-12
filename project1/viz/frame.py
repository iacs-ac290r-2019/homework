'''
frame.py: at frame x, plot the temperature field, velocity (ux, uy) field,
pressure field and save figures in /temp_frames /ux_frames /uy_frames /p_frames

'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import netCDF4

ny = 65 # subject to change 1025
nx = 65 # subject to change 2049

# open the Exodus file
path = 'test.exo' # subject to change on Odyssey
ray = netCDF4.Dataset(path)

# mapping between NetCDF variable names and actual meaning
# transform NetCDF variables into np arrays
p = np.array(ray.variables['vals_nod_var1'])
temp = np.array(ray.variables['vals_nod_var2'])
ux = np.array(ray.variables['vals_nod_var3'])
uy = np.array(ray.variables['vals_nod_var4'])
x = np.array(ray.variables['coordx']).reshape(ny,nx)
y = np.array(ray.variables['coordy']).reshape(ny,nx)
t = np.array(ray.variables['time_whole'])

# function that save the temperature profile at a frame
def save_temp_frame(frame_number):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(x,y,temp[frame_number].reshape(ny,nx),cmap=plt.cm.get_cmap('RdBu'))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title('Temperature Field at Frame %i' % frame_number)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame' + str(frame_number) + '.png')
    plt.close('all')

# function that save the ux profile at a frame
def save_ux_frame(frame_number):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(x,y,ux[frame_number].reshape(ny,nx),cmap=plt.cm.get_cmap('RdBu'))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title('Ux Field at Frame %i' % frame_number)
    ax.set_aspect(1.0)
    plt.savefig('ux_frames/frame' + str(frame_number) + '.png')
    plt.close('all')

# function that save the uy profile at a frame
def save_uy_frame(frame_number):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(x,y,uy[frame_number].reshape(ny,nx),cmap=plt.cm.get_cmap('RdBu'))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title('Uy Field Frame %i' % frame_number)
    ax.set_aspect(1.0)
    plt.savefig('uy_frames/frame' + str(frame_number) + '.png')
    plt.close('all')

# function that save the pressure profile at a frame
def save_p_frame(frame_number):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(x,y,p[frame_number].reshape(ny,nx),cmap=plt.cm.get_cmap('RdBu'))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title('Pressure Field at Frame %i' % frame_number)
    ax.set_aspect(1.0)
    plt.savefig('p_frames/frame' + str(frame_number) + '.png')
    plt.close('all')

# save all temperature fields
niters, ncells = temp.shape
for i in range(niters):
    save_temp_frame(i)
    save_ux_frame(i)
    save_uy_frame(i)
    save_p_frame(i)
    print("Figures saved at Frame " + str(i))