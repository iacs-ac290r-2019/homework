'''
frame.py: at frame x, plot the temperature field, velocity (ux, uy) field,
pressure field and save figures in /temp_frames /ux_frames /uy_frames /p_frames

this script is loading the variables frame by frame, hence in theory would not
crash if dealing with reasonable size of resolution

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
p = ray.variables['vals_nod_var1']
temp = ray.variables['vals_nod_var2']
ux = ray.variables['vals_nod_var3']
uy = ray.variables['vals_nod_var4']
x = ray.variables['coordx']
y = ray.variables['coordy']
t = ray.variables['time_whole']

# get number of iterations and cells
niters, ncells = p.shape
# reshape x,y,t for plotting purpose
x_2d = np.reshape(x[:], [ny, nx])
y_2d = np.reshape(y[:], [ny, nx])
t_1d = np.reshape(t[:], [1, niters])

# function that save the temperature profile at a frame
def save_temp_frame(i):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,temp_i,100,cmap=plt.cm.get_cmap('RdBu'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=1.0, alpha=0.6)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Temperature at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame' + str(frame_number) + '.png')
    plt.close('all')
    # plt.show()

# function that save the ux profile at a frame
def save_ux_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,ux_i,100,cmap=plt.cm.get_cmap('RdBu'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=1.0, alpha=0.6)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_x at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame' + str(frame_number) + '.png')
    plt.close('all')
    # plt.show()

# function that save the uy profile at a frame
def save_uy_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,uy_i,100,cmap=plt.cm.get_cmap('RdBu'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=1.0, alpha=0.6)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_y at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame' + str(frame_number) + '.png')
    plt.close('all')
    # plt.show()

def save_p_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,p_i,100,cmap=plt.cm.get_cmap('RdBu'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=1.0, alpha=0.6)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Pressure at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame' + str(frame_number) + '.png')
    plt.close('all')
    # plt.show()

# save all temperature fields
for i in range(niters):
    # load variables step by step
    p_i = p[i].reshape(ny,nx)
    temp_i = temp[i].reshape(ny,nx)
    ux_i = ux[i].reshape(ny,nx)
    uy_i = uy[i].reshape(ny,nx)
#     save_temp_frame(i)
#     save_ux_frame(i)
#     save_uy_frame(i)
#     save_p_frame(i)
    print("Figures saved at Frame " + str(i))