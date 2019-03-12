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

# set all contour lines to be solid
# u_x u_y fields some contour lines are dottted to indicate negativity
plt.rcParams['contour.negative_linestyle'] = 'solid'

ny = 1025 # subject to change 1025
nx = 2049 # subject to change 2049

# open the Exodus file
path = 'case1.exo' # subject to change on Odyssey
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
def save_temp_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    ax.pcolormesh(x_2d,y_2d,temp_i,cmap=plt.cm.get_cmap('jet'))
    # cs = ax.contourf(x_2d,y_2d,temp_i,100,cmap=plt.cm.get_cmap('jet'))
    # cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Temperature at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

# function that save the ux profile at a frame
def save_ux_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,ux_i,100,cmap=plt.cm.get_cmap('jet'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_x at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('ux_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

# function that save the uy profile at a frame
def save_uy_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,uy_i,100,cmap=plt.cm.get_cmap('jet'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_y at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('uy_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

def save_p_frame(frame_number):
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(111)
    cs = ax.contourf(x_2d,y_2d,p_i,100,cmap=plt.cm.get_cmap('jet'))
    cs2 = ax.contour(cs, levels=cs.levels[::10], colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Pressure at Frame %i' % i, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('p_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

# save all temperature fields
for i in range(0,50):
    # load variables step by step
    # p_i = p[i].reshape(ny,nx)
    # temp_i = temp[i].reshape(ny,nx)
    ux_i = ux[i].reshape(ny,nx)
    # uy_i = uy[i].reshape(ny,nx)
    # save_temp_frame(i)
    save_ux_frame(i)
#     save_uy_frame(i)
#     save_p_frame(i)
# color the temp by velocity or 
# color the streamline by velocity_magnitude/temp
    print("Figures saved at Frame %04d" % i)