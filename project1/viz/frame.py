'''
frame.py: at frame x, plot the temperature field, velocity (ux, uy) field,
pressure field and save figures in /temp_frames /ux_frames /uy_frames /p_frames

this script is loading the variables frame by frame, hence in theory would not
crash if dealing with reasonable size of resolution

'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import time
import netCDF4

# set all contour lines to be solid
# u_x u_y fields some contour lines are dottted to indicate negativity
plt.rcParams['contour.negative_linestyle'] = 'solid'

# shape of the grid
ny = 1025 # subject to change 1025
nx = 2049 # subject to change 2049

# open the Exodus file
path = '../data/case1.exo' # subject to change on Odyssey
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
def save_temp_frame_mesh(frame_number, t_max = 1.0):
    fig, ax = plt.subplots(figsize=(15,6))
    fig.subplots_adjust(top=0.9,right=0.9)

    cs = ax.pcolormesh(x_2d, y_2d, temp[frame_number].reshape(ny,nx), cmap=plt.cm.get_cmap('jet'))
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Temperature at Frame %04d' % frame_number, fontsize=20)
    ax.set_aspect(1.0)
    cs.set_clim(vmin=0.0, vmax=t_max)
    fig.colorbar(cs,ax=ax)
    plt.savefig('temp_frames_mesh/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

def save_temp_frame_contour(frame_number, t_max = 1.0):
    fig, ax = plt.subplots(figsize=(15,6))
    fig.subplots_adjust(top=0.9,right=0.9)

    cs = ax.contourf(x_2d,y_2d,temp[frame_number].reshape(ny,nx),levels=np.linspace(0, 1, 11),cmap=plt.cm.get_cmap('jet'))
    # cs.set_clim(0.0,1.0)
    ax.contour(cs, colors='k',linewidths=2.0)
    cs.set_clim(vmin=0.0,vmax=t_max)
    fig.colorbar(cs,ax=ax)

    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Temperature at Frame %04d' % frame_number, fontsize=20)
    ax.set_aspect(1.0)
    plt.savefig('temp_frames_contour/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()
    
# function that save the ux profile at a frame
def save_ux_frame(frame_number):
    fig, ax = plt.subplots(figsize=(15,6))
    fig.subplots_adjust(top=0.9,right=0.9)
    cs = ax.contourf(x_2d,y_2d,ux[frame_number].reshape(ny,nx),10,cmap=plt.cm.get_cmap('jet'))
    # cs.set_clim(0.0,1.0)
    ax.contour(cs, colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_x at Frame %04d' % frame_number, fontsize=20)
    # ax.set_aspect(1.0)
    # cs.set_clim(vmin=0.0,vmax=1.0)
    fig.colorbar(cs,ax=ax)
    plt.savefig('ux_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

# function that save the uy profile at a frame
def save_uy_frame(frame_number):
    fig, ax = plt.subplots(figsize=(15,6))
    fig.subplots_adjust(top=0.9,right=0.9)
    cs = ax.contourf(x_2d,y_2d,uy[frame_number].reshape(ny,nx),10,cmap=plt.cm.get_cmap('jet'))
    # cs.set_clim(0.0,1.0)
    ax.contour(cs, colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('U_y at Frame %04d' % frame_number, fontsize=20)
    ax.set_aspect(1.0)
    # cs.set_clim(vmin=0.0,vmax=1.0)
    fig.colorbar(cs,ax=ax)
    plt.savefig('uy_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

def save_p_frame(frame_number):
    fig, ax = plt.subplots(figsize=(15,6))
    fig.subplots_adjust(top=0.9,right=0.9)
    cs = ax.contourf(x_2d,y_2d,p[frame_number].reshape(ny,nx),10,cmap=plt.cm.get_cmap('jet'))
    # cs.set_clim(0.0,1.0)
    ax.contour(cs, colors='k',linewidths=2.0)
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_title('Pressure at Frame %04d' % frame_number, fontsize=20)
    ax.set_aspect(1.0)
    # cs.set_clim(vmin=0.0,vmax=1.0)
    fig.colorbar(cs,ax=ax)
    plt.savefig('p_frames/frame%04d.png' % frame_number)
    plt.close('all')
    # plt.show()

# save all fields
# start timer
t0 = time.time()
# Number of frames
iMax = len(p)
# Status
print(f'Processing {iMax} frames...')
for i in range(iMax):
    # Set maximum temperature for temp frames
    t_max = 0.25
    # Save the frames
    save_temp_frame_mesh(i, t_max)
    save_temp_frame_contour(i, t_max)
    save_ux_frame(i)
    save_uy_frame(i)
    save_p_frame(i)
    # compute elapsed time
    et = time.time() - t0
    # Compute ETA
    pace = (i+1) / et
    eta = (iMax - i - 1) / pace
    print(f"Figures saved at Frame {i:04d}; elapsed {int(et)} sec; ETA {int(eta)} sec.")
