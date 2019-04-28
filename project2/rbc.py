"""
Harvard AC 290R
Project 2: Blood Flow Simulation

Michael S. Emanuel
Sat Apr 27 15:58:59 2019
"""

import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# *********************************************************************************************************************
def plot_style() -> None:
    """Set plot style for the session."""
    # Set up math plot library to use TeX
    # https://matplotlib.org/users/usetex.html
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)
    # Set default font size to 20
    mpl.rcParams.update({'font.size': 20})

# *********************************************************************************************************************
def load_pos():
    """Load positions of grid mesh points"""
    # Load the point positions
    filename_point_pos = os.path.join(dir_np, 'point_pos.npy')
    point_pos = np.load(filename_point_pos)

    # Load the cell positions
    filename_cell_pos = os.path.join(dir_np, 'cell_pos.npy')
    cell_pos = np.load(filename_cell_pos)

    # Load the cell volumes
    filename_cell_vol = os.path.join(dir_np, 'cell_vol.npy')
    cell_vol= np.load(filename_cell_vol)

    # Return the three position arrays
    return point_pos, cell_pos, cell_vol

def load_frame(frame_num):
    """Load frame from cached numpy arrays (faster)"""
    # Table of arrays and corresponding file name prefixes
    prefixes = ['rho','vel', 'drug']
    # Table of arrays
    tbl = dict()
    for prefix in prefixes:
        basename = f'{prefix}_{frame_num:07}.npy'
        filename = os.path.join(dir_np, basename)
        tbl[prefix] = np.load(filename)
    # Return a tuple (rho_b, vel_b, rho_d, vel_d)
    return tuple(tbl.values())

def frames_avail(dir_np):
    """Get list of available frame numbers"""
    # All available files for blood density
    filenames = glob.glob(os.path.join(dir_np, 'rho_*.npy'))
    # Extract the frame numbers
    basenames = [os.path.basename(filename) for filename in filenames]
    frame_nums = [int(basename.split('.')[0][-7:]) for basename in basenames]
    return frame_nums

# *********************************************************************************************************************
def drug_delivery(frame_num, mask_sten):
    """Compute the amount of drug delivered to stenotic region at the frame number"""
    # Attributes at this frame
    rho, vel, drug = load_frame(frame_num)
    # Amount of drug in stenotic region
    drug_sten = np.sum(drug[mask_sten])
    # Total amount of drug in the system
    drug_tot = np.sum(drug)
    return (drug_sten, drug_tot)

def calc_drug(frame_nums, mask_sten):
    """Compute the drug quantity (stenotic and total) at each frame"""
    # Number of frame
    frame_count = len(frame_nums)
    # Preallocate storage for drug in stenotic region and whole system
    drug_sten = np.zeros(frame_count)
    drug_tot = np.zeros(frame_count)

    # Compute drug at each available frame and save it to arrays
    for i, frame_num in enumerate(frame_nums):
        # Drug delivery for this frame
        sten, tot = drug_delivery(frame_num, mask_sten)
        # Save to two arrays
        drug_tot[i] = tot
        drug_sten[i] = sten
    
    # Compute the drug fraction in the stenotic region
    drug_frac = drug_sten / drug_tot
    
    return drug_sten, drug_tot, drug_frac

def plot_drug(drug_tot, drug_frac, vol_frac, frame_time):
    """Create two plots with drug quantity over time"""
    # Compute volume fraction in the stenotic region
    vol_frac = vol_sten / vol_tot
    
    # Plot 1: Total Drug in System (Quality Check)
    fig, ax = plt.subplots(figsize=[16,9])
    ax.set_title('Total Drug vs. Time')
    ax.set_xlabel('Time in Milliseconds')
    ax.set_ylabel('Total Drug in System')
    ax.set_ylim(0, 120000)
    ax.plot(frame_time, drug_tot, color='b', linewidth=3.0, marker='o', markersize=10, label='Total Drug')
    ax.grid()
    # Save figure
    fname_tot = os.path.join(dir_fig, 'drug_tot.png')
    fig.savefig(fname_tot, bbox_inches='tight')
        
    # Plot 2: Fraction of Drug Delivered
    fig, ax = plt.subplots(figsize=[16,9])
    ax.set_title('Drug Delivery Fraction vs. Time')
    ax.set_xlabel('Time in Milliseconds')
    ax.set_ylabel('Drug Delivery \%')
    ax.set_ylim(0.0, 3.0)
    ax.plot(frame_time, drug_frac*100, color='b', linewidth=3.0, marker='o', markersize=10, label='Drug \%')
    # second series is a baseline; what fraction of the volume is in the stenotic region?
    vol_frac_vec = vol_frac*100*np.ones_like(frame_time)
    ax.plot(frame_time, vol_frac_vec, color='r', linewidth=3.0,label='Volume \%')
    ax.grid()
    ax.legend()
    # Save figure
    fname_frac = os.path.join(dir_fig, 'drug_frac.png')
    fig.savefig(fname_frac, bbox_inches='tight')

# *********************************************************************************************************************
def plot_speed_contour(z_plot, frame_num):
    """Create a 2D contour plot of the speed at longitude z=z_max and the given frame number"""
    # Mask for the selected z
    mask = (point_z == z_plot)
    
    # x and y on this cross section
    xx = point_x[mask]
    yy = point_y[mask]
    
    # Time in milliseconds
    t_ms = frame_num // 1000
    
    # Attributes at this frame
    _, vel, _ = load_frame(frame_num)
    # Reduce to mask
    vel = vel[mask]
    # Velocity in the x, y and z direction
    u_x = vel[:, 0]
    u_y = vel[:, 1]
    u_z = vel[:, 2]
    # Velocity magnitude
    uu = np.sqrt(u_x**2 + u_y**2 + u_z**2)
    
    # Interpolation points
    pi = np.linspace(3, 102, 100)
    num_pi = pi.shape[0]
    xi = pi
    yi = pi
    
    # Create a padded array on the full grid
    uu_pad = -1
    uu_aug = uu_pad * np.ones((num_pi, num_pi))
    # Fill in uu_aug with real data; use numpy indexing for speed
    uu_aug[xx.astype(np.int32), yy.astype(np.int32)] = uu
    
    # levels for contour
    num_levels = 8
    uu_min = np.min(uu)
    uu_max = np.max(uu)
    levels = np.linspace(uu_min, uu_max, num_levels)
    
    # color map for contour
    cmap = plt.cm.get_cmap('jet')
    
    # Shift x and y down by 2.5 to put them from the range (3, 102) to (0.5, 99.5)
    # This is just so plots will line up with sensible scales
    shift = 2.5
    xi = xi-shift
    yi = yi-shift
    
    # x and y limits for plots
    xMin = np.floor(np.min(xx-shift)/10)*10
    xMax = np.ceil(np.max(xx-shift)/10)*10
    yMin = np.floor(np.min(yy-shift)/10)*10
    yMax = np.ceil(np.max(yy-shift)/10)*10
    
    # Make axes
    fig, ax = plt.subplots(figsize=[15,12])
    ax.set_title(f'Contours of Speed $|U|$ at $z=${z_plot}, $T$={t_ms} ms')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)

    # Make the contour plot
    cs = ax.contourf(xi-2.5, yi-2.5, uu_aug, levels=levels, cmap=cmap)
    ax.contour(cs)
    ax.grid()
    fig.colorbar(cs, ax=ax)
    
    # Save file
    filename = os.path.join(dir_fig, f'contour_speed_z_{z_plot}_t_{t_ms}.png')
    fig.savefig(filename, bbox_inches='tight')

# *********************************************************************************************************************
# Geometry constants
# Radius of the artery in healthy region
R = 50
# Length of the artery; on z axis
L = 20 * R
# Half length; center point of stenosis
A = L//2
# Length of stenosis
S = L//10
# Cross section in the stenosis
# G = R/2.0
# Start of the bolus
# d = L/4
# Width of the bolus
# W = L/20

# Start and end of stenosis
sten_beg = A - S//2
sten_end = A + S//2

# The directory with numpy arrays
dir_np = os.path.join(os.getcwd(), 'data/RBC_0_Re10/Numpy')
# Directory to outout figures
dir_fig = os.path.join(os.getcwd(), 'figs')

# List of available frames
frame_nums = frames_avail(dir_np)

# Position data is shared by all frames
point_pos, cell_pos, cell_vol = load_pos()

# Extract x, y and z from cell_pos
# These are the centers of each cell
cell_x = cell_pos[:, 0]
cell_y = cell_pos[:, 1]
cell_z = cell_pos[:, 2]

# Extract x, y, and z from point_pos
point_x = point_pos[:, 0]
point_y = point_pos[:, 1]
point_z = point_pos[:, 2]

# Mask for stenotic region
mask_sten = (sten_beg <= cell_z) & (cell_z < sten_end)
# Volume in the stenotic region
vol_sten = np.sum(cell_vol[mask_sten])
# Total volume of the simulation
vol_tot = np.sum(cell_vol)
# Fraction of volume in the stenotic region
vol_frac = vol_sten / vol_tot

# Convert frame numbers to time in milliseconds
frame_time = 0.001 * np.array([frame_num for frame_num in frame_nums])
# Number of frames
frame_count = len(frame_nums)

def main():
    # Set plot style
    plot_style()
    
    # Calculate drug quantity vs. time
    drug_sten, drug_tot, drug_frac = calc_drug(frame_nums, mask_sten)
    
    # Plot the total amount of drug and delivery fraction vs. time
    plot_drug(drug_tot, drug_frac, vol_frac, frame_time)
    
    # Status update
    d_min = np.min(drug_tot)
    d_max = np.max(drug_tot)
    d_chng_pct = (d_max - d_min) / d_min * 100
    print(f'Change in total drug quantity over simulation: {d_chng_pct:0.2f}%')

# *********************************************************************************************************************

frame_num = 100000
z_plot = A


plot_speed_contour(A, 1000000)
