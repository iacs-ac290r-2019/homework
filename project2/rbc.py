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

def plot_style() -> None:
    """Set plot style for the session."""
    # Set up math plot library to use TeX
    # https://matplotlib.org/users/usetex.html
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)
    # Set default font size to 20
    mpl.rcParams.update({'font.size': 20})

def drug_delivery(frame_num):
    """Compute the amount of drug delivered to stenotic region at the frame number"""
    # Attributes at this frame
    rho, vel, drug = load_frame(frame_num)
    # Amount of drug in stenotic region
    drug_sten = np.sum(drug[mask_sten])
    # Total amount of drug in the system
    drug_tot = np.sum(drug)
    return (drug_sten, drug_tot)

def calc_drug(frame_nums):
    """Compute the drug quantity (stenotic and total) at each frame"""
    # Number of frame
    frame_count = len(frame_nums)
    # Preallocate storage for drug in stenotic region and whole system
    drug_sten = np.zeros(frame_count)
    drug_tot = np.zeros(frame_count)

    # Compute drug at each available frame and save it to arrays
    for i, frame_num in enumerate(frame_nums):
        # Drug delivery for this frame
        sten, tot = drug_delivery(frame_num)
        # Save to two arrays
        drug_tot[i] = tot
        drug_sten[i] = sten
    
    # Compute the drug fraction in the stenotic region
    drug_frac = drug_sten / drug_tot
    
    return drug_sten, drug_tot, drug_frac

def plot_drug(drug_tot, drug_frac, frame_time):
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

# Position data is shared by all frames
point_pos, cell_pos, cell_vol = load_pos()

# Extract x, y and z from cell_pos; we work primarily with cells, not points
# These are the centers of each cell
x = cell_pos[:, 0]
y = cell_pos[:, 1]
z = cell_pos[:, 2]

# Total volume of the simulation
vol_tot = np.sum(cell_vol)

# Mask for stenotic region
mask_sten = (sten_beg <= z) & (z < sten_end)

# Volume in the stenotic region
vol_sten = np.sum(cell_vol[mask_sten])

# List of available frames
frame_nums = frames_avail(dir_np)
# Convert frame numbers to time in milliseconds
frame_time = 0.001 * np.array([frame_num for frame_num in frame_nums])
# Number of frames
frame_count = len(frame_nums)

# Set plot style
plot_style()

# Calculate drug quantity vs. time
drug_sten, drug_tot, drug_frac = calc_drug(frame_nums)

# Plot the total amount of drug and delivery fraction vs. timef
plot_drug(drug_tot, drug_frac, frame_time)

# Status update
d_min = np.min(drug_tot)
d_max = np.max(drug_tot)
d_chng_pct = (d_max - d_min) / d_min * 100
print(f'Change in total drug quantity over simulation: {d_chng_pct:0.2f}%')
