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
import tqdm

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

def load_drug(frame_num):
    """Load the drug density keyed by cell for this frame number"""
    basename = f'drug_{frame_num:07}.npy'
    filename = os.path.join(dir_np, basename)
    return np.load(filename)

def load_drug_point(frame_num):
    """Load the drug density keyed by point for this frame number"""
    basename = f'drug_point_{frame_num:07}.npy'
    filename = os.path.join(dir_np, basename)
    return np.load(filename)

def load_vel(frame_num):
    """Load the velocity keyed by point for this frame number"""
    basename = f'vel_{frame_num:07}.npy'
    filename = os.path.join(dir_np, basename)
    return np.load(filename)

def load_rho(frame_num):
    """Load the density rho keyed by cell for this frame number"""
    basename = f'rho_{frame_num:07}.npy'
    filename = os.path.join(dir_np, basename)
    return np.load(filename)

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
    # Drug concentration keyed by cell
    drug = load_drug(frame_num)
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

def plot_drug_ts(drug_tot, drug_frac, vol_frac, frame_time):
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
    # Time in milliseconds
    t_ms = frame_num // 1000
    
    # Check if file exists; if so, quit early
    filename = os.path.join(dir_fig, f'contour_speed/contour_speed_z_{z_plot}_t_{t_ms}.png')
    if os.path.isfile(filename):
        return

    # Mask for the selected z
    mask = (point_z == z_plot)
    
    # x and y on this cross section
    xx = point_x[mask]
    yy = point_y[mask]
    
    # Attributes at this frame
    vel = load_vel(frame_num)
    # Reduce to mask
    vel = vel[mask]
    # Velocity in the x, y and z direction
    u_x = vel[:, 0]
    u_y = vel[:, 1]
    u_z = vel[:, 2]
    # Velocity magnitude
    uu = np.sqrt(u_x**2 + u_y**2 + u_z**2)
    
    # Interpolation points
    pi = np.arange(3, 103, 1, dtype=np.int32)
    num_pi = pi.shape[0]
    xi = pi
    yi = pi
    
    # Index for augmented arrays
    idx_xx = xx.astype(np.int32) - xi[0]
    idx_yy = yy.astype(np.int32) - yi[0]

    # Create a padded array on the full grid
    uu_pad = -1
    uu_aug = uu_pad * np.ones((num_pi, num_pi))
    # Fill in uu_aug with real data; use numpy indexing for speed
    uu_aug[idx_xx, idx_yy] = uu
    
    # levels for contour
    num_levels = 8
    uu_min = np.min(uu)
    uu_max = np.max(uu)
    levels = np.linspace(uu_min, uu_max, num_levels)
    
    # color map for contour
    cmap = plt.cm.get_cmap('jet')
    
    # Shift x and y down by 2.5 to put them from the range (3, 102) to (0.5, 99.5)
    # This is just so plots will line up with sensible scales
    # Back this off by 0.5 to account for the cell vs. point difference and center the plot.
    shift = 2.0
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
    ax.set_xticks(np.arange(xMin, xMax+10, 10))
    ax.set_yticks(np.arange(yMin, yMax+10, 10))

    # Make the contour plot
    cs = ax.contourf(xi, yi, uu_aug, levels=levels, cmap=cmap)
    ax.contour(cs)
    ax.grid()
    fig.colorbar(cs, ax=ax)
    
    # Save file
    fig.savefig(filename, bbox_inches='tight')
    plt.close()

# *********************************************************************************************************************
def plot_streamlines(z_plot, frame_num):
    """Create a 2D stream plot longitude z=z_plot and the given frame number"""
    # Time in milliseconds
    t_ms = frame_num // 1000
    
    # Check if file exists; if so, quit early
    filename = os.path.join(dir_fig, f'streamlines/streamlines_z_{z_plot}_t_{t_ms}.png')
    if os.path.isfile(filename):
        return

    # Mask for the selected z
    mask = (point_z == z_plot)
    
    # x and y on this cross section
    xx = point_x[mask]
    yy = point_y[mask]
    
    # Attributes at this frame
    vel = load_vel(frame_num)
    # Reduce to mask
    vel = vel[mask]
    # Velocity in the x and y directions
    u = vel[:, 0]
    v = vel[:, 1]
    # Speed in the xy plane; for line width
    speed = np.sqrt(u*u + v*v)
    # Line width as a vector; proportional to speed
    lw_vec = speed / np.mean(speed) * 2.0
    
    # Interpolation points
    pi = np.arange(3, 103, 1, dtype=np.int32)
    num_pi = pi.shape[0]
    xi = pi
    yi = pi
    
    # Index for augmented arrays
    idx_xx = xx.astype(np.int32) - xi[0]
    idx_yy = yy.astype(np.int32) - yi[0]
    
    # Create a padded array on the full grid for the x velocity, u
    u_aug = np.zeros((num_pi, num_pi))
    # Fill in uu_aug with real data; use numpy indexing for speed
    u_aug[idx_xx, idx_yy] = u
    
    # Create a padded array on the full grid for the x velocity, u
    v_aug = np.zeros((num_pi, num_pi))
    # Fill in uu_aug with real data; use numpy indexing for speed
    v_aug[idx_xx, idx_yy] = v
    
    # 2D array for linewidth
    lw = np.zeros((num_pi, num_pi))
    lw[idx_xx, idx_yy] = lw_vec
    
    # levels for contour
    density = 2.0
    
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
    fig, ax = plt.subplots(figsize=[12,12])
    ax.set_title(f'Streamlines at $z=${z_plot}, $T$={t_ms} ms')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    ax.set_xticks(np.arange(xMin, xMax+10, 10))
    ax.set_yticks(np.arange(yMin, yMax+10, 10))

    # Make the contour plot
    ax.streamplot(x=xi, y=yi, u=u_aug, v=v_aug, density=density, color='blue', linewidth=lw, arrowsize=1.5)
    ax.grid()
    # fig.colorbar(cs, ax=ax)
    
    # Save file
    fig.savefig(filename, bbox_inches='tight')
    plt.close()


# *********************************************************************************************************************
def plot_drug_conc(z_plot, frame_num):
    """Create a 2D contour plot of drug concentration"""
    # Time in milliseconds
    t_ms = frame_num // 1000
    
    # Check if file exists; if so, quit early
    filename = os.path.join(dir_fig, f'contour_drug/contour_drug_z_{z_plot}_t_{t_ms}.png')
    if os.path.isfile(filename):
        return

    # Mask for the selected z
    mask = (cell_z == z_plot + 0.5)
    
    # x and y on this cross section
    xx = cell_x[mask] - 0.5
    yy = cell_y[mask] - 0.5
    
    # Drug concentration at this frame
    drug = load_drug(frame_num)
    # Reduce to mask
    drug = drug[mask]
    
    # Interpolation points
    pi = np.arange(3, 103, 1, dtype=np.int32)
    num_pi = pi.shape[0]
    xi = pi
    yi = pi
    
    # Index for augmented arrays
    idx_xx = xx.astype(np.int32) - xi[0]
    idx_yy = yy.astype(np.int32) - yi[0]

    ## Create a padded array on the full grid
    drug_pad = -1
    drug_aug = drug_pad * np.ones((num_pi, num_pi))
    # Fill in drug_aug with real data; use numpy indexing for speed
    drug_aug[idx_xx, idx_yy] = drug

    # levels for contour
    num_levels = 8
    drug_min = np.min(drug)
    drug_max = np.max(drug)
    levels = np.linspace(drug_min, drug_max, num_levels)
    if drug_min == drug_max:
        eps = 1e-6
        levels = np.linspace(drug_min*(1-eps), drug_min*(1+eps), num_levels)
    
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
    ax.set_title(f'Contours of Drug Density at $z=${z_plot}, $T$={t_ms} ms')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    ax.set_xticks(np.arange(xMin, xMax+10, 10))
    ax.set_yticks(np.arange(yMin, yMax+10, 10))
    
    # Make the contour plot
    cs = ax.contourf(xi, yi, drug_aug, levels=levels, cmap=cmap)
    ax.contour(cs)
    ax.grid()
    fig.colorbar(cs, ax=ax)
    
    # Save file
    fig.savefig(filename, bbox_inches='tight')
    plt.close()

# *********************************************************************************************************************
def plot_drug_profile(frame_num):
    """Compute the mean drug profile vs. longitude z"""
    # Time in milliseconds
    t_ms = frame_num // 1000
    
    # Check if file exists; if so, quit early
    filename = os.path.join(dir_fig, f'drug_profile/drug_profile_t_{t_ms}.png')
    if os.path.isfile(filename):
        return

    # Drug concentration at this frame
    drug = load_drug(frame_num)
    
    zz = cell_sd_z
    drug_mean = np.zeros_like(zz)
    for i, z in enumerate(zz):
        # Mask for cells on this slice
        mask = (cell_z == zz[i])
        # Mean concentration on this slice
        drug_mean[i] = np.mean(drug[mask])
    
    # Smooth out "holes" by setting drug concentration to be minimum of drug mean and its two neighbors
    N = zz.shape[0]
    dm = drug_mean
    dm[1:N-1] = np.maximum(dm[0:N-2], dm[1:N-1], dm[2:N-0])
    
    # Plot smoothed mean drug concentration vs. z at this frame   
    fig, ax = plt.subplots(figsize=[16,9])
    ax.set_title(f'Cross Sectional Drug Concentration at T={t_ms} ms')
    ax.set_xlabel('Midpoint of $z$')
    ax.set_ylabel('Mean Drug Concentration')
    ax.plot(zz, dm, color='b', linewidth=3.0, label='Mean Conc')
    # Manually set lower y limit to 0
    ylim = ax.get_ylim()
    ax.set_ylim(0.0, ylim[1])
    ax.grid()

    # Save figure
    fig.savefig(filename, bbox_inches='tight')
    plt.close()

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
dir_np = os.path.join(os.getcwd(), '../data/RBC_0_Re10/Numpy')
# Directory to outout figures
dir_fig = os.path.join(os.getcwd(), '../figs')

# List of available frames
frame_nums = frames_avail(dir_np)

# Position data is shared by all frames
point_pos, cell_pos, cell_vol = load_pos()

# Extract x, y and z from cell_pos
# These are the centers of each cell
cell_x = cell_pos[:, 0]
cell_y = cell_pos[:, 1]
cell_z = cell_pos[:, 2]

# Sorted, distinct z
cell_sd_z = np.array(sorted(set(cell_z)))

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

# *********************************************************************************************************************
def main():
    # Set plot style
    plot_style()
    
    # Calculate drug quantity vs. time
    drug_sten, drug_tot, drug_frac = calc_drug(frame_nums, mask_sten)
    
    # Plot the total amount of drug and delivery fraction vs. time
    plot_drug_ts(drug_tot, drug_frac, vol_frac, frame_time)
    
    # Status update
    d_min = np.min(drug_tot)
    d_max = np.max(drug_tot)
    d_chng_pct = (d_max - d_min) / d_min * 100
    print(f'Change in total drug quantity over simulation: {d_chng_pct:0.2f}%')
    
    # Sample frame numbers
    frame_nums_plot = np.array([100, 200, 300, 400, 600, 800, 1000])*1000
    
    # Sample points for longitudinal plots
    z_plot = [200, 300, 400, 500, 600, 700, 800]
    
    # Plot the contour of speed at the middle of the stenosis
    print('Speed contour.')
    for frame_num in tqdm.tqdm(frame_nums_plot):
        for z in z_plot:
            plot_speed_contour(z, frame_num)

    # Plot the streamlines of the velocity field in the XY plane at the middle of the stenosis at selected times
    print('\nStreamlines.')
    for frame_num in tqdm.tqdm(frame_nums_plot):
        for z in z_plot:
            plot_streamlines(z, frame_num)    

    # Plot the drug concententration cross-sectionally for selected z and time
    print(f'\nDrug concentration cross section.')
    for frame_num in tqdm.tqdm(frame_nums_plot):
        for z in z_plot:
            plot_drug_conc(z, frame_num)

    # Plot the drug profile vs. z for all frames
    print(f'Drug profile vs. z.')
    for frame_num in tqdm.tqdm(frame_nums):
        plot_drug_profile(frame_num)

# *********************************************************************************************************************
# main()
