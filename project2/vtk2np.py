"""
Harvard AC 290R
Project 2: Blood Flow Simulation

Michael S. Emanuel
Sat Apr 27 09:47:17 2019
"""

import os
import glob
import numpy as np
import vtki
import tqdm

# *********************************************************************************************************************
def load_frame_pos(frame_num, dir_abs):
    """
    Load all available data for a frame number and specified directory (blood or drug).
    INPUTS:
    frame_num: the frame number to load, e.g. 100000
    dir_abs: the absolute directory of data to load, e.g. 
             'd:/iacs/ac290/project2/data/RBC_0_Re10/DIRDATE_BloodFlow/VTK'
    """
    # The filename of the .pvtu file for the requested frame
    filename = os.path.join(dir_abs, f'T{frame_num:010}.pvtu')
    # Read this file with VTKI
    ds = vtki.read(filename)
    # Extract the points
    points = np.array(ds.points)

    # Return the points grid
    return points

def load_frame_type(frame_num, dir_abs):
    """
    Load all available data for a frame number and specified directory (blood or drug).
    INPUTS:
    frame_num: the frame number to load, e.g. 100000
    dir_abs: the absolute directory of data to load, e.g. 
             'd:/iacs/ac290/project2/data/RBC_0_Re10/DIRDATE_BloodFlow/VTK'
    """
    # The filename of the .pvtu file for the requested frame
    filename = os.path.join(dir_abs, f'T{frame_num:010}.pvtu')
    # Read this file with VTKI
    ds = vtki.read(filename)
    # Extract the velocity and density
    density = ds.point_arrays['density']
    velocity = ds.point_arrays['velocity']

    # Return the points grid, density and velocity (mapped to points)
    return (density, velocity)

def load_frame_vtu(frame_num):
    """Load all data for this frame (blood and drug) from vtu files"""
    # Extract blood data for this frame; note that the position is the same for blood and drug
    rho_b, vel_b = load_frame_type(frame_num, dir_blood)
    # Extract drug data for the frame
    rho_d, vel_d = load_frame_type(frame_num, dir_drug)
    # Return a tuple of 4 elements: velocity and density for blood and drug
    return (vel_b, rho_d, vel_d, rho_d)

def save_frame(frame_num, dir_np):
    """Load frame from VTU, save it as numpy array"""
    # Convert from VTU to numpy arrays
    rho_b, vel_b, rho_d, vel_d = load_frame_vtu(frame_num)
    
    # Table of arrays and corresponding file name prefixes
    file_tbl = {
        'rho_b': rho_b,
        'vel_b': vel_b,
        'rho_d': rho_d,
        'vel_d': vel_d
        }
    
    # Save each grid as a numpy array
    for fname, arr in file_tbl.items():
        fname_prefix = os.path.join(dir_np, fname)
        fname = f'{fname_prefix}_{frame_num:07}.npy'
        np.save(fname, arr)

def save_frames(dir_np):
    """Save frames for all available VTU files as numpy arrays"""
    # Build frame table from blood file names; assume drug file names are parallel (otherwise program will fail)
    filenames = glob.glob(os.path.join(dir_blood, '*.pvtu'))
    num_frames = len(filenames)
    print(f'Identified {num_frames} frames.')
    
    # Extract the frame numbers
    basenames = [os.path.basename(filename) for filename in filenames]
    frame_nums = [int(basename.split('.')[0][-7:]) for basename in basenames]

    # Get the points from the first frame
    frame_num = frame_nums[0]
    pos = load_frame_pos(frame_num, dir_blood)
    # Save the points grid
    fname_pos = os.path.join(dir_np, 'pos.npy')
    np.save(fname_pos, pos)
    
    # Iterate over frame numbers, saving numpy arrays for each one
    for frame_num in tqdm.tqdm(frame_nums[0:1]):
        save_frame(frame_num, dir_np)
        print(f'Completed frame {frame_num:7}.')

# *********************************************************************************************************************
# The root directory
dir_root = os.path.join(os.getcwd(), 'data/RBC_0_Re10')
# The relative directory with blood data
dir_blood_rel = 'DIRDATA_BloodFlow/VTK'
# The relative directory with drug data
dir_drug_rel = 'DIRDATA_Bolus/VTK'
# Directory to save numpy arrays
dir_np_rel = 'Numpy'

# Absolute directories
dir_blood = os.path.join(dir_root, dir_blood_rel)   
dir_drug = os.path.join(dir_root, dir_drug_rel)   
dir_np = os.path.join(dir_root, dir_np_rel)   

# Save all the frames
save_frames(dir_np)
