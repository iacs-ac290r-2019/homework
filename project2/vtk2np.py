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
    Load position data for a frame number and specified directory (blood or drug).
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
    point_pos = np.array(ds.points)
    # Extract the center of each cell
    cell_pos = ds.cell_centers().points
    # Compute the volume of each cell
    cell_vol = ds.compute_cell_sizes().cell_arrays['Volume'].astype(np.float32)

    # Return the points grid
    return point_pos, cell_pos, cell_vol

def load_frame_type(frame_num, dir_abs):
    """
    Load the VTKI grid object for a frame number and specified directory (blood or drug).
    INPUTS:
    frame_num: the frame number to load, e.g. 100000
    dir_abs: the absolute directory of data to load, e.g. 
             'd:/iacs/ac290/project2/data/RBC_0_Re10/DIRDATE_BloodFlow/VTK'
    """
    # The filename of the .pvtu file for the requested frame
    filename = os.path.join(dir_abs, f'T{frame_num:010}.pvtu')
    # Read this file with VTKI
    ds = vtki.read(filename)
    return ds
    
def load_frame_vtu(frame_num):
    """Load field data for this frame (blood and drug) from vtu files"""
    # Get the blood data from VTKI
    ds_blood = load_frame_type(frame_num, dir_blood)

    # Get the velocity vel at POINTS
    vel = ds_blood.point_arrays['velocity']
    
    # Get the density rho at CELLS
    rho = ds_blood.cell_arrays['density']
    # vel = ds_blood.cell_arrays['velocity']

    # Extract drug data for the frame
    ds_drug = load_frame_type(frame_num, dir_drug)
    # Extract the density of the drug by cell; this is in the confusingly named field temperature
    drug = ds_drug.cell_arrays['temperature']

    # Extract the density of the drug at points as well for contour plots
    drug_point = ds_drug.point_arrays['density']
    
    # Return a tuple of 4 elements: velocity and density for blood and drug
    return (rho, vel, drug, drug_point)

def save_frame(frame_num, dir_np):
    """Load frame from VTU, save it as numpy array"""
    # Convert from VTU to numpy arrays
    rho, vel, drug, drug_point = load_frame_vtu(frame_num)
    
    # Table of arrays and corresponding file name prefixes
    file_tbl = {
        'rho': rho,
        'vel': vel,
        'drug': drug,
        'drug_point': drug_point
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

    # Get the position for points and cells from the first frame
    frame_num = frame_nums[0]
    point_pos, cell_pos, cell_vol = load_frame_pos(frame_num, dir_blood)

    # Save the point positions
    fname_point_pos = os.path.join(dir_np, 'point_pos.npy')
    np.save(fname_point_pos, point_pos)
    # Save the cell positions
    fname_cell_pos = os.path.join(dir_np, 'cell_pos.npy')
    np.save(fname_cell_pos, cell_pos)
    # Save the cell volumes
    fname_cell_vol = os.path.join(dir_np, 'cell_vol.npy')
    np.save(fname_cell_vol, cell_vol)
    
    # Set interval; can't take all frames, too slow
    interval = 10000
    
    # Iterate over frame numbers, saving numpy arrays for each one
    for frame_num in tqdm.tqdm(frame_nums):
        if frame_num % interval == 0:
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
    
if __name__ == '__main__':
    # Save all the frames
    save_frames(dir_np)

