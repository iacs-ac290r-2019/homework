#!/bin/sh

#SBATCH --partition=fas_gpu
#SBATCH --time=02:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --job-name=jacobi
#SBATCH --output=jacobi-320.out
#SBATCH --res=ac290r_gpu

# Can run commands here, or use ‘echo’ to look at CUDA_VISIBLE_DEVICES
hostname
echo $CUDA_VISIBLE_DEVICES
./parallel 320
