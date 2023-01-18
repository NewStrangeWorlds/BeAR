#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=daniel.kitzmann@csh.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=all

# Job name
#SBATCH --job-name="helios-r2 WASP-43b"

# Runtime and memory
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=15G

# Workdir
#SBATCH --chdir=/storage/homefs/kitzmann/code/helios-r2_retrievals/wasp43-b

# Partition
#SBATCH --partition=all

# For parallel jobs
#SBATCH --cpus-per-task=3
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4

# For GPU
#SBATCH --partition=gpu-invest
#SBATCH --qos=job_gpu_Heng
#SBATCH --gres=gpu:gtx1080ti:1

#SBATCH --signal=CONT@120


#SBATCH --output=WASP-43b/cluster_output/output.txt
#SBATCH --error=WASP-43b/cluster_output/error.txt

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10

#### Your shell commands below this line ####
#module load Boost
#module load GCC
#module load CUDA
srun ./helios-r WASP-43b/ -r
