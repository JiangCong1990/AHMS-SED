#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=2:00:00
#SBATCH --account=qxia1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=qxia1@uni-koeln.de
#SBATCH --output=runlog.txt
#SBATCH --error=errorlog.txt

####SBATCH account=AG-Shao
# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)

# module load compiler
#module load intel/15.0 intelmpi/4.1.0 hdf5/1.8.13 netcdf/4.1.3 szlib/2.1
module load gnu openmpi/1.6.5 hdf5/1.8.13 netcdf/4.1.3 szlib/2.1

export WRF_HYDRO=1
export HYDRO_D=0

ulimit -s unlimited

date

# MPI Running
#srun -n $SLURM_NTASKS ./wrf.exe
#./real.exe
#srun ./wrf.exe

# Serial Running
./wrf_hydro_NoahMP.exe

date

