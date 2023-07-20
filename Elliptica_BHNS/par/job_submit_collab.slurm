#! /bin/bash
### submit jobs with dependencies: X=1; jobid=`sbatch job_submit_expanse.slurm|tail -n1 | sed s'/Submitted batch job //'`; echo "Launching jobs:"; echo $jobid; for i in `seq ${X}`; do jobid=`sbatch --dependency=afterany:$jobid job_submit_expanse.slurm|tail -n1 | sed s'/Submitted batch job //'`;echo $jobid;done; SQ;

#SBATCH -t 48:00:00                               # time for simulation
###SBATCH --account=open
#SBATCH --account=dur566_c_gpu
###SBATCH --partition=sla-prio 
#SBATCH --partition=burst
#SBATCH --qos=burst2x 
#SBATCH --nodes=5                                # number of nodes to use
#SBATCH --ntasks-per-node=48                    # number of cores per node
#SBATCH --export=ALL                             # keep environment
#SBATCH --mail-type=ALL                          # get emails for anything that happens
#SBATCH --mail-user=pespino@berkeley.edu    # your email address
#SBATCH -J G2_LOR                             # Job name
#SBATCH -o G2_LOR.out                         # outfile
#SBATCH -e G2_LOR.err                         # errfile

# load modules
module load intel/2021.4.0
module load impi/2021.4.0
module load python/3.6.8
module load hdf5/1.12.1
module load fftw/3.3.10
module load boost/1.77.0
module load mkl/2021.4.0

# move to your working directory
cd /storage/home/ple5069/scratch/G2_BNS

# Export MPI/OMP environment variables
##for openmpi3
export OMPI_MCA_btl_openib_if_include="mlx5_2:1";export OMPI_MCA_btl=self,openib,vader;export CACTUS_NUM_THREADS=6; export OMP_NUM_THREADS=6; export CACTUS_NUM_PROCS=40
#export OMPI_MCA_btl_openib_if_include="mlx5_2:1";export OMPI_MCA_btl=self,openib,vader;export CACTUS_NUM_THREADS=1; export OMP_NUM_THREADS=1; export CACTUS_NUM_PROCS=128
#for openmpi4
#export OMPI_MCA_btl=self,vader; export CACTUS_NUM_THREADS=1; export OMP_NUM_THREADS=1; export CACTUS_NUM_PROCS=896

# Note start time and start job!
export CACTUS_STARTTIME=$(date +%s)
time srun --mpi=pmi2 --cpus-per-task=$OMP_NUM_THREADS --ntasks=$CACTUS_NUM_PROCS ./cactus_etk G2_LOR.par -reo | tee CCTK_Proc0.out
