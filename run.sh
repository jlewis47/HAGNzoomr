#!/bin/sh
#PBS -S /bin/sh
#PBS -N CosmoZoomTest
#PBS -j oe
#PBS -l nodes=2:ppn=128,walltime=48:00:00

module purge
module load inteloneapi/2022.2

date

cd /
mpirun -genv FI_PROVIDER=tcp -np 256 ./ramses3d ./cosmo.nml >& run_00000.log

date

exit 0
