#!/bin/sh
#PBS -S /bin/sh
#PBS -N CosmoZoomTest
#PBS -j oe
#PBS -l nodes=2:ppn=128,walltime=48:00:00

module purge
module load inteloneapi/2022.2

date


#avoid overwriting the last log file... this gets the highest number and adds 1
cur_nb=$(printf "%d" $(find . -maxdepth 1 -name 'run_*.log' -type f | cut -d _ -f 2 | cut -d . -f 1 | sort | tail -n 1))
nb=$(printf "%05d" $((cur_nb+1)))

cd /
mpirun -genv FI_PROVIDER=tcp -np 256 ./ramses3d ./cosmo.nml >& run_$nb.log

date

exit 0
