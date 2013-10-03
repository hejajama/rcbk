#!/bin/bash
#SBATCH -J rbk
#SBATCH -t 10:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
LC_ALL="C"
OMP_NUM_THREADS=4
export OMP_NUM_THREADS=4
module add intel
##parametrit: Qs^2 gamma x0 C^2 ec freeze output
echo Running on `hostname`
echo Modules:
echo `module list`
#if [ ! -f ${7} ]
#then
	./rbk -ic MV $1 $2 $3 $5 -alphas_scaling $4 -minr 1e-6 -rc BALITSKY -alphas_freeze_c ${6} -output ${7} -maxy 11 -ystep 0.2
#fi

