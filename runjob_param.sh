#!/bin/bash
#SBATCH -J rbk
#SBATCH -t 30:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
LC_ALL="C"
OMP_NUM_THREADS=6
#setenv OMP_NUM_THREADS 12
export OMP_NUM_THREADS=6
##parametrit: Qs^2 gamma x0 C^2 ec freeze output
#if [ ! -f ${7} ]
#then
	./rbk -ic MV $1 $2 $3 $5 -alphas_scaling $4 -minr 1e-6 -rc BALITSKY -alphas_freeze_c ${6} -output ${7} -maxy 16 -ystep 0.2
#fi

