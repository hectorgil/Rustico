#!/bin/bash
#SBATCH -p normal #partition (queue) ##Don't change this
#SBATCH -N 1 #number of nodes
#SBATCH -n 2 #number of cores
#SBATCH --mem 45000 #memory pool ##Mb of RAM
#SBATCH -t 20-20:00 # time limit (D-HH:MM)
#SBATCH -o /home/hector/rustico_v52/test.o #path file to store output
#SBATCH -e /home/hector/rustico_v52/test.e #path file to store error messages
#SBATCH --job-name test #job name

#cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=2

gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_linedata.c cubature.c ps_write.c order_algorithm.c -O3   -lm -fopenmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads/include/ -L/home/hector/fftw3_threads/lib/ -o file.out

time ./file.out params_boss.c
