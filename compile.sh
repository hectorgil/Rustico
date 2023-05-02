#compilation with GCC
gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -fopenmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads_gcc/include/ -L/home/hector/fftw3_threads_gcc/lib/ -o file.out

#Compilation with ICC (recomended if avaliable, the code may run faster). 
icc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -openmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads/include/ -L/home/hector/fftw3_threads/lib/ -o file.out

export OMP_NUM_THREADS=16 # number of threads you will require, otherwise it take the maximum available.
