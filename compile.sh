gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c mask.c structures.c -O3   -lm -fopenmp -lfftw3_omp  -lfftw3 -lm  -o file.out
