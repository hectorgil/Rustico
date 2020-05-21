RUSTICO-X Rapid foUrier STatIstics COde for X-terms

Author: Hector Gil Marin

Last Release Date: 21st May 2020

email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu

====Compilation====

gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -fopenmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads_gcc/include/ -L/home/hector/fftw3_threads_gcc/lib/ -o rustico_gcc


icc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -openmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads/include/ -L/home/hector/fftw3_threads/lib/ -o rustico_icc

====Run====

./rustico_icc param_file.txt

====Parameter file options====

#Type of Box (periodic/cutsky): 'periodic' for periodic boxes with boundary conditions. 'cutsky' for actual observations or mocks which include sky mask

#Type of file (ascii/gadget): 'ascii' is the option required for 'cutsky'. 'periodic' option allows 'ascii' files or 'gadget' files. Gadget units assumed kpc/h. See 'ascii file structure' below for the format of the file

#Number of gadget files(int): In case the gadget boxes are split in more than 1 gadget file

#RSD distorsion on gadget periodic box (yes/no): yes. For gadget boxes allow this option to distort particles along the z-axis for redshift space distortions. The values of redhisft and Omega matter, will be taken from the header of the gadget file.

#Size of the Box (double/double): Low and Upper limits, respectively, of the cubic box where the galaxies are placed.

#Type of Computation (DSE/DSY/FFT): Type of Computation for the Power Spectrum: Direct Sum Exact (DSE); Direct Sum Yamamoto (DSY); Fast Fourier Transform (FFT). For the bispectrum computation FFT is required.

#Binning for the Power Spectrum (linear/log10): Binning type for the power spectrum output. Linear or 10-base logarithmic.

#Size of the bin for the power spectrum (double). Size of the bin for the power spectrum. In case log10 is choosen as binning, the interval is provided in log10-scale.

#k-range for computation (double/double): Low and Upper limits, respectively, of the k-values choosen for printing the power spectrum.

#Do Bispectrum (yes/no): Whether the bispectrum should be computed by the code

Do Multigrid (yes/no): Option for the bispectrum computation. If enable, the bispectrum triangles will be split according to their k-values and associated to different grid-sizes for a more optimal computation (large scale modes do not requires small grid cell ressolution). However, each grid-size computation will requires to re-associate the particles to the grid cells, which a potential lose of optimality. We recomend enable such option when many triangle shapes are required and when the datasets do not consists of many particles. In practice, each specific case will requires testing for determing the best performance option. When the multigrid option is enable, we require the interlacing option to be also enabled (see below), with at least 2 interlacing steps.

#Triangle Shapes (ALL/EQU/ISO/SQU). Triangle shapes to be computed. All (ALL), equilateral (EQU), Isosceles (ISO), squeezed (SQU). We define the squeezed triangles as those |k2-k3|<=k1 and K1<=0.1 K2; where by definition K1<=K2<=K3. Note that this condition is applied to the center of bin k-values and not to the effective k-values.

#Size of the bin for the bispectrum (double): Size of the bispectrum bin

#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): This determine the way how the Bispectrum is normalized: FFT is through performing Eq. xxx of ... using Fourier Transforms; APR is using the approximate solution of Eq. xxx 8pi^2k1k2k3Dk^3 and EXA is using the full and exact analytic solution (only available for Equilateral triangles at the moment). SUM option is computing APR and EXA for each triangle shape in the k-bin, whereas EFF is computing APR and EXA only for the effective k-values in the k-bin.

#Write triangles in each bin(yes/no): This option enables writting each triangle shape inside each k-bin. This option is only available for squeezed triangles (SQU choice enabled).

#Path for triangles in each bin: Determines the path for writting the above triangles

#Path of data: path of data file

#Path of randoms: path for the random file. If periodic box enabled, write 'nothing'

#Path of output: path were the output files will be written

#Identifier of output: identification string for the output files

#Write header: option for writting the header for power spectrum and bispectrum output files.

#Number of Grid Cells power (int): Number of grid-cells-per-side power input. If the input is n, the number of grid-cells per side will be 2^n. The input is an integer in the range 4<n<15

#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): mass interpolation assignment scheme: nearest-grid-point NGC, cloud-in-cell CIC, triangular-shaped-cloud TSC, piecewise-cubic-spline, piecewise-quartic-spline P4S, piecewise-quintic-spline P5S.

#Type of Yamamoto (GridCenter/GridAverage): Option only for skycut option. (details to be referred in a paper).

#Number of interlacing steps (int): Number of interlacing steps.

#Do Grid Correction? (yes/no): Grid correction option. (Jing et al 2005)

#Redshift Range (double/double): Low and Upper limits, respectively, for the redshift cuts in the skycut option

#Omega matter value (double): Omega matter value used for converting redshifts to comoving distances in the skycut option.

#Area effective value in deg^2 (double): Value of the area used for the normalization and of the power spectrum and bispectrum in the skycut option

#Hexadecapole as (L4/L2L2): Hexadecapole type of computation for the skycut option. L4 refers to Eq. xxx of Bianchi et al. 2015 and L2L2 to Eq. xxx of Scoccimarro et al. 2015.

#Compute Normalization as (area/density): Compute the normalization of the power spectrum using either the area value of the number density column of the input. Only for skycut option

#Compute Normalization using (randoms/data):Compute the normalization of the power spectrum and bispectrum using either the data or random n(z) computed from the objects. Only for skycut option

#Compute Shot noise as (double): Shot noise factor parameter. See. Gil-Marin et al. 2014 Eq. ..... Only for skycut option

====Ascii File Structure=====

For 'periodic' option: 4-column entry

x-position (double), y-position (double), z-position (double), weight (double)

For example,

2.234479 8.066518 6.439348 1.0

1.406139 5.003024 3.664719 1.0

9.559097 7.672305 10.529753 1.0

16.126962 1.955036 1.807368 1.0

13.377709 2.506474 2.225270 1.0

11.555976 0.431673 3.149166 1.0

14.122240 16.198794 3.532218 1.0

14.214268 14.341477 4.141918 1.0

For 'skycut' option: data and random entry required

#data. 8-column entry:

Right Ascension (double), declination (double), redshift (double), weight_fkp (double), weight_colision (integer), weight_systematics (double), number_density (double), veto (integer)

For example

1.291761884477e+02 4.894649924515e+01 5.425299000000e-01 1.238355000000e-01 1 1.029931000000e+00 7.075228831797e-04 1

1.174169630420e+02 3.927675925314e+01 3.996815000000e-01 6.395705000000e-01 1 1.020382000000e+00 5.635492881551e-05 1

1.169127239443e+02 3.944331088030e+01 5.377024000000e-01 1.231396000000e-01 1 9.972784000000e-01 7.120864449779e-04 0

1.169501719710e+02 3.949076919700e+01 5.191724000000e-01 1.143980000000e-01 1 1.061827000000e+00 7.741411563139e-04 1

1.175284705313e+02 4.017649329509e+01 5.431913000000e-01 1.246438000000e-01 1 1.022436000000e+00 7.022861947405e-04 1

1.238161587640e+02 4.663678411199e+01 5.896081000000e-01 1.800036000000e-01 1 1.009149000000e+00 4.555444446667e-04 1

1.276012766986e+02 4.977593704841e+01 5.481966000000e-01 1.309460000000e-01 1 1.053429000000e+00 6.636735753669e-04 1

#randoms. 6-column entry:

Right Ascension (double), declination (double), redshift (double), weight_fkp (double), number_density (double), veto (int)

For example

1.350240900000e+02 4.260718600000e+01 6.651101708412e-01 3.819904327393e-01 8.089333060371e-03 1

1.849009090000e+02 5.155070000000e+01 4.627352356911e-01 1.447689384222e-01 2.953779556923e-02 1

2.178783680000e+02 1.166351500000e+01 4.882844388485e-01 1.163647100329e-01 3.796835353765e-02 1

1.413639550000e+02 5.575799200000e+01 5.019413828850e-01 1.094641387463e-01 4.067705969521e-02 0

2.276985360000e+02 1.560757100000e+01 5.646098852158e-01 1.411796659231e-01 3.041586507737e-02 1

====Output Structure====

The power spectrum output have the following format: k-eff, k-centerbin, Monopole-Pshotnoise, Quadrupole, Hexadecapole, number of modes, Pshotnoise

The bispectrum output have the following format: k1-eff, k1-centerbin, k2-eff, k2-centerbin, k3-eff, k3-centerbin, B0-Bshotnoise, Bshotnoise, Reduced Bispectrum, Reduced Bispectrum shot noise, number of triangles

For the skycut option, the code automatically generates two extra file for the number density of objects as a function of redshift for the data and random files with the following format,

#data

#Interval: xxx Mpc/h

#z < nobs > < wc nobs> < wc wfkp nobs>

#random

#Interval: yyy Mpc/h

#alpha< ns > alpha < wfkp ns >

where xxx and yyy are the interval choosen by the code to bin the data and randoms, respectively, < nobs > are the raw number density of observed objects, < wc nobs> is the number density of observed objects weighted by the collision weights, and < wc wfkp nobs> weighted by the collision weights and fkp weights. For the randoms, these number densities are scaled by alpha in order to match the data ones.

====Citation====

If you use this code for your published or unpublished work, please refer it to Gil-Marin, Hector in prep. 2017

====Disclaimer====

blah blah blah

