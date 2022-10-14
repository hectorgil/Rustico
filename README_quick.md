RUSTICO-X Rapid foUrier STatIstics COde for X-terms

Author: Hector Gil Marin

Last Release Date: 14th Oct 2022

email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu

====Compilation====

gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -fopenmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads_gcc/include/ -L/home/hector/fftw3_threads_gcc/lib/ -o rustico_gcc


icc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c structures.c mask.c -O3   -lm -openmp -lfftw3_omp  -lfftw3 -lm -I/home/hector/fftw3_threads/include/ -L/home/hector/fftw3_threads/lib/ -o rustico_icc

====Run====

./rustico_icc param_file.txt

==== Considerations ===

When running the code for the very first time with the bispectrum option enabled for, the code will create some "wisdom" files in order to speed up the computation of the bispectrum for further runs. These files will depend on the grid used: wisdom64, wisdom128, wisdom256 etc. Please, allow for the 1st time a single process to run, so these files can be created without an overlapping writting errors. Once these files are created, several processs can simultaneously run, as they will be required only to be read. 

====Parameter file options====

#Type of code (rustico/rusticoX): rustico option for auto-statistics, rusticoX option for cross statistics.

#Type of Box (periodic/cutsky): 'periodic' for periodic boxes with boundary conditions. 'cutsky' for actual observations or mocks which include sky mask

#Type of file (ascii/gadget)x2: 'ascii' is the option required for 'cutsky'. 'periodic' option allows 'ascii' files or 'gadget' files. Gadget units assumed kpc/h. See 'ascii file structure' below for the format of the file. In case of rusticoX option is selected 2 inputs are required, for the 2 object-field to cross-correlate. 

#Type of input (density/particles): 'particles' is the recommended option, if your input consist of particles/galaxies/haloes etc. In case you start from a density grid, then you should select .density'.

#Number of gadget files(int)x2: In case the gadget boxes are split in more than 1 gadget file. For rusticox two inputs are required, each for each input path above. 

#RSD distorsion on gadget periodic box (yes/no)x2: yes. For gadget boxes allow this option to distort particles along the z-axis for redshift space distortions. The values of redhisft and Omega matter, will be taken from the header of the gadget file. For rusticox two inputs are required, each for each input path above. 

#Size of the Box (double/double): Low and Upper limits, respectively, of the cubic box where the galaxies are placed.

#Type of Computation (DSE/DSY/FFT): Type of Computation for the Power Spectrum: Direct Sum Exact (DSE); Direct Sum Yamamoto (DSY); Fast Fourier Transform (FFT). For the bispectrum computation FFT is required. Recommended FFT.

#Binning for the Power Spectrum (linear/log10): Binning type for the power spectrum output. Linear or 10-base logarithmic.

#Size of the bin for the power spectrum (double). Size of the bin for the power spectrum. In case log10 is choosen as binning, the interval is provided in log10-scale.

#k-range for computation (double/double): Low and Upper limits, respectively, of the k-values choosen for printing the power spectrum.

#Do anisotropy signal (yes/no): Whether only the monopole or higher multipoles are computed (quadrupole and hexadecapole). By default set to "yes". 

#Do odd multipoles (yes/no): Whether Dipole and octoploe are computed. Recommended "no", unless you need those for a particular project. 

#Write kvectors in each bin(yes/no): This option will write all the k-vector (kx,ky,kz) for each k-bin. By default set to 'no'.

#Do mu-binning Power Spectrum (yes/no): no Bin P(k,mu), only available for periodic boxes. By default set to 'no'.  

#Number of mu-bins (int): 120 Number of mu bins between 0 and 1. Only relevant if Do mu-binning is set to yes. 

#Different files for mu-bin (yes/no): no Whether the mu bins are all writen in the same or different output files. 

#Do Bispectrum (yes/no): Whether the bispectrum should be computed by the code

#Do Bispectrum multipoles (yes/no): Whether the bispectrum quadrupole(s) 200, 020, 002 are computed (not available for rusticoX). 

#Do Multigrid (yes/no): Option for the bispectrum computation. If enable, the bispectrum triangles will be split according to their k-values and associated to different grid-sizes for a more optimal computation (large scale modes do not requires small grid cell ressolution). However, each grid-size computation will requires to re-associate the particles to the grid cells, which a potential lose of optimality. We recomend enable such option when many triangle shapes are required and when the datasets do not consists of many particles. In practice, each specific case will requires testing for determing the best performance option. When the multigrid option is enable, we require the interlacing option to be also enabled (see below), with at least 2 interlacing steps.

#Bispectrum optimization (lowmem/standard/himem): Different optimization techniques. By default set it to 'standard'. Himem option should run faster but will require higher RAM. If lot of RAM avaliable set it to himem. 

#Triangle Shapes (ALL/EQU/ISO/SQU). Triangle shapes to be computed. All (ALL), equilateral (EQU), Isosceles (ISO), squeezed (SQU). We define the squeezed triangles as those |k2-k3|<=k1 and K1<=0.1 K2; where by definition K1<=K2<=K3. Note that this condition is applied to the center of bin k-values and not to the effective k-values.

#Size of the bin for the bispectrum (double): Size of the bispectrum bin

#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): This determine the way how the Bispectrum is normalized: FFT is through performing Eq. xxx of ... using Fourier Transforms; APR is using the approximate solution of Eq. xxx 8pi^2k1k2k3Dk^3 and EXA is using the full and exact analytic solution (only available for Equilateral triangles at the moment). SUM option is computing APR and EXA for each triangle shape in the k-bin, whereas EFF is computing APR and EXA only for the effective k-values in the k-bin. By default set it to "FFT". 

#Write triangles in each bin(yes/no): This option enables writting each triangle shape inside each k-bin. This option is only available for squeezed triangles (SQU choice enabled). By default set it to "no".

#Path for triangles in each bin: Determines the path for writting the above triangles

#Path of datax2: path of data file. Two inputs required in case of rusticoX

#Path of randomsx2: path for the random file. If periodic box enabled, write 'none'. Two inputs required in case of rusticoX

#Path of output: path were the output files will be written

#Identifier of output: identification string for the output files

#Write header: option for writting the header for power spectrum and bispectrum output files.

#Write density: Option for writting out the density field in the grid. By default set it to "no". 

#Number of Grid Cells power (int): Number of grid-cells-per-side power input. If the input is n, the number of grid-cells per side will be 2^n. The input is an integer in the range 4<n<15. For example 8 stands for 2^8=256 grid per side. By default set it to "8"; then increase or reduce depending on your needs. Exceeding 10 or 11 may take a lot of time to run. 

#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): mass interpolation assignment scheme: nearest-grid-point NGC, cloud-in-cell CIC, triangular-shaped-cloud TSC, piecewise-cubic-spline, piecewise-quartic-spline P4S, piecewise-quintic-spline P5S. Recommended "PCS", but can increase it if it doesn't penalize with time (it really depends on the compiler). 

#Type of Yamamoto (GridCenter/GridAverage): Option only for skycut option. (details to be referred in a paper). 'GridCenter' strongly recommended (unless you understand what these two mean). 

#Number of interlacing steps (int): Number of interlacing steps. By default 2 it's enough. 

#Do Grid Correction? (yes/no): Grid correction option. (Jing et al 2005). By default set it to "yes". 

#Redshift Range (double/double)x2: Low and Upper limits, respectively, for the redshift cuts in the skycut option. Two pairs of inputs required in case of rusticoX

#Omega matter value (double): Omega matter value used for converting redshifts to comoving distances in the skycut option.

#Area effective value in deg^2 (double)x2: Value of the area used for the normalization and of the power spectrum and bispectrum in the skycut option. Two inputs required in case of rusticoX

#Smoothing scale for normalization in Mpc/h (double): "it sets the scale at which n(z) is smoothed to compute the normalization. If set to 0 an automahic value is choosen (but it requires more time). By default set it to 1. 

#Quadrupole as (L0L2/L1L1): Options for projecting the line-of-sight for the quadrupole. (see pdf for more details). Recommended L0L2

#Octopole as (L0L3/L1L2): Options for projecting the line-of-sight for the octopole. (see pdf for more details). Recommended L0L3

#Hexadecapole as (L0L4/L2L2/L1L3): Options for projecting the line-of-sight for the hexadecapole. (see pdf for more details). Recommended L2L2

#Compute Normalization as (area/density): Compute the normalization of the power spectrum using either the area value of the number density column of the input. Only for skycut option

#Compute Normalization using (randoms/data):Compute the normalization of the power spectrum and bispectrum using either the data or random n(z) computed from the objects. Only for skycut option

#Compute Shot noise as (double): Shot noise factor parameter. See. Gil-Marin et al. 2014 Eq. ..... Only for skycut option. By defalut set it to 1

#Shuffle randoms (no/redshift/radec/both): Whether to create or not a new random catalogue (of the same size of the inputed) using the radec positions of the data and the z of the randoms (redshift); using the redshift positions of the data and the radec of the randoms (radec); usig both radec and redshift from the data (both). By defaul set it to "no". 

#Write shuffled randoms (yes/no): Whether the shuffled random catalogue is written in an output. By default set it to "no". 

#Compute Window Selection function (yes/no): Whether the window function RR counts are computed (for rusticoX 3 different RR couts are performed). This is performed through a brute-force pair-counting and needs paralelization and takes a lot of time. By default set it to "no".

#Bin for window normalization (int) Bin used to nonrmalize W0=1

#DeltaS binning (double) 1.0 Size of the bin for the window. 

#Percentage of randoms selected in % (double) 10. Percentage of randoms selected to perform the RR counts

#Yamamoto aproximation (yes/no): Whether the LOS yamamoto is used (all LOS-mu weight in one member of the pair). 

In case your indput comes from a grid-cell density (and not particles) you must specify the following data

#Value of Pnoise (double) 1/n
#Value of Bnoise1 (double) 1/n
#Value of Bnoise2 (double) 1/n2
#I22 normalization (double): 1.
#I33 normalization (double): 1.

== Reading file==

in order to adapt the code to the format of your files please remove all the headers and edit the read_line.c file. 

====Citation====

If you use this code for your published or unpublished work, please refer it to Gil-Marin, Hector 2020

====Disclaimer====

The author assumes no responsibility or liability for any errors or omissions in the content of this code.  The information contained in this site is provided on an “asis” basis with no guarantees of completeness, accuracy, usefulness or timeliness.
