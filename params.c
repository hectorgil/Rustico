#Main parameters
#Type of code (rustico/rusticoX): rustico
#Type of Box (periodic/periodicFKP/cutsky): periodic
#Type of file (ascii/gadget)x2: ascii
#Type of input (density/particles): particles
#Number of gadget files(int)x2: 1
#RSD distorsion on gadget periodic box (yes/no)x2: yes
#Size of the Box (double/double): 0 2600.
#Type of Computation (DSE/DSY/FFT): FFT

#Power Spectrum options
#Binning for the Power Spectrum (linear/log10): linear
#Size of the bin for the power spectrum (double):  0.01
#k-range for computation (double/double): 0 0.61
#Do anisotropy signal (yes/no): yes
#Do odd multipoles (yes/no): no
#Write kvectors in each bin(yes/no): no
#Do mu-binning Power Spectrum (yes/no): no
#Number of mu-bins (int): 120
#Different files for mu-bin (yes/no): no

#Bispectrum parameters
#Do Bispectrum (yes/no): no
#Do Bispectrum multipoles (yes/no): no
#Do Multigrid (yes/no): yes
#Bispectrum optimization (lowmem/standard/himem): standard
#Triangle Shapes (ALL/EQU/ISO/SQU): ALL
#Size of the bin for the bispectrum (double): 0.0144996584
#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): FFT
#Write triangles in each bin(yes/no): no
#Path for triangles in each bin: ./power_spectra/triangles

#Read inout parameters
#Path of datax2: /home/DATA/hector/challenge_mocks/NSeriesCutsky/Nseries_boxes/BoxN1_w.mock
#Path of randomsx2: none
#Path of output: ./test_weights
#Identifier of output: redshift_space_Dk6_box1_w2
#Write header: yes
#Write density: no

#FFT parameters
#Number of Grid Cells power (int): 9
#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): PCS
#Type of Yamamoto (GridCenter/GridAverage): GridCenter
#Number of interlacing steps (int): 2
#Do Grid Correction? (yes/no): yes

#Cutsky parameters
#Redshift Range (double/double)x2: 0.0 20.0
#Omega matter value (double): 0.307115
#Area effective value in deg^2 (double)x2: 1550
#Smoothing scale for normalization in Mpc/h (double): 1. 
#Quadrupole as (L0L2/L1L1): L0L2
#Octopole as (L0L3/L1L2): L1L2
#Hexadecapole as (L0L4/L2L2/L1L3): L2L2
#Compute Normalization as (area/density): density
#Compute Normalization using (randoms/data): data
#Compute Shot noise as (double): 1.0
#Shuffle randoms (no/redshift/radec/both): no
#Write shuffled randoms (yes/no): no

#Window function paircounts
#Compute Window Selection function (yes/no): no
#Bin for window normalization (int) 1
#DeltaS binning (double) 1.0
#Percentage of randoms selected in % (double) 1.
#Yamamoto aproximation (yes/no): yes

#Density input options
#Value of Pnoise (double) 1000.
#Value of Bnoise1 (double) 1000.
#Value of Bnoise2 (double) 1000.
#I22 normalization (double): 1.
#I33 normalization (double): 1.
