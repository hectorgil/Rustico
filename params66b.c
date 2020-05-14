#Main parameters
#Type of code (rustico/rusticoX): rustico
#Type of Box (periodic/cutsky): cutsky
#Type of file (ascii/gadget)x2: ascii
#Number of gadget files(int)x2: 16
#RSD distorsion on gadget periodic box (yes/no)x2: yes
#Size of the Box (double/double): -1800. +1700.
#Type of Computation (DSE/DSY/FFT): FFT

#Power Spectrum options
#Binning for the Power Spectrum (linear/log10): linear
#Size of the bin for the power spectrum (double):  0.01
#k-range for computation (double/double): 0 0.45
#Do anisotropy signal (yes/no): yes
#Do odd multipoles (yes/no): yes
#Do mu-binning Power Spectrum (yes/no): no
#Number of mu-bins (int): 120
#Different files for mu-bin (yes/no): no

#Bispectrum parameters
#Do Bispectrum (yes/no): no
#Do Bispectrum multipoles (yes/no): no
#Do Multigrid (yes/no): no
#Triangle Shapes (ALL/EQU/ISO/SQU): ALL
#Size of the bin for the bispectrum (double): 0.015707963
#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): FFT
#Write triangles in each bin(yes/no): no
#Path for triangles in each bin: ./power_spectra/triangles

#Read inout parameters
#Path of datax2: /DATA/hector/challenge_mocks/NSeriesCutsky/CutskyN1.rdzw
#Path of randomsx2: test/Randoms_Nseries_1_randomredshift.txt
#Path of output: ./test
#Identifier of output: Nseries_1_randomredshift_check
#Write header: yes

#FFT parameters
#Number of Grid Cells power (int): 8
#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): PCS
#Type of Yamamoto (GridCenter/GridAverage): GridCenter
#Number of interlacing steps (int): 2
#Do Grid Correction? (yes/no): yes

#Cutsky parameters
#Redshift Range (double/double)x2: 0.43 0.70
#Omega matter value (double): 0.31
#Area effective value in deg^2 (double)x2: 7341.000000
#Quadrupole as (L0L2/L1L1): L0L2
#Octopole as (L0L3/L1L2): L1L2
#Hexadecapole as (L0L4/L2L2/L1L3): L2L2
#Compute Normalization as (area/density): area
#Compute Normalization using (randoms/data): data
#Compute Shot noise as (double): 1.0
#Shuffle randoms (no/redshift/radec/both): no
#Write shuffled randoms (yes/no): no

#Window function paircounts
#Compute Window Selection function (yes/no): no
#Bin for window normalization (int) 5
#DeltaS binning (double) 1.0
#Percentage of randoms selected in % (double) 1.
#Yamamoto aproximation (yes/no): yes
