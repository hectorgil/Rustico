#Main parameters
#Type of Box (periodic/cutsky): periodic
#Type of file (ascii/gadget): gadget
#Number of gadget files(int): 16
#RSD distorsion on gadget periodic box (yes/no): yes
#Size of the Box (double/double): 0 +2400.
#Type of Computation (DSE/DSY/FFT): FFT

#Power Spectrum options
#Binning for the Power Spectrum (linear/log10): linear
#Size of the bin for the power spectrum (double):  0.01
#k-range for computation (double/double): 0 0.33
#Do anisotropy signal (yes/no): yes
#Do odd multipoles (yes/no): yes
#Do mu-binning Power Spectrum (yes/no): no
#Number of mu-bins (int): 120
#Different files for mu-bin (yes/no): no

#Bispectrum parameters
#Do Bispectrum (yes/no): yes
#Do Bispectrum multipoles (yes/no): yes
#Do Multigrid (yes/no): no
#Triangle Shapes (ALL/EQU/ISO/SQU): ALL
#Size of the bin for the bispectrum (double): 0.015707963
#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): FFT
#Write triangles in each bin(yes/no): no
#Path for triangles in each bin: ./power_spectra/triangles

#Read inout parameters
#Path of data: /DATA/cwagner/newcov/1/snapshot_004
#Path of randoms: none
#Path of output: ./test
#Identifier of output: z0_Dk6_run1
#Write header: yes

#FFT parameters
#Number of Grid Cells power (int): 8
#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): P5S
#Type of Yamamoto (GridCenter/GridAverage): GridCenter
#Number of interlacing steps (int): 2
#Do Grid Correction? (yes/no): yes

#Cutsky parameters
#Redshift Range (double/double): 0.8 2.2
#Omega matter value (double): 0.31
#Area effective value in deg^2 (double): 2843.514974043143
#Quadrupole as (L2/L1L1): L2
#Octopole as (L3/L1L2): L1L2
#Hexadecapole as (L4/L2L2/L1L3): L2L2
#Compute Normalization as (area/density): density
#Compute Normalization using (randoms/data): data
#Compute Shot noise as (double): 1.0
#Shuffle randoms (no/redshift/radec/both): no
#Write shuffled randoms (yes/no): no
#Compute Window Selection function (yes/no): no
