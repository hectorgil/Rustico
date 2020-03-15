/*
RUSTICO Rapid foUrier STatIstics COde
All rights reserved
Author: Hector Gil Marin
Date: 1st Feb 2019
email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "mass_assignment.h"
#include "read_positions.h"
#include "functions.h"
#include "fftw_compute.h"
#include "ps_write.h"
//#include "order_algorithm.h"
#include "bispectrum.h"

typedef struct{
          double OMEGA_M;
}f_params;


int main(int argc, char *argv[])
{

FILE *f,*g;
int i,tid;

char name_file_ini[2000];//name of inizialization file
//Main Parameters read in
double L1,L2;//Box limit in Mpc/h
double kf,kny;
char type_of_survey[50];//type of survey: Periodic or Cutsky
char type_of_computation[10];//type of computation: DSY /  FFT / DSE
int power_grid;//power in number of grid cells: integer between 6 and 15 
int ngrid;//Number of grid cells: pow(2,power_grid)
char binning_type[10];//Type of binning for the power spectrum: linear or log
double bin_ps;//Size of the bin for the power spectrum in h/Mpc
int mubin;//number of bins for mu
char do_mu_bins[10];//perform mu bining (only for periodic boxes)
char do_anisotropy[10];// compute P2,P4 and maybe P1 P3
char do_odd_multipoles[10];//compute P1,P3
char file_for_mu[10];
char header[10];//write header yes/no
char type_of_file[10];//ascii or gadget(only for periodic boxes)
int gadget_files;//number of gadget files per realization
int snapshot_num;//individual gadget file
char RSD[10];//RSD distorsion on Gadget boxes?

//Bispectrum parameters
char do_bispectrum[10];//Do Bispectrum at all?
char do_bispectrum2[10];//Do bispectrum multipoles?
double Deltakbis;//Triangle Bin size in terms of k-fundamental
double kmin,kmax;//Minimum and maximum k
char triangles_num[2000];//FFT,APR_SUM,APR_EFF,EXA_EFF
char write_triangles[2000];
char path_for_triangles[2000];
char triangles_id[2000];
char do_multigrid[10];
char triangle_shapes[10];

//Read inout parameters
char name_data_in[2000];//Path and name of data source
char name_gadget_file[2000];//Full path for gadget like file
long int Ndata;//Number of lines of data source
long int Ndata2;//Number of lines of data used
char name_randoms_in[2000];//Path and name of random source
long int Nrand;//Number of lines of random source
long int Nrand2;//Number of lines of random used
//long int Nrand22;//Number of lines of random used (duplicated)

char name_path_out[2000];//Path where to write output
char name_id[2000];//String to identify output
char name_ps_out[2000];//Final name of the output for the power spectrum
char name_ps_out2[2000];//Final name of the output for the power spectrum for mu bin filess
char name_bs_out[2000];//Final name of the output for the bispectrum
char name_bs002_out[2000];//Final name of the output for the bispectrum002
char name_bs020_out[2000];//Final name of the output for the bispectrum020
char name_bs200_out[2000];//Final name of the output for the bispectrum200

char name_den_out[2000];//Final name of the output for the density
char name_wink_out[2000];//Final name of the output for the power spectrum


//FFT parameters
char type_of_mass_assigment[10];//Type of mass assignment: NGC, CIC, TSC, PCS, P4S, P5S
int Ninterlacing;//Number of interlacing steps (1 for no-interlacing steps)
char grid_correction_string[10];// Do Grid Correction: yes/no
int mode_correction;//Correction factor power
char type_of_yamamoto[20];

//Cutsky parameters
double z_min,z_max;//Minimum and maximum (excluding) redshift cuts
double Omega_m;//Value of Omega matter;
double Area_survey;//Value of Area of the survey in deg^2
char Hexadecapole_type[20];//L4 or L2L2 or L1L3
char Octopole_type[20];//L3 or L1L2
char Quadrupole_type[20];//L2 or L1L1
char type_normalization_mode[20];//Normalize according to the area of the survey or the n(z) value: Area / Density
char type_normalization_mode2[20];//Normalize according to the n(z) of data or randoms file
double Shot_noise_factor;// Factor between 0 and 1. 0 Correspond ....
char shuffle[10];//Generate n(z) for the randoms
char window_function[10];//Compute Window Selection function
char write_shuffled_randoms[10];

//Positions pointers
  double* pos_x;
  double* pos_y;
  double* pos_z;
  double* weight;
  double* radata;
  double* decdata;
  double* zdata;
  double* wcoldata;
  double* wsysdata;
  double* wfkpdata;
  double* nzdata;
 

  double* pos_x_rand;
  double* pos_y_rand;
  double* pos_z_rand;
  double* weight_rand;
//  double DeltaR;  

 //Maximum and minimum values for positions of galaxies and randoms.
  double max,min;
//Parameters relative to shot noise, effective redshifts, number of particles and normalization
double Psn_1a, Psn_1b, Psn_2a,Psn_2b, z_effective_data,I_norm_data,z_effective_rand,I_norm_rand,alpha_data,alpha_data1,alpha_rand,alpha_rand1,alpha,alpha1,I22,I_norm_data2,I_norm_data3,I_norm_data4, I_norm_rand2,I_norm_rand3,I_norm_rand4;
double P_shot_noise1,P_shot_noise2,P_shot_noise,P_shot_noise_win;
double Bsn1a,Bsn2a,Bsn1b,Bsn2b,Bsn1,Bsn2,IN1,IN2,IN11,IN22,I3_norm_data,I3_norm_data2,I3_norm_data3,I3_norm_data4,I3_norm_rand,I3_norm_rand2,I3_norm_rand3,I3_norm_rand4,I33,Bsn,IN;
double num_effective,num_effective_rand;
double num_effective2,num_effective2_rand;
double num_effective3,num_effective3_rand;

    int n_lines_parallel;//Number of parallel threads



//Read Inizialization parameters
sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_survey);
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_file);
fscanf(f,"%*s %*s %*s %*s %d\n",&gadget_files);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",RSD);
printf("== Inizialization Parameters ==\n");
printf("Type of Survey: %s\n",type_of_survey);
printf("Type of File: %s\n",type_of_file);
if(strcmp(type_of_file, "gadget") == 0){printf("Number of gadget files: %d\n",gadget_files);}
if(strcmp(type_of_file, "gadget") == 0){printf("RSD distortion: %s\n",RSD);}
fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&L1,&L2);
printf("Box edges at %lf Mpc/h and %lf Mpc/h; Size of the Box %lf Mpc/h\n",L1,L2,L2-L1);
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_computation);
printf("Type of Computation: %s\n",type_of_computation);


fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n",binning_type);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&bin_ps);
fscanf(f,"%*s %*s %*s %*s %lf %lf\n\n",&kmin,&kmax);
printf("Binning for the Power Spectrum: %s; Size of bin: %lf; k-range %lf < k[h/Mpc] < %lf\n\n",binning_type,bin_ps,kmin,kmax);
fscanf(f,"%*s %*s %*s %*s %s\n",do_anisotropy);
fscanf(f,"%*s %*s %*s %*s %s\n",do_odd_multipoles);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",do_mu_bins);
fscanf(f,"%*s %*s %*s %*s %d\n",&mubin);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",file_for_mu);
if( strcmp(do_mu_bins, "yes") == 0){printf("Mu-binned in %d parts\n",mubin);}
if( strcmp(do_mu_bins, "no") == 0){mubin=1;}

fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",do_bispectrum);
fscanf(f,"%*s %*s %*s %*s %s\n",do_bispectrum2);
fscanf(f,"%*s %*s %*s %s\n",do_multigrid);
fscanf(f,"%*s %*s %*s %s\n",triangle_shapes);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&Deltakbis);
fscanf(f,"%*s %*s %*s %s\n",triangles_num);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",write_triangles);
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n\n",path_for_triangles);

printf("== Bispectrum Parameters ==\n");
printf("Do Bispectrum? %s\n",do_bispectrum);
printf("Do Bispectrum multipoles? %s\n",do_bispectrum2);
if(strcmp(do_bispectrum, "yes") == 0){printf("Bispectrum bins %lf\nMultigrid Computation:%s\nTriangle Shapes: %s\nTriangle normalization: %s\nWrite individual Triangles:%s \n\n",Deltakbis,do_multigrid,triangle_shapes,triangles_num,write_triangles);}

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",name_data_in);
if(strcmp(type_of_file, "ascii") == 0)
{
g=fopen(name_data_in,"r");
if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
fclose(g);
Ndata=countlines(name_data_in);
}

printf("== Read in/out options ==\n");
printf("Reading data file %s; %ld lines\n",name_data_in,Ndata);

fscanf(f,"%*s %*s %*s %s\n",name_randoms_in);

if(strcmp(type_of_survey, "cutsky") == 0)
{
g=fopen(name_randoms_in,"r");
if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_randoms_in);return 0;}
fclose(g);

Nrand=countlines(name_randoms_in);
printf("Reading randoms file %s; %ld lines\n",name_randoms_in,Nrand);
}

if(strcmp(type_of_file, "gadget") == 0)
{ 
  Ndata=0;
  for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {
     sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);
     g=fopen(name_gadget_file,"r");
     if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
     fclose(g);
     Ndata+=count_particles_gadget(name_gadget_file);
  }
P_shot_noise=pow(L2-L1,3)/Ndata*1.;
}

fscanf(f,"%*s %*s %*s %s\n",name_path_out);
printf("Output files at %s\n",name_path_out);
fscanf(f,"%*s %*s %*s %s\n",name_id);
printf("Output Id %s\n",name_id);
fscanf(f,"%*s %*s %s\n",header);
printf("Write header? %s\n\n",header);
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %d\n",&power_grid);
ngrid=pow(2,power_grid);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",type_of_mass_assigment);
printf("== FFT options ==\n");
printf("Number of k-modes and grid cells per side: %d\n",ngrid);
printf("Type of mass assingment %s\n",type_of_mass_assigment);

if(strcmp(type_of_mass_assigment, "NGC") == 0){mode_correction=1;}
if(strcmp(type_of_mass_assigment, "CIC") == 0){mode_correction=2;}
if(strcmp(type_of_mass_assigment, "TSC") == 0){mode_correction=3;}
if(strcmp(type_of_mass_assigment, "PCS") == 0){mode_correction=4;}
if(strcmp(type_of_mass_assigment, "P4S") == 0){mode_correction=5;}
if(strcmp(type_of_mass_assigment, "P5S") == 0){mode_correction=6;}
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_yamamoto);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Type of Yamamoto: %s\n",type_of_yamamoto);}
fscanf(f,"%*s %*s %*s %*s %*s %d\n",&Ninterlacing);
printf("Number of Interlacing steps %d\n",Ninterlacing);
fscanf(f,"%*s %*s %*s %*s %s\n",grid_correction_string);
printf("Do grid correction? %s\n\n",grid_correction_string);
if(strcmp(grid_correction_string, "no") == 0){mode_correction=0;}
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %lf %lf\n",&z_min,&z_max);
if(strcmp(type_of_survey, "cutsky") == 0){printf("== Cutsky options ==\n");}
if(strcmp(type_of_survey, "cutsky") == 0){printf("Redshift cuts: %lf < z < %lf\n",z_min,z_max);}
fscanf(f,"%*s %*s %*s %*s %lf\n",&Omega_m);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Omega_m: %lf\n",Omega_m);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %lf\n",&Area_survey);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Area of the survey %lf deg2\n",Area_survey);}
fscanf(f,"%*s %*s %*s %s\n",Quadrupole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Quadrupole as %s\n",Quadrupole_type);}
fscanf(f,"%*s %*s %*s %s\n",Octopole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles, "yes") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Octopole as %s\n",Octopole_type);}
fscanf(f,"%*s %*s %*s %s\n",Hexadecapole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Hexadecapole as %s\n",Hexadecapole_type);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s\n",type_normalization_mode);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode2);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s file n(z)\n",type_normalization_mode2);}
fscanf(f,"%*s %*s %*s %*s %*s %lf\n",&Shot_noise_factor);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Shot noise factor set to %lf\n",Shot_noise_factor);}
fscanf(f,"%*s %*s %*s %s\n",shuffle);
if( strcmp(shuffle, "no") ==0 ){printf("Random catalogue provided used.\n");}
if( strcmp(shuffle, "radec") ==0 ){printf("RA-dec shuffling.\n");}
if( strcmp(shuffle, "redshift") ==0 ){printf("Redshift shuffling\n");}
if( strcmp(shuffle, "both") ==0 ){printf("Both RA-dec & redshift shuffling (random catalogue only used for setting the size of the new random catalogue).\n");}

fscanf(f,"%*s %*s %*s %*s %s\n",write_shuffled_randoms);
if( strcmp(write_shuffled_randoms, "yes") == 0 ||  strcmp(write_shuffled_randoms, "no") ==0 ){printf("Write Shuffled randoms: %s\n",write_shuffled_randoms);}

fscanf(f,"%*s %*s %*s %*s %*s %s\n",window_function);
if( strcmp(window_function, "yes") == 0 ||  strcmp(window_function, "no") ==0 ){printf("Compute window function: %s\n",window_function);}
fclose(f);

//Error conditions
if(strcmp(type_of_file, "gadget") != 0 && strcmp(type_of_file, "ascii") != 0){printf("File type must be either 'gadget' or 'ascii'. Entry read %s. Exiting now...\n",type_of_file);return 0;}
if( strcmp(type_of_survey, "cutsky") != 0 && strcmp(type_of_survey, "periodic") != 0){printf("Survey type must be either 'cutsky' or 'periodic'. Entry read %s. Exiting now...\n",type_of_survey);return 0;}
if( strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_file, "gadget") == 0 ){printf("Warning. Cutsky+gadget option not available. Exiting now...\n");return 0;}
if(gadget_files<1){printf("Warning. gadget files entry must be >0. Entered value %d. Exiting now...\n",gadget_files);return 0;}
if( strcmp(RSD, "yes") != 0 && strcmp(RSD, "no") != 0 ){printf("Warning. RSD option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",RSD);return 0;}
if(L2<=L1){printf("Error: L2 parameter has to be large then L1. Exiting now...\n");return 0;}
if( strcmp(type_of_computation, "DSE") !=0 &&  strcmp(type_of_computation, "DSY") !=0 &&  strcmp(type_of_computation, "FFT") !=0){printf("Type of computation only accepts 'DSE', 'DSY' or 'FFT' options. Entry read %s. Exiting now...\n",type_of_computation);return 0;}
if( strcmp(type_of_computation, "FFT") !=0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. Bispectrum computation only accepts FFT computation option. Entry read %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(binning_type, "log10") !=0 && strcmp(binning_type, "linear") !=0){printf("Warning. Binning type must be either 'log10' or 'linear'. Entry read %s. Exiting now...\n",binning_type);return 0;}
if( strcmp(binning_type, "log10") ==0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. 'log10' type of binning not available for bispectrum computation on this version. Exiting now...\n"); return 0;}
if(bin_ps<=0){printf("Error: Size of the power spectrum bin has to be greater than 0: Entry read %lf. Exiting now...\n",bin_ps);return 0;}
if(kmin<0 || kmax<=0 || kmin>kmax){printf("Warning: Unusual values for maximum and/or minimum k-values for the bispectrum computation: kmin=%lf, kmax=%lf. Exiting now...\n",kmin, kmax);return 0;}
if(kmin==0 && strcmp(binning_type, "log10") == 0){printf("Cannot set kmin=0 and log-k binning. Exiting now...\n");return 0;}
if( strcmp(do_bispectrum, "yes") != 0 && strcmp(do_bispectrum, "no") != 0){printf("Error. Bispectrum entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(do_multigrid, "yes") != 0 && strcmp(do_multigrid, "no") != 0){printf("Error. Multigrid entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_multigrid);return 0;}
if( strcmp(do_multigrid, "yes") == 0 && Ninterlacing<2 ){printf("Warning. Multigrid option requires a number of interlacing steps >1. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"NGC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"CIC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"TSC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(triangle_shapes,"ALL") !=0 && strcmp(triangle_shapes,"EQU") !=0 && strcmp(triangle_shapes,"ISO") !=0 && strcmp(triangle_shapes,"SQU") !=0 ){printf("Error. Triangle shapes entry only accepts 'ALL', 'ISO', 'EQU' or 'SQU'. Read entry %s. Exiting now...\n",triangle_shapes);return 0;}
if(Deltakbis<=0){printf("Error: Size of the bispectrum bin has to be greater than 0: %lf kf. Exiting now...\n",Deltakbis);return 0;}
if( strcmp(triangles_num, "FFT") != 0 && strcmp(triangles_num, "APR_SUM") !=0 && strcmp(triangles_num, "EXA_SUM") !=0  && strcmp(triangles_num, "APR_EFF") !=0 && strcmp(triangles_num, "EXA_EFF") !=0){printf("Error. Number of Triangles normalization option only accepts: 'FFT', 'APR_SUM', 'EXA_SUM', 'APR_EFF', 'EXA_EFF'. Entry read %s. Exiting now...\n",triangles_num);return 0;}
if( strcmp(triangles_num, "EXA_SUM") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_SUM' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(triangles_num, "EXA_EFF") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_EFF' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(write_triangles, "yes") != 0 &&  strcmp(write_triangles, "no") !=0){printf("Error. write triangle entry must be either 'yes' or 'now'. Exiting now...\n");return 0;}
if( strcmp(write_triangles, "yes") == 0 && strcmp(triangle_shapes,"SQU") !=0){printf("Waring. Write triangle option 'yes' is only recomended for squeezed triangle shapes 'SQU'. Exiting now...\n");return 0;}
if( strcmp(header, "yes") !=0 && strcmp(header, "no") !=0){printf("Waring. Write header option must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",header);return 0;}
if(power_grid<4 || power_grid>15){printf("Warning: Unusual value for number of grid cells per side: 2^%d=%d. Exiting now...\n",power_grid,ngrid);return 0;}
if(strcmp(type_of_mass_assigment,"NGC") !=0 && strcmp(type_of_mass_assigment,"CIC") !=0 && strcmp(type_of_mass_assigment,"TSC") !=0 && strcmp(type_of_mass_assigment,"PCS") !=0 && strcmp(type_of_mass_assigment,"P4S") !=0 &&  strcmp(type_of_mass_assigment,"P5S") !=0){printf("Error. Type of mass assigment must be either 'NGC', 'CIC', 'TSC', 'PCS', 'P4S' or 'P5S'. Entry read %s. Exiting now...\n",type_of_mass_assigment);return 0;}
if( strcmp(type_of_yamamoto, "GridCenter") != 0 && strcmp(type_of_yamamoto, "GridAverage") != 0){printf("Error. Type of Yamamoto option must be either 'GridCenter' or 'GridAverage'. Entry read %s. Exiting now...\n",type_of_yamamoto);return 0;}
if(Ninterlacing<=0){printf("Error: Number of interglacing steps has to be equal or larger than 1. %d\n",Ninterlacing);return 0;}
if( strcmp(grid_correction_string, "yes") !=0 && strcmp(grid_correction_string, "no") !=0){printf("Grid correction input must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",grid_correction_string);return 0;}

if(strcmp(type_of_survey, "cutsky") == 0){
if(z_min>=z_max){printf("Error. Minimum value for redshift is larger than the maximum: z_min=%lf; z_max=%lf. Exiting now...\n",z_min,z_max);return 0;}
if(Omega_m<=0 || Omega_m>1){printf("Warning. Unusual value for Omega_m, Omega_m=%lf. Exiting now...\n",Omega_m);return 0;}
if(Area_survey<=0){printf("Warning. Usual value for the Area of the survey: %lf. Exiting now...\n",Area_survey);return 0;}
if( strcmp(Hexadecapole_type, "L4") !=0 && strcmp(Hexadecapole_type, "L2L2") !=0 && strcmp(Hexadecapole_type, "L1L3") !=0){printf("Hexadecapole option must be either 'L2L2' or 'L4' or 'L1L3'. Entry read %s. Exiting now...\n",Hexadecapole_type);return 0;}
if( strcmp(Octopole_type, "L3") !=0 && strcmp(Octopole_type, "L1L2") !=0){printf("Octopole option must be either 'L1L2' or 'L3'. Entry read %s. Exiting now...\n",Octopole_type);return 0;}
if( strcmp(Quadrupole_type, "L2") !=0 && strcmp(Quadrupole_type, "L1L1") !=0){printf("Quadrupole option must be either 'L1L1' or 'L2'. Entry read %s. Exiting now...\n",Quadrupole_type);return 0;}

if( strcmp(type_normalization_mode, "area") !=0 && strcmp(type_normalization_mode, "density") !=0){printf("Error. Normalisation type must be either 'area' or 'density'. Entry read %s. Exiting now...\n",type_normalization_mode);return 0;}
if( strcmp(type_normalization_mode2, "data") !=0 && strcmp(type_normalization_mode2, "randoms") !=0){printf("Error. Normalisation type must be either 'data' or 'randoms'. Entry read %s. Exiting now...\n",type_normalization_mode2);return 0;}
if(Shot_noise_factor>1 || Shot_noise_factor<0){printf("Warning. Usual value for the Shot noise factor: %lf. Exiting now...\n",Shot_noise_factor);return 0;}
if( strcmp(shuffle, "radec") != 0 && strcmp(shuffle, "no") != 0 && strcmp(shuffle, "redshift") != 0 && strcmp(shuffle, "both") != 0){printf("Warning. Shuffling option not reconised: %s. Exiting now...\n",shuffle);return 0;}
if( strcmp(window_function, "yes") != 0 && strcmp(window_function, "no") != 0 ){printf("Warning. Window function option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",window_function);return 0;}
if( strcmp(write_shuffled_randoms, "yes") != 0 && strcmp(write_shuffled_randoms, "no") != 0 ){printf("Warning. Write shuffled randoms option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",write_shuffled_randoms);return 0;}

if( strcmp(do_bispectrum, "no") == 0 && strcmp(do_bispectrum2, "yes") == 0 ){printf("Warning. Bispectrum multipoles computation requires a bispectrum calculation. Exiting now...\n");return 0;}

if( strcmp(do_mu_bins, "yes") != 0 &&  strcmp(do_mu_bins, "no") != 0){printf("Warning. Mu-binning option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_mu_bins);return 0;}

if( strcmp(do_anisotropy, "yes") != 0 &&  strcmp(do_anisotropy, "no") != 0){printf("Warning. Anisotropy option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_odd_multipoles);return 0;}

if( strcmp(do_odd_multipoles, "yes") != 0 &&  strcmp(do_odd_multipoles, "no") != 0){printf("Warning. Odd multipole option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_odd_multipoles);return 0;}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(type_of_survey, "periodic") != 0){printf("Warning. Mu-binning option only available for non-varying line-of-sight samples. Exiting now...\n");return 0;}

if(mubin<0){printf("Warning. The numer of mu-bins has to be a positive integer. Exiting now...\n");return 0;}

if( strcmp(do_mu_bins, "yes") == 0){
if( strcmp(file_for_mu, "yes") != 0 &&  strcmp(file_for_mu, "no") != 0){printf("Warning. Different files when Mu-binning option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_mu_bins);return 0;}
}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(do_odd_multipoles, "yes") == 0){printf("Warning. Cannot do odd multipoles if mu-binning option is also selected. Enable either mu-binning or odd-multipoles. Exiting now...\n");return 0;}

if( strcmp(do_anisotropy, "no") == 0 &&  strcmp(do_odd_multipoles, "yes") == 0){printf("Warning. Cannot do odd multipoles if anisotropy option is not selected. Enable either anisotropy or disable odd-multipoles. Exiting now...\n");return 0;}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(do_anisotropy, "no") == 0){printf("Warning. Cannot do mu-binning if anisotropy option is not selected. Enable either anisotropy or disable mu-binning. Exiting now...\n");return 0;}

if( strcmp(type_of_computation, "DSE") ==0 && strcmp(do_anisotropy, "no") == 0){printf("Warning. There's no point of doing an DSE computation for an isotropic signal, as this is equivalent to a DSY or FFT calculation. Exiting now...\n");exit(0);}

//if( strcmp(do_bispectrum, "yes") == 0 && strcmp(type_normalization_mode, "density") == 0 ){printf("Warning. Bispectrum computation requires a normalization by 'area' and not by 'density'. Exiting now...\n");return 0;}
}
//etc strings.....

//Determine number of processors available for openmp  
        #pragma omp parallel for private(i,tid) shared(n_lines_parallel,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){n_lines_parallel=omp_get_num_threads();}
        }
        i=fftw_init_threads();
        fftw_plan_with_nthreads(n_lines_parallel);
        printf("Number of processors used: %d\n",n_lines_parallel);


//Reading files.
double parameter_value[28];
parameter_value[0]=Omega_m;
parameter_value[1]=z_min;
parameter_value[2]=z_max;
parameter_value[3]=Ndata;
parameter_value[13]=Area_survey;

if(strcmp(type_of_survey, "cutsky") == 0)
{
                Ndata2=get_number_used_lines_data(name_data_in,parameter_value);if(Ndata2<1000){printf("Warning, unusual low value for Ndata=%ld\n",Ndata2);if(Ndata2==0){exit(0);}}
alpha_data1=parameter_value[12];
}
if(strcmp(type_of_survey, "periodic") == 0)
{
Ndata2=Ndata;
}

if(strcmp(type_of_file, "ascii") == 0)//these are only kept stored during all the process for ascii files. Gadget files keep the name of the file and read it each time they need
{
                pos_x = (double*) calloc(Ndata2, sizeof(double));
		pos_y = (double*) calloc(Ndata2, sizeof(double));
	        pos_z = (double*) calloc(Ndata2, sizeof(double));
		weight = (double*) calloc(Ndata2, sizeof(double));
}

if(strcmp(type_of_survey, "cutsky") == 0)
{
//pos_x,pos_y,pos_z,weight are loaded with the position of particles
//Ndata is uploaded to the number of particles used
//Psn_1a,Psn_1b,Psn_2a,Psn_2b are uploaded with information on the shot noise
//z_efffective is uploaded with the effective redshift of the sample
//num_effective is uploaded with the effective number of particles
//I_norm is uploaded with information relative to the normalization based on density
//alpha is uploeaded with information relative to the effective ratio between data and randoms
printf("Reading %s...",name_data_in);

                radata = (double*) calloc(Ndata2, sizeof(double));
                decdata = (double*) calloc(Ndata2, sizeof(double));
                zdata = (double*) calloc(Ndata2, sizeof(double));
                wcoldata = (double*) calloc(Ndata2, sizeof(double));
                wsysdata = (double*) calloc(Ndata2, sizeof(double));
                wfkpdata = (double*) calloc(Ndata2, sizeof(double));
                nzdata = (double*) calloc(Ndata2, sizeof(double));
                

get_skycuts_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value,type_normalization_mode,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,shuffle);


Ndata2=parameter_value[3];//doesn't change
//for(i=0;i<Ndata2;i++){if(zdata[i]>2.2 || zdata[i]<0.8){printf("%lf\n",zdata[i]);}}exit(0);
Psn_1a=parameter_value[4];
Psn_1b=parameter_value[5];
I3_norm_data=parameter_value[6];
z_effective_data=parameter_value[7];
num_effective=parameter_value[8];
I_norm_data=parameter_value[9];
min=parameter_value[10];
max=parameter_value[11];
alpha_data=parameter_value[12];

I_norm_data2=parameter_value[14];
I_norm_data3=parameter_value[15];
I_norm_data4=parameter_value[16];
num_effective2=parameter_value[17];
num_effective3=parameter_value[13];

I3_norm_data2=parameter_value[18];
I3_norm_data3=parameter_value[19];
I3_norm_data4=parameter_value[20];

Bsn1a=parameter_value[21];
Bsn1b=parameter_value[22];
IN1=parameter_value[23];
IN2=parameter_value[24];
IN11=parameter_value[26];
IN22=parameter_value[27];

parameter_value[3]=Ndata;
sprintf(name_den_out,"%s/Density_galaxies_%s.txt",name_path_out,name_id);
get_skycuts_write_density_data(name_data_in, parameter_value,name_den_out);
parameter_value[3]=Ndata2;

printf("Ok!\n");

parameter_value[0]=Omega_m;
parameter_value[1]=z_min;
parameter_value[2]=z_max;
parameter_value[3]=Nrand;


if( strcmp(shuffle, "no") == 0 )
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
}
else
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
Nrand2=Nrand;
}
alpha_rand1=parameter_value[12];
alpha1=alpha_data1/alpha_rand1;


        pos_x_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand2, sizeof(double));
        weight_rand = (double*) calloc(Nrand2, sizeof(double));

/*
Nrand22=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand22<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand22);}
parameter_value[4]=Nrand22;

        pos_x_rand = (double*) calloc(Nrand22, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand22, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand22, sizeof(double));
        weight_rand = (double*) calloc(Nrand22, sizeof(double));
*/
								
printf("Reading %s...",name_randoms_in);
parameter_value[12]=alpha_data;
sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);

get_skycuts_randoms(name_path_out,name_id,name_randoms_in, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, parameter_value,type_normalization_mode, type_normalization_mode2,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,Ndata2,alpha1,shuffle,write_shuffled_randoms);

Nrand2=parameter_value[3];
Psn_2a=parameter_value[4];
Psn_2b=parameter_value[5];
z_effective_rand=parameter_value[7];
num_effective_rand=parameter_value[8];
num_effective2_rand=parameter_value[20];
num_effective3_rand=parameter_value[13];

if(min>parameter_value[10]){min=parameter_value[10];}
if(max<parameter_value[11]){max=parameter_value[11];}
alpha_rand=parameter_value[12];

I_norm_rand=parameter_value[9];
I_norm_rand2=parameter_value[14];
I_norm_rand3=parameter_value[15];
I_norm_rand4=parameter_value[16];

I3_norm_rand=parameter_value[6];
I3_norm_rand2=parameter_value[17];
I3_norm_rand3=parameter_value[18];
I3_norm_rand4=parameter_value[19];
    
    Bsn2a=parameter_value[21];
    Bsn2b=parameter_value[22];

alpha=alpha_data/alpha_rand;
I_norm_rand=I_norm_rand*alpha;
I3_norm_rand=I3_norm_rand*alpha;

    if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "density") == 0 ){I22=I_norm_rand;I33=I3_norm_rand;}
    if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "density") == 0){I22=I_norm_data;I33=I3_norm_data;}

if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "area") == 0 ){I22=I_norm_rand2;I33=I3_norm_rand2;}
if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "area") == 0){I22=I_norm_data2;I33=I3_norm_data2;}


P_shot_noise1=(Psn_1a+alpha*alpha*Psn_2a)/I22;
P_shot_noise2=(Psn_1b+alpha*alpha*Psn_2b)/I22;
P_shot_noise=P_shot_noise1*Shot_noise_factor+P_shot_noise2*(1.-Shot_noise_factor);
P_shot_noise_win=alpha*alpha*(Psn_2a*Shot_noise_factor+Psn_2b*(1.-Shot_noise_factor))/I22;
Bsn1=(Bsn1a-alpha*alpha*alpha*Bsn2a)/I33;
Bsn2=(Bsn1b-alpha*alpha*alpha*Bsn2b)/I33;
    
Bsn=Bsn1*Shot_noise_factor+(1.-Shot_noise_factor)*Bsn2;
if( strcmp(type_normalization_mode, "area") == 0 ){IN=(IN1*Shot_noise_factor+IN2*(1.-Shot_noise_factor))/I33;}
if( strcmp(type_normalization_mode, "density") == 0 ){IN=(IN11*Shot_noise_factor+IN22*(1.-Shot_noise_factor))/I33;}

parameter_value[3]=Nrand;
//parameter_value[4]=Nrand22;

sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);
get_skycuts_write_density_randoms(name_randoms_in, parameter_value,alpha,name_den_out,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,Ndata2,alpha1,shuffle);
//parameter_value[3]=Nrand2;//useless?

free(radata);
free(decdata);
free(zdata);
free(wcoldata);
free(wsysdata);
free(wfkpdata);
free(nzdata);
printf("Ok!\n\n");


}

if(strcmp(type_of_survey, "periodic") == 0 && strcmp(type_of_file, "ascii") == 0)
{
printf("Reading %s...",name_data_in);
get_periodic_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value);
P_shot_noise=pow(L2-L1,3)/Ndata*1.;
min=parameter_value[10];
max=parameter_value[11];
printf("Ok!\n\n");
}

if(max>L2 || min<L1){printf("Warning: Limits of the box are exceeded by the data or random galaxies: data particles found at the limits %lf and %lf. Exiting now...\n",min,max);return 0;}

if(strcmp(window_function, "yes") == 0){sprintf(name_wink_out,"%s/Windowk_%s.txt",name_path_out,name_id);}

if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){sprintf(name_ps_out,"%s/Power_Spectrum_%s",name_path_out,name_id);}
else{sprintf(name_ps_out,"%s/Power_Spectrum_%s.txt",name_path_out,name_id);}


sprintf(name_bs_out,"%s/Bispectrum_%s.txt",name_path_out,name_id);

sprintf(name_bs002_out,"%s/Bispectrum_Quadrupole002_%s.txt",name_path_out,name_id);
sprintf(name_bs020_out,"%s/Bispectrum_Quadrupole020_%s.txt",name_path_out,name_id);
sprintf(name_bs200_out,"%s/Bispectrum_Quadrupole200_%s.txt",name_path_out,name_id);

sprintf(triangles_id,"%s/Triangles_%s",path_for_triangles,name_id);

if(strcmp(window_function, "yes") == 0){

f=fopen(name_wink_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}

if( strcmp(header, "yes") == 0)
{
fprintf(f,"# Keff \t Kav \t W0 \t W2 \t W4 \t W6 \t W8\n");
}
fclose(f);

}

for(i=0;i<mubin;i++){

sprintf(name_ps_out2,"%s_%d.txt",name_ps_out,i);
if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){f=fopen(name_ps_out2,"w");}
else{f=fopen(name_ps_out,"w");}

if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_ps_out);return 0;}

//write header for the power spectrum file
if( strcmp(header, "yes") == 0)
{
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s\n", type_normalization_mode );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I22);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Shot noise factor %lf\n",Shot_noise_factor);}
fprintf(f,"#Shot noise value %lf\n",P_shot_noise);
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0 ){fprintf(f,"#Quadrupole as %s\n",Quadrupole_type);}
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles,"yes") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Octopole as %s\n",Octopole_type);}
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Hexadecapole as %s\n",Hexadecapole_type);}

if(strcmp(do_anisotropy,"yes") == 0){
if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"no") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Quadrupole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Dipole\t Quadrupole\t Octopole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
if(strcmp(do_mu_bins, "yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t mu-centerbin\t mu-eff\t P(k,mu)-Pshotnoise\t number of modes\t Pshotnoise\n");}
}else{
fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t number of modes\t Pshotnoise\n");
}

}
fclose(f);

}

//Write header for the Bispectrum file
if( strcmp(header, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0)
{
f=fopen(name_bs_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-eff\t k1-centerbin\t k2-eff\t k2-centerbin\t k3-eff\t k3-centerbin\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

fclose(f);
}
if( strcmp(header, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 &&  strcmp(do_bispectrum2, "yes") == 0)
{
f=fopen(name_bs002_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs002_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-eff\t k1-centerbin\t k2-eff\t k2-centerbin\t k3-eff\t k3-centerbin\t Bispectrum002-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);
////
f=fopen(name_bs020_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs020_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-eff\t k1-centerbin\t k2-eff\t k2-centerbin\t k3-eff\t k3-centerbin\t Bispectrum020-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);
///
f=fopen(name_bs200_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs200_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-eff\t k1-centerbin\t k2-eff\t k2-centerbin\t k3-eff\t k3-centerbin\t Bispectrum200-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);

}

//Start Power Spectrum Computation for Cutsky
printf("== Computing the Power Spectrum ==\n");
/*
//Determine number of processors available for openmp  
        #pragma omp parallel for private(i,tid) shared(n_lines_parallel,ngrid)
	for(i=0;i<ngrid;i++)
	{
		tid=omp_get_thread_num();
		if(tid==0 && i==0){n_lines_parallel=omp_get_num_threads();}
	}
	i=fftw_init_threads();
	fftw_plan_with_nthreads(n_lines_parallel);
	printf("Number of processors used: %d\n",n_lines_parallel);
*/
//Compute and write the Power Spectrum for FFT+Skycut type of survey.
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0)
{
parameter_value[0]=L1;
parameter_value[1]=L2;
check_box_for_yamamoto(parameter_value,ngrid);

L1=parameter_value[0];
L2=parameter_value[1];

if(strcmp(type_of_yamamoto, "GridCenter") == 0){

if(strcmp(window_function, "yes") == 0){
kf=2.*(4.*atan(1.))/(L2-L1);
kny=2.*(4.*atan(1.))/(L2-L1)*ngrid/2.;
loop_interlacing_skycut(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise_win, kf, I22, alpha, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_wink_out, type_of_mass_assigment,do_bispectrum,1);
}

loop_interlacing_skycut(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise, bin_ps, I22, alpha, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out, type_of_mass_assigment,do_bispectrum,0);

}

if(strcmp(type_of_yamamoto, "GridAverage") == 0){

if(strcmp(window_function, "yes") == 0){
kf=2.*(4.*atan(1.))/(L2-L1);
kny=2.*(4.*atan(1.))/(L2-L1)*ngrid/2.;
loop_interlacing_skycut2(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise_win, kf, I22, alpha, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_wink_out, type_of_mass_assigment,do_bispectrum, 1);
}

loop_interlacing_skycut2(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise, bin_ps, I22, alpha, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out, type_of_mass_assigment,do_bispectrum, 0);


}

}

if(strcmp(type_of_computation, "DSY") == 0 || strcmp(type_of_computation, "DSE") == 0)
{

     if( strcmp(type_of_survey, "cutsky") == 0)
     {
         loop_directsum_yamamoto_skycut_caller(kmin,kmax,pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, P_shot_noise, bin_ps, I22, alpha, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out,type_of_computation);         
     }
     if(strcmp(type_of_survey, "periodic") == 0)
     {
          printf("No Direct Sum for periodic box at the moment. Exiting now...\n");return 0;
     }

}

//Compute and write the Power Spectrum for FFT+Box with constant line of sight along z.
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
if(strcmp(type_of_file, "ascii") == 0)
{
loop_interlacing_periodic(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight, Ndata, L1, L2, ngrid, P_shot_noise, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out, type_of_mass_assigment,do_odd_multipoles,do_anisotropy,do_bispectrum,mubin,file_for_mu);
}
if(strcmp(type_of_file, "gadget") == 0)
{
loop_interlacing_periodic_gadget(kmin,kmax,Ninterlacing, name_data_in ,gadget_files, L1, L2, ngrid, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out, type_of_mass_assigment,Shot_noise_factor,grid_correction_string,RSD,do_odd_multipoles,do_anisotropy,mubin,file_for_mu);
}
}

if(strcmp(do_bispectrum, "no") == 0){
printf("Computation of Power Spectrum finished sucessfully!\n\n");
return 0;
}

printf("== Computing the Bispectrum ==\n");

if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0)
{
//write bispectrum header for cutsky
loop_bispectrum_skycut_caller(kmin, kmax, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, ngrid, P_shot_noise, Deltakbis, I33,I22, IN, Bsn, alpha, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, triangle_shapes,do_bispectrum2);
}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
if(strcmp(type_of_file, "ascii") == 0)
{
//write bispectrum header for periodic
loop_bispectrum_skycut_caller(kmin, kmax, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, ngrid, P_shot_noise, Deltakbis, 0,0, 0, 0, 0, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, triangle_shapes,do_bispectrum2);
}
if(strcmp(type_of_file, "gadget") == 0)
{
loop_bispectrum_periodic_for_gadget_caller(kmin,kmax,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_data_in,gadget_files, do_multigrid, triangle_shapes,RSD,do_bispectrum2);

}

}
if(strcmp(type_of_computation, "FFT") != 0)
{
printf("No direct Sum for the Bispectrum yet!\n");
return 0;
}



//Improve: sorting algorithm. Avoid system calls

printf("Computation of the Bispectrum finished sucessfully!\n\n");
return 0;
}
