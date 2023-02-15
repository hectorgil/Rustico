#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "cubature.h"
#include "read_line.h"
#include "functions.h"
//#include "structures.h"


long int get_number_used_lines_periodic(char *filename, double parameter_value[],int type)
{
FILE *f;
long int npar=(long int)(parameter_value[3]);
long int i;
double weight;
int veto;
long int npar_used=0;
double params[9];
f=fopen(filename,"r");
for(i=0;i<npar;i++)
{

get_line_periodic(f, params,type);
weight=params[3];
veto=(int)(params[7]);

if(weight!=0 && veto>0){npar_used++;}

}
fclose(f);
return npar_used;
}

double get_number_used_lines_weighted_periodic(char *filename, double parameter_value[],int type)
{
FILE *f;
long int npar=(long int)(parameter_value[3]);
long int i;
double weight;
int veto;
double npar_used=0;
double params[9];
f=fopen(filename,"r");
//printf("%ld\n",npar);
for(i=0;i<npar;i++)
{

get_line_periodic(f, params,type);
weight=(double)(params[3]);
veto=(int)(params[7]);

if(weight!=0 && veto>0){npar_used=npar_used+weight;}

}
fclose(f);
return npar_used;
}


long int get_number_used_lines_data(char *filename, double parameter_value[])
{
FILE *f;
long int npar=(long int)(parameter_value[3]);
double z_min=parameter_value[1];
double z_max=parameter_value[2];
long int i;
int j;
double  weight_col;
double weight_sys;
double redshift;
long int npar_used=0;
double alpha=0;
double distance;
int veto;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;
double MIN[1], MAX[1],radial,nuissance;
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;

f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
get_line(f, params,0);
redshift=params[2];
distance=params[8];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);

if(redshift==-1 && distance==-1){printf("Error, line %ld in data returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_number_used_lines_data\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_number_used_lines_data\n",redshift);}
}

if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
{
npar_used++;
alpha=alpha+weight_col*weight_sys;
}
}
fclose(f);

if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}

parameter_value[12]=alpha;
return npar_used;
}

long int get_number_used_lines_randoms(char *filename, double parameter_value[])
{
FILE *f;
double redshift,distance;
long int npar=(long int)(parameter_value[3]);
double z_min=parameter_value[1];
double z_max=parameter_value[2];
long int i;
int j;
double weight_col,weight_sys;
int veto;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;
long int npar_used=0;
double alpha=0;
double MIN[1], MAX[1],radial,nuissance;
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;

f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
get_line(f, params,1);
veto=(int)(params[7]);
weight_col=params[5];
weight_sys=params[6];
redshift=params[2];
distance=params[8];

if(redshift==-1 && distance==-1){printf("Error, line %ld in data returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0
}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_number_used_lines_randoms\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_number_used_lines_randoms\n",redshift);}
}

if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
{
npar_used++;
alpha=alpha+weight_col*weight_sys;
}

}
fclose(f);

if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}

parameter_value[12]=alpha;
return npar_used;
}

void get_skycuts_write_density_data(char *filename, double parameter_value[], char *name_den_out)
{
double Pi=(4.*atan(1.));
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;

double MIN[1];
double MAX[1];
double r_min,r_max;
long int npar;
double weight_col,wnoz,wcp;
long int i;
int j;
FILE *f,*g;
double redshift,weight_fkp,weight_sys, radial,nuissance,distance;
npar=(long int)(parameter_value[3])-3;
double Area;
long int n_bin_r;
double DeltaR=parameter_value[25];
long int index_radial;
double *radial_cell, *radial_all_weight_cell,*radial_fkp_cell, *radial_weight_cell;
double *z_cell;
int veto,obj;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;

Area=parameter_value[13]*pow(Pi/180.,2);

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);

if(Omega_m+Omega_L==1){r_min=r_min*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_min=sin(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_min=sinh(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

if(Omega_m+Omega_L==1){r_max=r_max*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_max=sin(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_max=sinh(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

  n_bin_r=(long int)((r_max-r_min)/DeltaR);//printf("\n %lf %d\n",DeltaR,n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_all_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  z_cell= (double*) calloc(n_bin_r, sizeof(double));

f=fopen(filename,"r");
for(i=0;i<npar;i++)
{

get_line(f, params,0);
redshift=params[2];
weight_fkp=params[3];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);
distance=params[8];

if(redshift==-1 && distance==-1){printf("Error, line %ld in data returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_data\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_data\n",redshift);}
}


if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
  {

  MAX[0]=redshift;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

      index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("\n Error bins data1 (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%ld) %ld %lf. Exiting now...\n",radial,redshift,r_min,r_max,i,n_bin_r,DeltaR);exit(0);}
      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;
      radial_weight_cell[index_radial]+=weight_col*weight_sys;
      radial_all_weight_cell[index_radial]+=weight_col*weight_sys*weight_fkp;
      z_cell[index_radial]+=redshift;

}

}
printf("\nWriting %s...",name_den_out);
        f=fopen(name_den_out,"w");
        if(f==NULL){printf("File %s could not be created. Exiting now...\n",name_den_out);exit(0);}
        fprintf(f,"#Interval: %lf Mpc/h, Area= %lf deg^2\n",DeltaR,Area*pow(Pi/180.,-2));
        fprintf(f,"# z <nobs> <wc nobs> <wc wfkp nobs> N_objects\n");
        for(i=0;i<n_bin_r;i++)
        {
                if(radial_cell[i]!=0)
                {
                        z_cell[i]=z_cell[i]/radial_cell[i];
                        obj=(int)(radial_cell[i]);
                        fprintf(f,"%lf %e %e %e %d\n",z_cell[i],radial_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.),  radial_weight_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.), radial_all_weight_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.),obj );
                }
        }
        fclose(f);
printf("Ok!\n");
free(radial_cell);
free(radial_all_weight_cell);
free(radial_fkp_cell);
free(radial_weight_cell);
free(z_cell);
free(function_parameters);
if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}
}


void get_skycuts_data(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[], char *type_normalization_mode, double radata[], double decdata[], double zdata[], double wcoldata[], double wsysdata[],double wfkpdata[], double nzdata[],char *shuffle,double Rsmooth)
{
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
double z_min=parameter_value[1];
double z_max=parameter_value[2];

f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;

double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
long int npar;
double weight_col;
long int i;
int j;
FILE *f;
double RA,dec,redshift,weight_fkp,weight_sys,n_z,theta, radial,nuissance,max,min,distance;
npar=(long int)(parameter_value[3]);
double Area;
long int n_bin_r;
double r_min,r_max;
Area=parameter_value[13]*pow(Pi/180.,2);
double normalization,normalizationbis;
double zeff;
double num,numzeff;
double num2,num3;
long int npar_used;
double Psn_1a;
double Psn_1b;
double Bsn_1a,Bsn_1b;
double alpha_data;
double I22,I22_w_data,I22_w_data_will;
double I33,I33_w_data,I33_w_data_will;
double IN1,IN2;
double DeltaR;
double DeltaR_min;
long int index_radial;
double *radial_cell, *radial_all_weight_cell,*radial_fkp_cell, *radial_weight_cell;
double *radial_weight_cell_N1, *radial_weight_cell_N2, *z_cell;
int i_DeltaR;
double I22_min,I22_w_data_min,I22_w_data_will_min;
double I33_min,I33_w_data_min,I33_w_data_will_min;
double IN1_min,IN2_min,IN11,IN22;
int veto;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;

I22_min=0;
I22_w_data_min=0;
I22_w_data_will_min=0;
I33_min=0;
I33_w_data_min=0;
I33_w_data_will_min=0;
IN1_min=0;
IN2_min=0;

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);

if(Omega_m+Omega_L==1){r_min=r_min*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_min=sin(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_min=sinh(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

if(Omega_m+Omega_L==1){r_max=r_max*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_max=sin(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_max=sinh(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

i_DeltaR=0;
do
{
normalization=0;
normalizationbis=0;
zeff=0;
num=0;
num2=0;
num3=0;
numzeff=0;
npar_used=0;
Psn_1a=0;
Psn_1b=0;
Bsn_1a=0;
Bsn_1b=0;
IN11=0;
IN22=0;
alpha_data=0;

i_DeltaR++;
DeltaR=i_DeltaR*0.5;
if(Rsmooth>0){DeltaR=Rsmooth;}

  n_bin_r=(long int)((r_max-r_min)/DeltaR);//printf("\n %lf %d\n",DeltaR,n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_all_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_N1= (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_N2= (double*) calloc(n_bin_r, sizeof(double));
  z_cell= (double*) calloc(n_bin_r, sizeof(double));

  max=-9999999;
  min=9999999;
f=fopen(filename,"r");

for(i=0;i<npar;i++)
{
get_line(f, params,0);

RA=params[0];
dec=params[1];
redshift=params[2];
weight_fkp=params[3];
n_z=params[4];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);
distance=params[8];

if(redshift==-1 && distance==-1){printf("Error, line %ld in data returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_data\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_data\n",redshift);}
}

theta=90.-dec;
if(redshift>z_min && redshift<z_max  && veto==1 && weight_col>0)
  {

//              normalization+=n_z*pow(weight_fkp*weight_sys*weight_col,2);//effective normalization from data. correct if nz provided is nz with no weights
              normalization+=n_z*weight_fkp*weight_fkp*weight_sys*weight_col;//effective normalization from data. correct if nz provided is wcol*wsys weighted
//              normalizationbis+=n_z*n_z*pow(weight_fkp*weight_sys*weight_col,3);//effective normalization from data. correct if nz provided is nz with no weights
              normalizationbis+=n_z*n_z*weight_sys*weight_col*pow(weight_fkp,3);//effective normalization from data. correct if nz provided is nz provided is wcol*wsys weighted


  zeff+=redshift*pow(weight_col*weight_sys*weight_fkp,2);
  num+=weight_col*weight_sys;//Effective number of objects with wsys (real number)
  numzeff+=pow(weight_col*weight_sys*weight_fkp,2);//normalizes zeff
  num2+=weight_col;//Effective number of objects without wsys
  num3+=weight_col*weight_sys*weight_fkp;//Effective number of objects without wsys

 radata[npar_used]=RA;
 decdata[npar_used]=dec; 
 zdata[npar_used]=redshift;
 wcoldata[npar_used]=weight_col;
 wsysdata[npar_used]=weight_sys;
 wfkpdata[npar_used]=weight_fkp;
 nzdata[npar_used]=n_z;

  //from deg to radiants
  RA=RA*Pi/180.;
  dec=dec*Pi/180.;
  theta=theta*Pi/180.;

  //Determination of comoving distance given the redshifts and the value of omegamatter
  MAX[0]=redshift;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

      index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("\n Error bins data2 (radial) radial=%lf z=(%lf < %lf < %lf)  r=(%lf,%lf) (line=%ld) %ld %lf. Om=%lf. Exiting now...\n",radial,z_min,redshift,z_max,r_min,r_max,i,n_bin_r,DeltaR,Omega_m);exit(0);}   

      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;
      radial_weight_cell[index_radial]+=weight_col*weight_sys;
      radial_all_weight_cell[index_radial]+=weight_col*weight_sys*weight_fkp;
      radial_weight_cell_N1[index_radial]+=pow(weight_col*weight_sys,2);
      radial_weight_cell_N2[index_radial]+=weight_col*pow(weight_sys,2);
      z_cell[index_radial]+=redshift;


  //From polar to cartesian coordinates
  pos_x[npar_used]=radial*sin(theta)*cos(RA);
  pos_y[npar_used]=radial*sin(theta)*sin(RA);
  pos_z[npar_used]=radial*cos(theta);
 
  weight[npar_used]=weight_fkp*weight_col*weight_sys;//total weight


  Psn_1a+=pow(weight_fkp*weight_col*weight_sys,2);//false pairs
  Psn_1b+=weight_fkp*weight_col*weight_sys*weight_fkp*weight_sys;// true pairs
  Bsn_1a+=pow(weight_fkp*weight_col*weight_sys,3);//false pairs
  Bsn_1b+=pow(weight_fkp*weight_sys,3)*weight_col;//true pairs
  IN11+=n_z*pow(weight_fkp,3)*pow(weight_col*weight_sys,2);//false pairs
  IN22+=n_z*pow(weight_fkp,3)*pow(weight_sys,2)*weight_col;//true pairs

  if(pos_x[npar_used]>max || npar_used==0){max=pos_x[npar_used];}
  if(pos_y[npar_used]>max){max=pos_y[npar_used];}
  if(pos_z[npar_used]>max){max=pos_z[npar_used];}
  if(pos_x[npar_used]<min || npar_used==0){min=pos_x[npar_used];}
  if(pos_y[npar_used]<min){min=pos_y[npar_used];}
  if(pos_z[npar_used]<min){min=pos_z[npar_used];}
 
 alpha_data+=weight_col*weight_sys*weight_fkp;
  npar_used++;

  }
  
  }
fclose(f);
//exit(0);
//for(i=0;i<10;i++){printf("%e %e %e\n",radata[i],decdata[i],zdata[i]);}

//if(strcmp(shuffle, "both") == 0){reshuffle(radata,decdata,wcoldata,wsysdata,npar_used);}//this re-shuffles the indices on angular-data wrt the unchanged radial-data in case of both-shuffling option.

//for(i=0;i<npar_used;i++){printf("%e %e %e\n",radata[i],decdata[i],zdata[i]);}
//exit(0);

        I22=0;I33=0;
        I22_w_data=0;I33_w_data=0;
        I22_w_data_will=0;I33_w_data_will=0;
        IN1=0;IN2=0;

        for(i=0;i<n_bin_r;i++)
        {
         I22+=Area*pow(radial_all_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
         I33+=Area*pow(radial_all_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

          if(radial_cell[i]!=0)
          {
//Normalization of the Power Spectrum
             I22_w_data+=Area*pow(radial_fkp_cell[i]/radial_cell[i]*1.,2)*pow(radial_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
             I22_w_data_will+=Area*pow(radial_all_weight_cell[i]/radial_cell[i]*1.,2)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

//Normalization of the Bispectrum
             I33_w_data+=Area*pow(radial_fkp_cell[i]/radial_cell[i]*1.,3)*pow(radial_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

             I33_w_data_will+=Area*pow(radial_all_weight_cell[i]/radial_cell[i]*1.,3)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

   
//Shot noise of the Bispectrum
                        IN2+=Area*pow(radial_fkp_cell[i]/radial_cell[i]*1.,3)*(radial_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*(radial_weight_cell_N2[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

                        IN1+=Area*pow(radial_fkp_cell[i]/radial_cell[i]*1.,3)*(radial_weight_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*(radial_weight_cell_N1[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;


                }
        }

free(radial_cell);
free(radial_all_weight_cell);
free(radial_fkp_cell);
free(radial_weight_cell);
free(radial_weight_cell_N1);
free(radial_weight_cell_N2);
free(z_cell);


if(i_DeltaR==1)//
{
I22_min=I22;
I22_w_data_min=I22_w_data;
I22_w_data_will_min=I22_w_data_will;
I33_min=I33;
I33_w_data_min=I33_w_data;
I33_w_data_will_min=I33_w_data_will;
IN1_min=IN1;
IN2_min=IN2;
DeltaR_min=i_DeltaR;
}

if(I22<I22_min){I22_min=I22;DeltaR_min=DeltaR;}
if(I22_w_data<I22_w_data_min){I22_w_data_min=I22_w_data;}
if(I22_w_data_will<I22_w_data_will_min){I22_w_data_will_min=I22_w_data_will;}

if(I33<I33_min){I33_min=I33;}
if(I33_w_data<I33_w_data_min){I33_w_data_min=I33_w_data;}
if(I33_w_data_will<I33_w_data_will_min){I33_w_data_will_min=I33_w_data_will;}

if(IN1<IN1_min){IN1_min=IN1;}
if(IN2<IN2_min){IN2_min=IN2;}

//printf("%lf %lf %.10lf %lf %lf\n",DeltaR,I22,I33,IN1,IN2);
//printf("%lf %lf %.10lf %lf %lf\n",DeltaR,I22_min,I33_min,IN1_min,IN2_min);


}while(DeltaR<40 && strcmp(type_normalization_mode, "area") == 0 && Rsmooth==0);

if(strcmp(type_normalization_mode, "density") == 0 && Rsmooth==0){DeltaR_min=10.;}//Returns a fixed value if normalized by density
if(strcmp(type_normalization_mode, "density") == 0 && Rsmooth>0){DeltaR_min=Rsmooth;}//Returns a fixed value if normalized by density

free(function_parameters);
//Copy needed information 
parameter_value[3]=npar_used*1.;
parameter_value[4]=Psn_1a;
parameter_value[5]=Psn_1b;
parameter_value[6]=normalizationbis;
parameter_value[7]=zeff/numzeff*1.;
parameter_value[8]=num2;
parameter_value[9]=normalization;
parameter_value[10]=min;
parameter_value[11]=max;
parameter_value[12]=alpha_data;
parameter_value[28]=num3;

parameter_value[14]=I22_min;
parameter_value[15]=I22_w_data_min;
parameter_value[16]=I22_w_data_will_min;
parameter_value[17]=num;
parameter_value[18]=I33_min;
parameter_value[19]=I33_w_data_min;
parameter_value[20]=I33_w_data_will_min;
parameter_value[21]=Bsn_1a;
parameter_value[22]=Bsn_1b;
parameter_value[23]=IN1_min;
parameter_value[24]=IN2_min;
parameter_value[25]=DeltaR_min;
parameter_value[26]=IN11;
parameter_value[27]=IN22;

if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}

}



void get_skycuts_write_density_randoms(char *filename, double parameter_value[], double alpha, char *name_den_out, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata, long int Ndata,double alpha_true, char *shuffle)
{
long int *L;
L =  (long int*) calloc(Ndata, sizeof(long int));
long int interval;

long int seed;
double P0;
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;
double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
double Area=parameter_value[13]*pow(Pi/180.,2);
double DeltaR=parameter_value[25];
long int i,npar,j,l,lmax,ieff,i2;//,npar2,npar3;
long int random;
double randomd;
FILE *f;
double redshift,weight_fkp, radial,nuissance,weight_col,weight_sys,weight_noz,weight_cp,distance;
double *radial_cell,*radial_fkp_cell,*z_cell,*radial_weight_cell,*radial_all_weight_cell;
double r_min,r_max;
long int n_bin_r;
long int index_radial;
int veto,obj;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;

npar=(long int)(parameter_value[3]);

//npar3=(long int)(parameter_value[3]);//total number of randoms
//npar2=(long int)(parameter_value[4]);//randoms which pass the z and wcol cuts
//if(strcmp(shuffle, "no") == 0 || strcmp(shuffle,"redshift") ==  0){npar=npar3;} //no, redshift
//if(strcmp(shuffle, "both") == 0 || strcmp(shuffle,"radec") ==  0){npar=npar2;} //both, radec

//long int chunks=(long int)(npar/Ndata*1.);
//printf("npar=%ld, Ndata=%ld\n",npar,Ndata);

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);

if(Omega_m+Omega_L==1){r_min=r_min*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_min=sin(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_min=sinh(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

if(Omega_m+Omega_L==1){r_max=r_max*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_max=sin(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_max=sinh(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

  n_bin_r=(long int)((r_max-r_min)/DeltaR);//printf("\n %d\n",n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_all_weight_cell = (double*) calloc(n_bin_r, sizeof(double));
  z_cell = (double*) calloc(n_bin_r, sizeof(double));

    

f=fopen(filename,"r");

if(strcmp(shuffle, "no") == 0)
{

for(i=0;i<npar;i++)
{

get_line(f, params,1);

redshift=params[2];
weight_fkp=params[3];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);
distance=params[8];

if(redshift==-1 && distance==-1){printf("Error, line %ld in randoms returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_randoms\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_randoms\n",redshift);}
}if(redshift==-1 && distance==-1){printf("Error, line %ld in randoms returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_randoms\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_randoms\n",redshift);}
}


        if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
        {
            MAX[0]=redshift;
            MIN[0]=0;
            adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

 index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins randoms1 (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%ld) %ld %lf\n",radial,redshift,r_min,r_max,i,n_bin_r,DeltaR);}
  
            radial_cell[index_radial]+=1.;
            radial_fkp_cell[index_radial]+=weight_fkp;
            radial_weight_cell[index_radial]+=weight_col*weight_sys;
            radial_all_weight_cell[index_radial]+=weight_col*weight_sys*weight_fkp;
            z_cell[index_radial]+=redshift;
            
        }
    
}


}

if(strcmp(shuffle, "radec") == 0 || strcmp(shuffle,"redshift") == 0 || strcmp(shuffle, "both") == 0)
{
seed=(long int)(radata[0]*decdata[0]*100000);//deterministic seed
//srand48(seed*2);
if(strcmp(shuffle, "both") == 0)
{
printf("\nShuffling RA-dec wrt redshifts...\n");
reshuffle(radata,decdata,wcoldata,wsysdata,Ndata,seed);
printf("Done\n");
}

printf("Shuffling data-catalogue ordering...\n");

shuffle_lines(Ndata,L,2*seed);
/*
i=0;
do
{
randomd=drand48();
random=(long int)(randomd*Ndata);
if(random>=Ndata || random<0){random=0;}
L[i]=random;
interval=1;

if(random==i){interval=0;}
else{

      if(i>0)
      {
          for(j=0;j<i;j++)//this can be parallelized
          {
            if(L[i]==L[j]){interval=0;break;}
          }
      }
}

i=i+interval;
}while(i<Ndata);
printf("Done\n");
*/


i=0;
j=0;
for(l=0;l<npar;l++)
{

        get_line(f, params,1); 
        redshift=params[2];
        weight_fkp=params[3];
        weight_col=params[5]; 
        weight_sys=params[6]; 
        veto=(int)(params[7]);
        distance=params[8];

if(redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_write_density_randoms\n",redshift);}
}

        
        if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
        {


i2=L[i];

if( strcmp(shuffle, "radec") == 0){

        
        redshift=zdata[i2];//get z from stored data
//        weight_fkp=wfkpdata[i2];//get fkp from data
        P0=(1./wfkpdata[i2]-1.)/nzdata[i2];
        weight_fkp=1./(1+nzdata[i2]*P0);

//        weight_col=1;//wcoldata[i];
//        weight_sys=1;//wsysdata[i];        
      
}


if( strcmp(shuffle, "redshift") == 0){
       
        weight_col=wcoldata[i2];
        weight_sys=wsysdata[i2];        
}


if( strcmp(shuffle, "both") == 0){

        //ieff=i2+j;
        ieff=i+j;
        if(ieff>=Ndata){ieff=ieff-Ndata;}
        ieff=L[ieff];
        weight_col=wcoldata[ieff];
        weight_sys=wsysdata[ieff];
        redshift=zdata[i2];//get z from stored data
//        weight_fkp=wfkpdata[i2];//get fkp from data
        P0=(1./wfkpdata[i2]-1.)/nzdata[i2];
        weight_fkp=1./(1+nzdata[i2]*P0);
        
}

i++;
if(i==Ndata){i=0;j++;}


            MAX[0]=redshift;
            MIN[0]=0;
            adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

       index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins randoms1 (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%ld) %ld %lf\n",radial,redshift,r_min,r_max,i,n_bin_r,DeltaR);}

            radial_cell[index_radial]+=1.;
            radial_fkp_cell[index_radial]+=weight_fkp;
            radial_weight_cell[index_radial]+=weight_col*weight_sys;
            radial_all_weight_cell[index_radial]+=weight_col*weight_sys*weight_fkp;
            z_cell[index_radial]+=redshift;
}

}


}



 fclose(f);



    printf("\nWriting %s...",name_den_out);
    f=fopen(name_den_out,"w");
    if(f==NULL){printf("File %s could not be created. Exiting now...\n",name_den_out);exit(0);}
fprintf(f,"#Interval: %lf Mpc/h, Area= %lf deg^2\n",DeltaR,Area*pow(Pi/180.,-2));
fprintf(f,"# z <nobs> <wc nobs> <wc wfkp nobs> N_obj\n");
for(i=0;i<n_bin_r;i++)
{
    if(radial_cell[i]!=0)
    {
        z_cell[i]=z_cell[i]/radial_cell[i];
        obj=(int)(radial_cell[i]);
        fprintf(f,"%lf %e %e %e %d\n",z_cell[i],radial_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.),  radial_weight_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.), radial_all_weight_cell[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.),obj );
    }
}

fclose(f);
printf("Ok!\n");
free(radial_cell);
free(radial_all_weight_cell);
free(radial_fkp_cell);
free(radial_weight_cell);
free(z_cell);
free(function_parameters);
free(L);

if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}

}

void get_skycuts_randoms(char *path, char *id, char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],char *type_normalization_mode, char *type_normalization_mode2, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata,long int Ndata, double alpha_true, char *shuffle, char *write_shuffled_randoms,char *name_out_randoms,double Rsmooth)
{

long int *L;
L =  (long int*) calloc(Ndata, sizeof(long int));
long int interval;
double P0;
long int seed;
double Omega_m=parameter_value[0];
double Omega_L=parameter_value[29];
double speed_c=parameter_value[30];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
//char name_out_randoms[200];
//sprintf(name_out_randoms,"%s/Randoms_%s.dat",path,id);

f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
(*function_parameters).OMEGA_L=Omega_L;

double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
double Area=parameter_value[13]*pow(Pi/180.,2);

long int i,i2,j,l,lmax,npar,ieff;//,npar2,npar3;
FILE *f,*g;
double RA,dec,redshift,distance,redshift_aux,weight_fkp,weight_col,weight_col_aux,weight_cp,weight_noz,weight_sys,n_z,theta, radial,nuissance,max,min;

npar=(long int)(parameter_value[3]);//total number of randoms 

//npar3=(long int)(parameter_value[3]);//total number of randoms
//npar2=(long int)(parameter_value[4]);//randoms which pass the z and wcol cuts
//npar=0;
//if(strcmp(shuffle, "no") == 0 || strcmp(shuffle,"redshift") ==  0){npar=npar3;} //no, redshift
//if(strcmp(shuffle, "both") == 0 || strcmp(shuffle,"radec") ==  0){npar=npar2;} //both, radec

//long int chunks;

//chunks=(long int)(npar/Ndata*1.);

//printf("npar2=%ld,npar3=%ld, Ndata=%ld\n",npar2,npar3,Ndata);

double normalization,normalizationbis;
double zeff,num,num2,num3,numeff;
long int npar_used;
double alpha_data;
double alpha;//alpha count from data

double I22_w_randoms;
double I33_w_randoms;
double DeltaR;
double DeltaR_min;
long int index_radial;
double *radial_cell,*radial_fkp_cell;
double r_min,r_max;
int i_DeltaR;
double I22_w_randoms_min;
double I33_w_randoms_min;
long int n_bin_r;
int veto;
double Psn_2a,Psn_2b,Bsn_2a,Bsn_2b;
long int random;
double randomd;
double params[9];
double zlow=-1;
double zhigh=20;
int N_z_table=10000;
int free_table_dist_z=0;
double *redshift_table,*dist_table;    

I22_w_randoms_min=0;
I33_w_randoms_min=0;

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);

if(Omega_m+Omega_L==1){r_min=r_min*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_min=sin(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_min=sinh(r_min*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

if(Omega_m+Omega_L==1){r_max=r_max*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){r_max=sin(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){r_max=sinh(r_max*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

i_DeltaR=0;
do
{
normalization=0;
normalizationbis=0;
zeff=0;
num=0;
num2=0;
num3=0;
numeff=0;
npar_used=0;
alpha_data=0;
alpha=parameter_value[12];

i_DeltaR++;
DeltaR=i_DeltaR*0.5;
if(Rsmooth>0){DeltaR=Rsmooth;}

  n_bin_r=(long int)((r_max-r_min)/DeltaR);//printf("\n %d\n",n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  max=-9999999;
  min=9999999;
    Psn_2a=0;
    Psn_2b=0;
    Bsn_2a=0;
    Bsn_2b=0;


f=fopen(filename,"r");

if(strcmp(shuffle, "no") == 0)
{

for(i=0;i<npar;i++)
{

get_line(f, params,1);

RA=params[0];
dec=params[1];
redshift=params[2];
weight_fkp=params[3];
n_z=params[4];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);
distance=params[8];

if(redshift==-1 && distance==-1){printf("Error, line %ld in randoms returns no physical redshift or distance. Exiting now...\n",i);exit(0);}

if(i==0 && redshift==-1 && distance!=-1)
{//generate table n interpolate

redshift_table = (double*) calloc(N_z_table, sizeof(double));
dist_table = (double*) calloc(N_z_table, sizeof(double));
free_table_dist_z=1;
for(j=0;j<N_z_table;j++)
{
redshift_table[j]=(zhigh-(zlow))/N_z_table*1.*(j-0.5);//interpolation between z= zhigh and zlow
MAX[0]=redshift_table[j];
MIN[0]=0;
adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);
//dist_table[j]=radial;

if(Omega_m+Omega_L==1){dist_table[j]=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){dist_table[j]=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){dist_table[j]=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

}

redshift=Interpol(distance, dist_table, redshift_table, N_z_table);
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_randoms\n",redshift);}
}

if(i>0 && redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_randoms\n",redshift);}
}


	theta=90.-dec;
	if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
        {

//              normalization+=n_z*pow(weight_fkp*weight_sys*weight_col,2);//effective normalization from data. correct if nz provided is nz with no weights
              normalization+=n_z*weight_fkp*weight_fkp*weight_sys*weight_col;//effective normalization from data. correct if nz provided is wcol*wsys weighted
//              normalizationbis+=n_z*n_z*pow(weight_fkp*weight_sys*weight_col,3);//effective normalization from data. correct if nz provided is nz with no weights
              normalizationbis+=n_z*n_z*weight_sys*weight_col*pow(weight_fkp,3);//effective normalization from data. correct if nz provided is nz provided is wcol*wsys weighted
              zeff+=redshift*pow(weight_col*weight_sys*weight_fkp,2);
              num+=weight_col*weight_sys;//Effective number of objects with wsys (real number)
              numeff+=pow(weight_col*weight_sys*weight_fkp,2);//Effective number of objects with wsys (real number)
              num2+=weight_col;//Effective number of objects without wsys
              num3+=weight_col*weight_sys*weight_fkp;
             
			//From deg to rad
			RA=RA*Pi/180.;
			dec=dec*Pi/180.;
			theta=theta*Pi/180.;

            //Determination of comoving distance given the redshifts and the value of omegamatter
            MAX[0]=redshift;
            MIN[0]=0;
	        adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

           //From polar to cartesian coordinates
		   pos_x[npar_used]=radial*sin(theta)*cos(RA);
		   pos_y[npar_used]=radial*sin(theta)*sin(RA);
		   pos_z[npar_used]=radial*cos(theta);

//printf("%e %e %e %e %e %e %e\n",RA,dec,redshift,weight_fkp,n_z,weight_col,weight_sys);

      index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins randoms22 (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%ld) %ld %lf\n",radial,redshift,r_min,r_max,i,n_bin_r,DeltaR);}
      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;

              weight[npar_used]=weight_fkp*weight_col*weight_sys;

		    if(pos_x[npar_used]>max || npar_used==0){max=pos_x[npar_used];}
			if(pos_y[npar_used]>max){max=pos_y[npar_used];}
			if(pos_z[npar_used]>max){max=pos_z[npar_used];}
			if(pos_x[npar_used]<min || npar_used==0){min=pos_x[npar_used];}
			if(pos_y[npar_used]<min){min=pos_y[npar_used];}
			if(pos_z[npar_used]<min){min=pos_z[npar_used];}

			alpha_data+=weight_fkp*weight_col*weight_sys;
            Psn_2a+=pow(weight_fkp*weight_col*weight_sys,2);//false pairs
            Psn_2b+=weight_fkp*weight_col*weight_sys*weight_fkp*weight_sys;// true pairs
            Bsn_2a+=pow(weight_fkp*weight_col*weight_sys,3);//false pairs
            Bsn_2b+=pow(weight_fkp*weight_sys,3)*weight_col;//true pairs

					    npar_used++;
						  }

	  }
}

if(strcmp(shuffle, "radec") == 0 || strcmp(shuffle, "redshift") == 0 || strcmp(shuffle, "both") == 0)
{
seed=(long int)(radata[0]*decdata[0]*100000);//deterministic seed
//srand48(2*seed);

if(strcmp(shuffle, "both") == 0)
{
printf("\nShuffling RA-dec wrt redshifts...\n");
reshuffle(radata,decdata,wcoldata,wsysdata,Ndata,seed);
printf("Done\n");

}

printf("Shuffling data-catalogue ordering...\n");


shuffle_lines(Ndata,L,2*seed);
/*
i=0;
do
{
randomd=drand48();
random=(long int)(randomd*Ndata);
if(random>=Ndata || random<0){random=0;}
L[i]=random;
//if(random==Ndata){printf("Error 222\n");exit(0);}
interval=1;

if(random==i){interval=0;}
else{

      if(i>0)
      {
          for(j=0;j<i;j++)
          {
            if(L[i]==L[j]){interval=0;break;}
          }
      }
}

i=i+interval;
}while(i<Ndata);
printf("Done\n");
*/


if(strcmp(write_shuffled_randoms, "yes") == 0){g=fopen(name_out_randoms,"w");fclose(g);}
if(strcmp(write_shuffled_randoms, "yes") == 0){g=fopen(name_out_randoms,"a");}


i=0;
j=0;
for(l=0;l<npar;l++)
{

     //i2=0;
     //for(i=0;i2<Ndata;i++)
     //{

get_line(f, params,1);
RA=params[0];
dec=params[1];
redshift=params[2];
weight_fkp=params[3];
n_z=params[4];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);
distance=params[8];

if(redshift==-1 && distance!=-1){//interpolate
redshift=Interpol(distance, dist_table, redshift_table, N_z_table);//z-dist monotonically growing
if(redshift<zlow || redshift>zhigh){printf("Warning, extrapolation made when converting dist into z (z=%lf). Please enlarge zhigh-zlow range in read_positions.c, get_skycuts_randoms\n",redshift);}
}
        
        if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
        {

        i2=L[i];

//if(L[i]<0 || L[i]>=Ndata){printf("Error 111b %ld %ld\n",i,L[i]);exit(0);}


if( strcmp(shuffle, "radec") == 0){

        
        theta=90.-dec;
        redshift=zdata[i2];//get z from stored data
//       weight_fkp=wfkpdata[i2];//get fkp from data
        P0=(1./wfkpdata[i2]-1.)/nzdata[i2];
        weight_fkp=1./(1+nzdata[i2]*P0);
    //    weight_col=1;//wcoldata[i];
    //    weight_sys=1;//wsysdata[i];
        n_z=nzdata[i2]/alpha_true;

}


if( strcmp(shuffle, "redshift") == 0){

        RA=radata[i2];
        dec=decdata[i2];
        theta=90.-dec;
        weight_col=wcoldata[i2];
        weight_sys=wsysdata[i2];

}


if( strcmp(shuffle, "both") == 0){

        //ieff=i2+j;
        ieff=i+j;
        if(ieff>=Ndata){ieff=ieff-Ndata;}
        ieff=L[ieff];
        RA=radata[ieff];
        dec=decdata[ieff];
        theta=90.-dec;
        weight_col=wcoldata[ieff];
        weight_sys=wsysdata[ieff];
        redshift=zdata[i2];//get z from stored data
        weight_fkp=wfkpdata[i2];//get fkp from data
        P0=(1./wfkpdata[i2]-1.)/nzdata[i2]; 
        weight_fkp=1./(1+nzdata[i2]*P0);
        n_z=nzdata[i2]/alpha_true;
        
}

//if(i2==Ndata){printf("Error 333, chunk=%ld\n",j);exit(0);}
//printf("%e %e %e %e %e\n",radata[i2],decdata[i2],zdata[i2],wfkpdata[i2],nzdata[i2]);

i++;
if(i==Ndata){i=0;j++;}


//              normalization+=n_z*pow(weight_fkp*weight_sys*weight_col,2);//effective normalization from data. correct if nz provided is nz with no weights
              normalization+=n_z*weight_fkp*weight_fkp*weight_sys*weight_col;//effective normalization from data. correct if nz provided is wcol*wsys weighted
//              normalizationbis+=n_z*n_z*pow(weight_fkp*weight_sys*weight_col,3);//effective normalization from data. correct if nz provided is nz with no weights
              normalizationbis+=n_z*n_z*weight_sys*weight_col*pow(weight_fkp,3);//effective normalization from data. correct if nz provided is nz provided is wcol*wsys weighted


    zeff+=redshift*pow(weight_col*weight_sys*weight_fkp,2);
    num+=weight_col*weight_sys;//Effective number of objects with wsys (real number)
    numeff+=pow(weight_col*weight_sys*weight_fkp,2);//Effective number of objects with wsys (real number)
    num2+=weight_col;//Effective number of objects without wsys
    num3+=weight_col*weight_sys*weight_fkp;

if(strcmp(write_shuffled_randoms, "yes") == 0){fprintf(g,"%e %e %e %e %e %e %e\n",RA,dec,redshift,weight_fkp,n_z,weight_col,weight_sys);}

                        RA=RA*Pi/180.;
                        dec=dec*Pi/180.;
                        theta=theta*Pi/180.;


     MAX[0]=redshift;
     MIN[0]=0;
     adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

if(Omega_m+Omega_L==1){radial=radial*speed_c/100.;}//flat
if(1-Omega_m-Omega_L<0){radial=sin(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok<0
if(1-Omega_m-Omega_L>0){radial=sinh(radial*pow(fabs(1-Omega_m-Omega_L),0.5))*speed_c/100.*pow(fabs(1-Omega_m-Omega_L),-0.5);}//Ok>0

                   pos_x[npar_used]=radial*sin(theta)*cos(RA);
                   pos_y[npar_used]=radial*sin(theta)*sin(RA);
                   pos_z[npar_used]=radial*cos(theta);

      index_radial=(long int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins randoms23 (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%ld) %ld %lf\n",radial,redshift,r_min,r_max,i,n_bin_r,DeltaR);}
      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;

              weight[npar_used]=weight_fkp*weight_col*weight_sys;

                      if(pos_x[npar_used]>max || npar_used==0){max=pos_x[npar_used];}
                        if(pos_y[npar_used]>max){max=pos_y[npar_used];}
                        if(pos_z[npar_used]>max){max=pos_z[npar_used];}
                        if(pos_x[npar_used]<min || npar_used==0){min=pos_x[npar_used];}
                        if(pos_y[npar_used]<min){min=pos_y[npar_used];}
                        if(pos_z[npar_used]<min){min=pos_z[npar_used];}

                        alpha_data+=weight_fkp*weight_col*weight_sys;
            Psn_2a+=pow(weight_fkp*weight_col*weight_sys,2);//false pairs
            Psn_2b+=weight_fkp*weight_col*weight_sys*weight_fkp*weight_sys;// true pairs
            Bsn_2a+=pow(weight_fkp*weight_col*weight_sys,3);//false pairs
            Bsn_2b+=pow(weight_fkp*weight_sys,3)*weight_col;//true pairs

                                            npar_used++;
}

     }

//if(strcmp(shuffle, "both") == 0){printf("chunk %ld/%ld\n",j,chunks);reshuffle(radata,decdata,wcoldata,wsysdata,Ndata);}//this re-shuffles the indices on angular-data wrt the unchanged radial-data in case of both-shuffling option.

if(strcmp(write_shuffled_randoms, "yes") == 0){fclose(g);}


}


fclose(f);


alpha*=1./alpha_data;
I22_w_randoms=0;
I33_w_randoms=0;

        for(i=0;i<n_bin_r;i++)
        {         
                if(radial_cell[i]!=0)
                {
                   I22_w_randoms+=pow(radial_fkp_cell[i]*1./radial_cell[i]*1.,2)*Area*pow(alpha,2)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
                   I33_w_randoms+=pow(radial_fkp_cell[i]*1./radial_cell[i]*1.,3)*Area*pow(alpha,3)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

                }
         }

if(i_DeltaR==1)//
{
I22_w_randoms_min=I22_w_randoms;
I33_w_randoms_min=I33_w_randoms;
DeltaR_min=i_DeltaR;
}

if(I22_w_randoms<I22_w_randoms_min){I22_w_randoms_min=I22_w_randoms;DeltaR_min=DeltaR;}
if(I33_w_randoms<I33_w_randoms_min){I33_w_randoms_min=I33_w_randoms;}

free(radial_cell);
free(radial_fkp_cell);

}while(DeltaR<40 && strcmp(type_normalization_mode, "area") == 0 && strcmp(type_normalization_mode2, "randoms") == 0 && Rsmooth==0);

if(strcmp(type_normalization_mode, "density") == 0 && Rsmooth==0){DeltaR_min=10.;}//If no area is used for normalization, keep deltaR by data
if(strcmp(type_normalization_mode, "density") == 0 && Rsmooth>0){DeltaR_min=Rsmooth;}//If no area is used for normalization, keep deltaR by data

if(strcmp(type_normalization_mode, "area") == 0 && strcmp(type_normalization_mode2, "data") == 0){DeltaR_min=parameter_value[25];}//if no randoms are explored, keep deltaR by data


//Copy needed information 
parameter_value[3]=npar_used*1.;
parameter_value[4]=Psn_2a;
parameter_value[5]=Psn_2b;
parameter_value[6]=normalizationbis;
parameter_value[7]=zeff/numeff*1.;
parameter_value[8]=num2;
parameter_value[20]=num;
parameter_value[9]=normalization;
parameter_value[10]=min;
parameter_value[11]=max;
parameter_value[12]=alpha_data;
parameter_value[28]=num3;

parameter_value[14]=I22_w_randoms_min;
parameter_value[15]=I22_w_randoms_min;
parameter_value[16]=I22_w_randoms_min;

parameter_value[17]=I33_w_randoms_min;
parameter_value[18]=I33_w_randoms_min;
parameter_value[19]=I33_w_randoms_min;
parameter_value[25]=DeltaR_min;
    
parameter_value[21]=Bsn_2a;
parameter_value[22]=Bsn_2b;

free(function_parameters);
free(L);
//exit(0);

if(free_table_dist_z==1)
{
free(redshift_table);
free(dist_table);
}

}


void get_periodic_data(char *filename_data, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],int type)
{
long int i,j,npar;
double max,min;
double wght,veto;
double params[9];
double sumweighted,Psn_1a,Bsn_1a;
FILE *f;
npar=(int)(parameter_value[3]);
f=fopen(filename_data,"r");
sumweighted=0;
Psn_1a=0;
Bsn_1a=0;
j=0;
for(i=0;i<npar;i++)
{
get_line_periodic(f, params,type);
veto=params[7];
wght=params[3];
if(veto>0 && wght!=0)
{
pos_x[j]=params[0];
pos_y[j]=params[1];
pos_z[j]=params[2];
weight[j]=params[3];
sumweighted=sumweighted+weight[j];

  Psn_1a+=pow(weight[j],2);//false pairs

  Bsn_1a+=pow(weight[j],3);//false pairs

//if(pos_z[i]-pos_z[i]!=0){printf("Error %ld %lf (%lf,%lf)\n",i,pos_z[i],pos_x[i],pos_y[i]);exit(0);}
//if(pos_z[i]/pos_z[i]!=1){printf("Error %ld %lf (%lf,%lf)\n",i,pos_z[i],pos_x[i],pos_y[i]);exit(0);}

      if(pos_x[j]>max || i==0){max=pos_x[j];}
      if(pos_y[j]>max){max=pos_y[j];}
      if(pos_z[j]>max){max=pos_z[j];}
      if(pos_x[j]<min || i==0){min=pos_x[j];}
      if(pos_y[j]<min){min=pos_y[j];}
      if(pos_z[j]<min){min=pos_z[j];}

j++;
}

}
fclose(f);
parameter_value[4]=Psn_1a;
parameter_value[21]=Bsn_1a;

parameter_value[10]=min;
parameter_value[11]=max;
parameter_value[12]=sumweighted;//alpha_data1  alpha_data1B alpha_rand1  alpha_rand1B -> alpha1 alpha1B
}

