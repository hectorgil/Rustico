#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "cubature.h"
#include "read_line.h"
#include "functions.h"
//#include "structures.h"



//void window_mask_function_RRcount(char *name_wink_out, char *name_randoms_in, int window_norm_bin,double deltaS_window, double percentage_randoms_window, char *yamamoto4window, double *parameter_value, char *header,double L,int nparallel)

void window_mask_function_RRcount(char *name_wink_out,char *name_winkBB_out,char *name_winkAB_out,char *name_randoms_in,char *name_randomsB_in, int window_norm_bin,double deltaS_window,double percentage_randoms_window,char *yamamoto4window,double *parameter_value,double *parameter_valueB,char *header,double L,int nparallel, char *type_of_code,char *shuffle, double *posX, double *posY, double *posZ, double *posW, double *posXB, double *posYB, double *posZB, double *posWB, char *Quadrupole_type, char *Hexadecapole_type, char *type_of_survey)
{
double Omega_m=parameter_value[0];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
double L2,L4,L6,L8;
double Omega_mB;
double z_minB;
double z_maxB;
    double decision_rand;
    if(strcmp(type_of_code, "rusticoX") == 0){
Omega_mB=parameter_valueB[0];
z_minB=parameter_valueB[1];
z_maxB=parameter_valueB[2];
    }
    
    
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
    
    if(strcmp(type_of_code, "rusticoX") == 0){
    if(Omega_m!=Omega_mB){printf("Warning Om different for both tracers. Exiting now...\n");exit(0);}}

double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
//double Area=parameter_value[13]*pow(Pi/180.,2);
double params[8];
double RA,dec,redshift,weight_fkp,weight_col,weight_sys,n_z,theta, radial,nuissance;//,max,min;
int veto;
long int i,j,l,npar_ran,npar_ranB,npar_used,npar_used_100,npar_usedB,npar_used_100B;
int tid;
FILE *f;
srand48(time(NULL));
double random;
double *s_x_ran,*s_y_ran,*s_z_ran,*weight_ran;
double *s_x_ranB,*s_y_ranB,*s_z_ranB,*weight_ranB;
double sumw,norm,sumwB,normB,sumwAB,normAB;
double precision=1e-6;
double **s_eff,**W0,**W2,**W4,**W6,**W8,**num_eff;
int N_bin;
double s,xlos,ylos,zlos,mu1,mu2,mu,weight_ij;
double xlos_i,xlos_j,ylos_i,ylos_j,zlos_i,zlos_j;
int index_s;

//if shuffle == 0
if(strcmp(shuffle,"no") == 0){

npar_ran=(long int)(parameter_value[3]);//total number of randoms
if(strcmp(type_of_code, "rusticoX") == 0){
npar_ranB=(long int)(parameter_valueB[3]);//total number of randoms
}

        s_x_ran = (double*) calloc(npar_ran, sizeof(double));
        s_y_ran = (double*) calloc(npar_ran, sizeof(double));
        s_z_ran = (double*) calloc(npar_ran, sizeof(double));
        weight_ran = (double*) calloc(npar_ran, sizeof(double));
    
    if(strcmp(type_of_code, "rusticoX") == 0){

        s_x_ranB = (double*) calloc(npar_ranB, sizeof(double));
        s_y_ranB = (double*) calloc(npar_ranB, sizeof(double));
        s_z_ranB = (double*) calloc(npar_ranB, sizeof(double));
        weight_ranB = (double*) calloc(npar_ranB, sizeof(double));
        
    }



npar_used=0;
npar_used_100=0;
sumw=0;
f=fopen(name_randoms_in,"r");

if( strcmp(type_of_survey,"cutsky") == 0){

for(i=0;i<npar_ran;i++)
{
get_line(f, params,1);
RA=params[0];
dec=params[1];
redshift=params[2];
weight_fkp=params[3];
//n_z=params[4];
weight_col=params[5];
weight_sys=params[6];
veto=(int)(params[7]);

        theta=90.-dec;
        if(redshift>z_min && redshift<z_max && veto==1 && weight_col>0)
        {
            random=drand48();
             npar_used_100++;
            if(percentage_randoms_window/100.>random){


                        RA=RA*Pi/180.;
                        dec=dec*Pi/180.;
                        theta=theta*Pi/180.;

                  MAX[0]=redshift;
                  MIN[0]=0;
                  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

                   s_x_ran[npar_used]=radial*sin(theta)*cos(RA);
                   s_y_ran[npar_used]=radial*sin(theta)*sin(RA);
                   s_z_ran[npar_used]=radial*cos(theta);
                   weight_ran[npar_used]=weight_fkp*weight_col*weight_sys;
                   sumw=sumw+weight_ran[npar_used];
                   npar_used++;
                        }

      }

}
}
if( strcmp(type_of_survey,"periodicFKP") == 0){
for(i=0;i<npar_ran;i++)
{
get_line_periodic(f, params,1);

            random=drand48();
            npar_used_100++;
            if(percentage_randoms_window/100.>random){

                   s_x_ran[npar_used]=params[0];
                   s_y_ran[npar_used]=params[1];
                   s_z_ran[npar_used]=params[2];
                   weight_ran[npar_used]=params[3];
                   sumw=sumw+weight_ran[npar_used];
                   npar_used++;
              }

}

}

fclose(f);
printf("%ld out of %ld randoms used for calculation\n",npar_used,npar_used_100);
    
    //------
    if(strcmp(type_of_code, "rusticoX") == 0){

    npar_usedB=0;
    npar_used_100B=0;
    sumwB=0;

    f=fopen(name_randomsB_in,"r");
    if(strcmp(type_of_survey, "cutsky") == 0){
    for(i=0;i<npar_ranB;i++)
    {
    get_line(f, params,1);
    RA=params[0];
    dec=params[1];
    redshift=params[2];
    weight_fkp=params[3];
    //n_z=params[4];
    weight_col=params[5];
    weight_sys=params[6];
    veto=(int)(params[7]);

            theta=90.-dec;
            if(redshift>z_minB && redshift<z_maxB && veto==1 && weight_col>0)
            {
                random=drand48();
                 npar_used_100B++;
                if(percentage_randoms_window/100.>random){


                            RA=RA*Pi/180.;
                            dec=dec*Pi/180.;
                            theta=theta*Pi/180.;

                      MAX[0]=redshift;
                      MIN[0]=0;
                      adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

                       s_x_ranB[npar_usedB]=radial*sin(theta)*cos(RA);
                       s_y_ranB[npar_usedB]=radial*sin(theta)*sin(RA);
                       s_z_ranB[npar_usedB]=radial*cos(theta);
                       weight_ranB[npar_usedB]=weight_fkp*weight_col*weight_sys;
                       sumwB=sumwB+weight_ranB[npar_usedB];
                       npar_usedB++;
                            }
          }
    }
}
    if(strcmp(type_of_survey, "periodicFKP") == 0){
    for(i=0;i<npar_ranB;i++)
    {
    get_line_periodic(f, params,1);

                 random=drand48();
                 npar_used_100B++;
                if(percentage_randoms_window/100.>random){


                       s_x_ranB[npar_usedB]=params[0];
                       s_y_ranB[npar_usedB]=params[1];
                       s_z_ranB[npar_usedB]=params[2];
                       weight_ranB[npar_usedB]=params[3];
                       sumwB=sumwB+weight_ranB[npar_usedB];
                       npar_usedB++;
                }

    }
}
    fclose(f);
    printf("%ld out of %ld randoms used for calculation\n",npar_usedB,npar_used_100B);
    }
    //-----

}//shuffle loop
else{
npar_used=(long int)(parameter_value[3]);
        s_x_ran = (double*) calloc(npar_used, sizeof(double));
        s_y_ran = (double*) calloc(npar_used, sizeof(double));
        s_z_ran = (double*) calloc(npar_used, sizeof(double));
        weight_ran = (double*) calloc(npar_used, sizeof(double));
        *s_x_ran=*posX;*s_y_ran=*posY;*s_z_ran=*posZ;*weight_ran=*posW;
    
    if(strcmp(type_of_code, "rusticoX") == 0){

npar_usedB=(long int)(parameter_valueB[3]);
        s_x_ranB = (double*) calloc(npar_usedB, sizeof(double));
        s_y_ranB = (double*) calloc(npar_usedB, sizeof(double));
        s_z_ranB = (double*) calloc(npar_usedB, sizeof(double));
        weight_ranB = (double*) calloc(npar_usedB, sizeof(double));
        *s_x_ranB=*posXB;*s_y_ranB=*posYB;*s_z_ranB=*posZB;*weight_ranB=*posWB;

    }



}    
    
    N_bin=(int)(L/deltaS_window);

W0=  (double**) calloc(N_bin, sizeof(double*));
W2= (double**) calloc(N_bin, sizeof(double*));
W4= (double**) calloc(N_bin, sizeof(double*));
W6= (double**) calloc(N_bin, sizeof(double*));
W8= (double**) calloc(N_bin, sizeof(double*));
s_eff= (double**) calloc(N_bin, sizeof(double*));
num_eff= (double**) calloc(N_bin, sizeof(double*));
for(i=0;i<N_bin;i++)
{
W0[i] = (double*) calloc(nparallel, sizeof(double));
W2[i] = (double*) calloc(nparallel, sizeof(double));
W4[i] = (double*) calloc(nparallel, sizeof(double));
W6[i] = (double*) calloc(nparallel, sizeof(double));
W8[i] = (double*) calloc(nparallel, sizeof(double));
s_eff[i] = (double*) calloc(nparallel, sizeof(double));
num_eff[i] = (double*) calloc(nparallel, sizeof(double));
}

printf("Starting the parallel loop...\n");
#pragma omp parallel for private(i,j,tid,s,xlos,ylos,zlos,xlos_i,xlos_j,ylos_i,ylos_j,zlos_i,zlos_j,mu,mu1,mu2,index_s,weight_ij,L2,L4,L6,L8,decision_rand) shared(s_x_ran,s_y_ran,s_z_ran,npar_used,W0,W2,W4,W6,W8,deltaS_window,N_bin,s_eff,weight_ran,num_eff,yamamoto4window,Quadrupole_type,Hexadecapole_type,type_of_survey)
for(i=0;i<npar_used;i++)
{
tid=omp_get_thread_num();//thread number
//density1=density[i];
for(j=i;j<npar_used;j++)
{
//density2=density[j];
weight_ij=weight_ran[i]*weight_ran[j];
s=sqrt( pow(s_x_ran[i]-s_x_ran[j],2)+pow(s_y_ran[i]-s_y_ran[j],2)+pow(s_z_ran[i]-s_z_ran[j],2) );

if(strcmp(type_of_survey,"periodicFKP")==0){

mu=(s_z_ran[i]-s_z_ran[j])/s;
if(i==j){mu=0;}
L2=Leg2(mu);
L4=Leg4(mu);
L6=Leg6(mu);
L8=Leg8(mu);

}
if(strcmp(yamamoto4window, "no") == 0  && strcmp(type_of_survey,"cutsky")==0){
xlos=s_x_ran[j]+(s_x_ran[i]-s_x_ran[j])/2.;
ylos=s_y_ran[j]+(s_y_ran[i]-s_y_ran[j])/2.;
zlos=s_z_ran[j]+(s_z_ran[i]-s_z_ran[j])/2.;
mu=((s_x_ran[i]-s_x_ran[j])*xlos+(s_y_ran[i]-s_y_ran[j])*ylos+(s_z_ran[i]-s_z_ran[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
if(i==j){mu=0;}
L2=Leg2(mu);
L4=Leg4(mu);
L6=Leg6(mu);
L8=Leg8(mu);
}
if(strcmp(yamamoto4window, "yes") == 0 && strcmp(type_of_survey,"cutsky")==0 ){


xlos_j=s_x_ran[j];
ylos_j=s_y_ran[j];
zlos_j=s_z_ran[j];

xlos_i=s_x_ran[i];
ylos_i=s_y_ran[i];
zlos_i=s_z_ran[i];

        decision_rand=drand48();
        if(decision_rand>=0.5){
mu1=((s_x_ran[i]-s_x_ran[j])*xlos_i+(s_y_ran[i]-s_y_ran[j])*ylos_i+(s_z_ran[i]-s_z_ran[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
mu2=((s_x_ran[i]-s_x_ran[j])*xlos_j+(s_y_ran[i]-s_y_ran[j])*ylos_j+(s_z_ran[i]-s_z_ran[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu1=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu2=0;}


}
else{
mu2=((s_x_ran[i]-s_x_ran[j])*xlos_i+(s_y_ran[i]-s_y_ran[j])*ylos_i+(s_z_ran[i]-s_z_ran[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
mu1=((s_x_ran[i]-s_x_ran[j])*xlos_j+(s_y_ran[i]-s_y_ran[j])*ylos_j+(s_z_ran[i]-s_z_ran[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu2=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu1=0;}


}

if(i==j){mu1=0;mu2=0;}

if( strcmp(Quadrupole_type,"L0L2") == 0){L2=Leg2(mu1);}//L0L2
if( strcmp(Quadrupole_type,"L1L1") == 0){L2=0.5*(3.*mu1*mu2-1.);}//L1L1

if( strcmp(Hexadecapole_type,"L0L4") == 0){L4=Leg4(mu1);}//L0L4
if( strcmp(Hexadecapole_type,"L2L2") == 0){L4=1./8.*(35.*mu1*mu1*mu2*mu2-30.*mu1*mu1+3.);}//L2L2
if( strcmp(Hexadecapole_type,"L1L3") == 0){L4=1./8.*(35.*mu1*mu1*mu1*mu2-30.*mu1*mu2+3.);}//L1L3

L6=Leg6(mu1);//L0L6 default option
L8=Leg8(mu1);//L0L8 default option

}

//mu=((s_x_ran[i]-s_x_ran[j])*xlos+(s_y_ran[i]-s_y_ran[j])*ylos+(s_z_ran[i]-s_z_ran[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
//if(i==j){mu=0;}
index_s=(int)(s/deltaS_window);
if(index_s<N_bin)
{
num_eff[index_s][tid]=num_eff[index_s][tid]+1.;
W0[index_s][tid]=W0[index_s][tid]+weight_ij/(s*s*deltaS_window);
//W2[index_s][tid]=W2[index_s][tid]+weight_ij*Leg2(mu)*5./(s*s*deltaS_window);
//W4[index_s][tid]=W4[index_s][tid]+weight_ij*Leg4(mu)*9./(s*s*deltaS_window);
//W6[index_s][tid]=W6[index_s][tid]+weight_ij*Leg6(mu)*13./(s*s*deltaS_window);
//W8[index_s][tid]=W8[index_s][tid]+weight_ij*Leg8(mu)*17./(s*s*deltaS_window);
W2[index_s][tid]=W2[index_s][tid]+weight_ij*L2*5./(s*s*deltaS_window);
W4[index_s][tid]=W4[index_s][tid]+weight_ij*L4*9./(s*s*deltaS_window);
W6[index_s][tid]=W6[index_s][tid]+weight_ij*L6*13./(s*s*deltaS_window);
W8[index_s][tid]=W8[index_s][tid]+weight_ij*L8*17./(s*s*deltaS_window);

s_eff[index_s][tid]=s_eff[index_s][tid]+s;
}

}

}

printf("End of parallel loop.... Reassigning sectors\n");

        for(i=0;i<N_bin;i++)
        {
                for(j=1;j<nparallel;j++)
                {
                        W0[i][0]=W0[i][0]+W0[i][j];
                        W2[i][0]=W2[i][0]+W2[i][j];
                        W4[i][0]=W4[i][0]+W4[i][j];
                        W6[i][0]=W6[i][0]+W6[i][j];
                        W8[i][0]=W8[i][0]+W8[i][j];
                        s_eff[i][0]=s_eff[i][0]+s_eff[i][j];
                        num_eff[i][0]=num_eff[i][0]+num_eff[i][j];

                }
       }

printf("Complete!\n");

norm=W0[window_norm_bin-1][0];


//call window
f=fopen(name_wink_out,"w");
if( strcmp(header, "yes") == 0)
{
fprintf(f,"#N= %e Sumwtot= %e\n",norm,sumw);
fprintf(f,"#seff \t sav \t W0 \t W2 \t W4 \t W6 \t W8\n");
}

for(i=1;i<N_bin;i++)
{
if(num_eff[i][0]==0){break;}//Once an empty bin is found break;

fprintf(f,"%e %e %e %e %e %e %e\n",(i+0.5)*deltaS_window,(s_eff[i][0]/num_eff[i][0]),W0[i][0]/norm,W2[i][0]/norm,W4[i][0]/norm,W6[i][0]/norm,W8[i][0]/norm);

}
fclose(f);

freeTokens(W0, N_bin);
freeTokens(W2, N_bin);
freeTokens(W4, N_bin);
freeTokens(W6, N_bin);
freeTokens(W8, N_bin);
freeTokens(s_eff,N_bin);
freeTokens(num_eff,N_bin);

//-----
    if(strcmp(type_of_code, "rusticoX") == 0){

    W0=  (double**) calloc(N_bin, sizeof(double*));
    W2= (double**) calloc(N_bin, sizeof(double*));
    W4= (double**) calloc(N_bin, sizeof(double*));
    W6= (double**) calloc(N_bin, sizeof(double*));
    W8= (double**) calloc(N_bin, sizeof(double*));
    s_eff= (double**) calloc(N_bin, sizeof(double*));
    num_eff= (double**) calloc(N_bin, sizeof(double*));
    for(i=0;i<N_bin;i++)
    {
    W0[i] = (double*) calloc(nparallel, sizeof(double));
    W2[i] = (double*) calloc(nparallel, sizeof(double));
    W4[i] = (double*) calloc(nparallel, sizeof(double));
    W6[i] = (double*) calloc(nparallel, sizeof(double));
    W8[i] = (double*) calloc(nparallel, sizeof(double));
    s_eff[i] = (double*) calloc(nparallel, sizeof(double));
    num_eff[i] = (double*) calloc(nparallel, sizeof(double));
    }

    printf("Starting the parallel loop...\n");
    #pragma omp parallel for private(i,j,tid,s,xlos,ylos,zlos,mu,mu1,mu2,xlos_i,ylos_i,zlos_i,xlos_j,ylos_j,zlos_j,index_s,weight_ij,L2,L4,L6,L8,decision_rand) shared(s_x_ranB,s_y_ranB,s_z_ranB,npar_usedB,W0,W2,W4,W6,W8,deltaS_window,N_bin,s_eff,weight_ranB,num_eff,yamamoto4window,Quadrupole_type,Hexadecapole_type,type_of_survey)
    for(i=0;i<npar_usedB;i++)
    {
    tid=omp_get_thread_num();//thread number
    //density1=density[i];
    for(j=i;j<npar_usedB;j++)
    {
    //density2=density[j];
    weight_ij=weight_ranB[i]*weight_ranB[j];
    s=sqrt( pow(s_x_ranB[i]-s_x_ranB[j],2)+pow(s_y_ranB[i]-s_y_ranB[j],2)+pow(s_z_ranB[i]-s_z_ranB[j],2) );

if(strcmp(type_of_survey,"periodicFKP")==0){

mu=(s_z_ranB[i]-s_z_ranB[j])/s;
if(i==j){mu=0;}
L2=Leg2(mu);
L4=Leg4(mu);
L6=Leg6(mu);
L8=Leg8(mu);

}

    if(strcmp(yamamoto4window, "no") == 0 && strcmp(type_of_survey,"cutsky") == 0){
    xlos=s_x_ranB[j]+(s_x_ranB[i]-s_x_ranB[j])/2.;
    ylos=s_y_ranB[j]+(s_y_ranB[i]-s_y_ranB[j])/2.;
    zlos=s_z_ranB[j]+(s_z_ranB[i]-s_z_ranB[j])/2.;
    mu=((s_x_ranB[i]-s_x_ranB[j])*xlos+(s_y_ranB[i]-s_y_ranB[j])*ylos+(s_z_ranB[i]-s_z_ranB[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
    if(i==j){mu=0;}

    L2=Leg2(mu);
    L4=Leg4(mu);
    L6=Leg6(mu);
    L8=Leg8(mu);

    }
    if(strcmp(yamamoto4window, "yes") == 0  && strcmp(type_of_survey,"cutsky") == 0){
//    xlos=s_x_ranB[j];
//    ylos=s_y_ranB[j];
//    zlos=s_z_ranB[j];

xlos_j=s_x_ranB[j];
ylos_j=s_y_ranB[j];
zlos_j=s_z_ranB[j];

xlos_i=s_x_ranB[i];
ylos_i=s_y_ranB[i];
zlos_i=s_z_ranB[i];

        decision_rand=drand48();
        if(decision_rand>=0.5){
mu1=((s_x_ranB[i]-s_x_ranB[j])*xlos_i+(s_y_ranB[i]-s_y_ranB[j])*ylos_i+(s_z_ranB[i]-s_z_ranB[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
mu2=((s_x_ranB[i]-s_x_ranB[j])*xlos_j+(s_y_ranB[i]-s_y_ranB[j])*ylos_j+(s_z_ranB[i]-s_z_ranB[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu1=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu2=0;}


}
else{
mu2=((s_x_ranB[i]-s_x_ranB[j])*xlos_i+(s_y_ranB[i]-s_y_ranB[j])*ylos_i+(s_z_ranB[i]-s_z_ranB[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
mu1=((s_x_ranB[i]-s_x_ranB[j])*xlos_j+(s_y_ranB[i]-s_y_ranB[j])*ylos_j+(s_z_ranB[i]-s_z_ranB[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu2=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu1=0;}

}
if(i==j){mu1=0;mu2=0;}

if( strcmp(Quadrupole_type,"L0L2") == 0){L2=Leg2(mu1);}//L0L2
if( strcmp(Quadrupole_type,"L1L1") == 0){L2=0.5*(3.*mu1*mu2-1.);}//L1L1

if( strcmp(Hexadecapole_type,"L0L4") == 0){L4=Leg4(mu1);}//L0L4
if( strcmp(Hexadecapole_type,"L2L2") == 0){L4=1./8.*(35.*mu1*mu1*mu2*mu2-30.*mu1*mu1+3.);}//L2L2
if( strcmp(Hexadecapole_type,"L1L3") == 0){L4=1./8.*(35.*mu1*mu1*mu1*mu2-30.*mu1*mu2+3.);}//L1L3

L6=Leg6(mu1);//L0L6 default option
L8=Leg8(mu1);//L0L8 default option



    }

//    mu=((s_x_ranB[i]-s_x_ranB[j])*xlos+(s_y_ranB[i]-s_y_ranB[j])*ylos+(s_z_ranB[i]-s_z_ranB[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
//    if(i==j){mu=0;}
    index_s=(int)(s/deltaS_window);
    if(index_s<N_bin)
    {
    num_eff[index_s][tid]=num_eff[index_s][tid]+1.;
    W0[index_s][tid]=W0[index_s][tid]+weight_ij/(s*s*deltaS_window);
  //  W2[index_s][tid]=W2[index_s][tid]+weight_ij*Leg2(mu)*5./(s*s*deltaS_window);
  //  W4[index_s][tid]=W4[index_s][tid]+weight_ij*Leg4(mu)*9./(s*s*deltaS_window);
  //  W6[index_s][tid]=W6[index_s][tid]+weight_ij*Leg6(mu)*13./(s*s*deltaS_window);
  //  W8[index_s][tid]=W8[index_s][tid]+weight_ij*Leg8(mu)*17./(s*s*deltaS_window);
    W2[index_s][tid]=W2[index_s][tid]+weight_ij*L2*5./(s*s*deltaS_window);
    W4[index_s][tid]=W4[index_s][tid]+weight_ij*L4*9./(s*s*deltaS_window);
    W6[index_s][tid]=W6[index_s][tid]+weight_ij*L6*13./(s*s*deltaS_window);
    W8[index_s][tid]=W8[index_s][tid]+weight_ij*L8*17./(s*s*deltaS_window);

    s_eff[index_s][tid]=s_eff[index_s][tid]+s;
    }

    }

    }

    printf("End of parallel loop.... Reassigning sectors\n");

            for(i=0;i<N_bin;i++)
            {
                    for(j=1;j<nparallel;j++)
                    {
                            W0[i][0]=W0[i][0]+W0[i][j];
                            W2[i][0]=W2[i][0]+W2[i][j];
                            W4[i][0]=W4[i][0]+W4[i][j];
                            W6[i][0]=W6[i][0]+W6[i][j];
                            W8[i][0]=W8[i][0]+W8[i][j];
                            s_eff[i][0]=s_eff[i][0]+s_eff[i][j];
                            num_eff[i][0]=num_eff[i][0]+num_eff[i][j];

                    }
           }

    printf("Complete!\n");

    normB=W0[window_norm_bin-1][0];


    //call window
    f=fopen(name_winkBB_out,"w");
    if( strcmp(header, "yes") == 0)
    {
    fprintf(f,"#N= %e Sumwtot= %e\n",normB,sumwB);
    fprintf(f,"#seff \t sav \t W0 \t W2 \t W4 \t W6 \t W8\n");
    }

    for(i=1;i<N_bin;i++)
    {
    if(num_eff[i][0]==0){break;}//Once an empty bin is found break;

    fprintf(f,"%e %e %e %e %e %e %e\n",(i+0.5)*deltaS_window,(s_eff[i][0]/num_eff[i][0]),W0[i][0]/normB,W2[i][0]/normB,W4[i][0]/normB,W6[i][0]/normB,W8[i][0]/normB);

    }
    fclose(f);

    freeTokens(W0, N_bin);
    freeTokens(W2, N_bin);
    freeTokens(W4, N_bin);
    freeTokens(W6, N_bin);
    freeTokens(W8, N_bin);
    freeTokens(s_eff,N_bin);
    freeTokens(num_eff,N_bin);
    
    //----
        
        W0=  (double**) calloc(N_bin, sizeof(double*));
        W2= (double**) calloc(N_bin, sizeof(double*));
        W4= (double**) calloc(N_bin, sizeof(double*));
        W6= (double**) calloc(N_bin, sizeof(double*));
        W8= (double**) calloc(N_bin, sizeof(double*));
        s_eff= (double**) calloc(N_bin, sizeof(double*));
        num_eff= (double**) calloc(N_bin, sizeof(double*));
        for(i=0;i<N_bin;i++)
        {
        W0[i] = (double*) calloc(nparallel, sizeof(double));
        W2[i] = (double*) calloc(nparallel, sizeof(double));
        W4[i] = (double*) calloc(nparallel, sizeof(double));
        W6[i] = (double*) calloc(nparallel, sizeof(double));
        W8[i] = (double*) calloc(nparallel, sizeof(double));
        s_eff[i] = (double*) calloc(nparallel, sizeof(double));
        num_eff[i] = (double*) calloc(nparallel, sizeof(double));
        }

        printf("Starting the parallel loop...\n");
        #pragma omp parallel for private(i,j,tid,s,xlos,ylos,zlos,xlos_i,ylos_i,zlos_i,xlos_j,ylos_j,zlos_j,mu,mu1,mu2,index_s,weight_ij,decision_rand,L2,L4,L6,L8) shared(s_x_ran,s_y_ran,s_z_ran,npar_used,s_x_ranB,s_y_ranB,s_z_ranB,npar_usedB,W0,W2,W4,W6,W8,deltaS_window,N_bin,s_eff,weight_ran,weight_ranB,num_eff,yamamoto4window,Quadrupole_type,Hexadecapole_type,type_of_survey)
        for(i=0;i<npar_used;i++)
        {
        tid=omp_get_thread_num();//thread number
        //density1=density[i];
        for(j=0;j<npar_usedB;j++)
        {
        //density2=density[j];
        weight_ij=weight_ran[i]*weight_ranB[j];
        s=sqrt( pow(s_x_ran[i]-s_x_ranB[j],2)+pow(s_y_ran[i]-s_y_ranB[j],2)+pow(s_z_ran[i]-s_z_ranB[j],2) );

if(strcmp(type_of_survey,"periodicFKP")==0){

mu=(s_z_ran[i]-s_z_ranB[j])/s;
if(s==0){mu=0;}
L2=Leg2(mu);
L4=Leg4(mu);
L6=Leg6(mu);
L8=Leg8(mu);

}


        if(strcmp(yamamoto4window, "no") == 0  && strcmp(type_of_survey,"cutsky") == 0){
        xlos=s_x_ranB[j]+(s_x_ran[i]-s_x_ranB[j])/2.;
        ylos=s_y_ranB[j]+(s_y_ran[i]-s_y_ranB[j])/2.;
        zlos=s_z_ranB[j]+(s_z_ran[i]-s_z_ranB[j])/2.;

        mu=((s_x_ran[i]-s_x_ranB[j])*xlos+(s_y_ran[i]-s_y_ranB[j])*ylos+(s_z_ran[i]-s_z_ranB[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
        if(sqrt(xlos*xlos+ylos*ylos+zlos*zlos)==0 || s==0){mu=0;}


        }
        if(strcmp(yamamoto4window, "yes") == 0  && strcmp(type_of_survey,"cutsky") == 0){
            
//        decision_rand=drand48();
//        if(decision_rand>=0.5){
//        xlos=s_x_ranB[j];
//        ylos=s_y_ranB[j];
//        zlos=s_z_ranB[j];}
//            else{
//        xlos=s_x_ran[i];
//        ylos=s_y_ran[i];
//        zlos=s_z_ran[i];
//            }

xlos_j=s_x_ranB[j];
ylos_j=s_y_ranB[j];
zlos_j=s_z_ranB[j];
  
xlos_i=s_x_ran[i];
ylos_i=s_y_ran[i];
zlos_i=s_z_ran[i];

        decision_rand=drand48();
        if(decision_rand>=0.5){
  mu1=((s_x_ran[i]-s_x_ranB[j])*xlos_i+(s_y_ran[i]-s_y_ranB[j])*ylos_i+(s_z_ran[i]-s_z_ranB[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
  mu2=((s_x_ran[i]-s_x_ranB[j])*xlos_j+(s_y_ran[i]-s_y_ranB[j])*ylos_j+(s_z_ran[i]-s_z_ranB[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu1=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu2=0;}

}
else{

  mu2=((s_x_ran[i]-s_x_ranB[j])*xlos_i+(s_y_ran[i]-s_y_ranB[j])*ylos_i+(s_z_ran[i]-s_z_ranB[j])*zlos_i)/(s*sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i));
  mu1=((s_x_ran[i]-s_x_ranB[j])*xlos_j+(s_y_ran[i]-s_y_ranB[j])*ylos_j+(s_z_ran[i]-s_z_ranB[j])*zlos_j)/(s*sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j));

        if(sqrt(xlos_i*xlos_i+ylos_i*ylos_i+zlos_i*zlos_i)==0 || s==0){mu2=0;}
        if(sqrt(xlos_j*xlos_j+ylos_j*ylos_j+zlos_j*zlos_j)==0 || s==0){mu1=0;}

}

if( strcmp(Quadrupole_type,"L0L2") == 0){L2=Leg2(mu1);}//L0L2
if( strcmp(Quadrupole_type,"L1L1") == 0){L2=0.5*(3.*mu1*mu2-1.);}//L1L1

if( strcmp(Hexadecapole_type,"L0L4") == 0){L4=Leg4(mu1);}//L0L4
if( strcmp(Hexadecapole_type,"L2L2") == 0){L4=1./8.*(35.*mu1*mu1*mu2*mu2-30.*mu1*mu1+3.);}//L2L2
if( strcmp(Hexadecapole_type,"L1L3") == 0){L4=1./8.*(35.*mu1*mu1*mu1*mu2-30.*mu1*mu2+3.);}//L1L3

L6=Leg6(mu1);//L0L6 default option
L8=Leg8(mu1);//L0L8 default option

          
        }

//        mu=((s_x_ran[i]-s_x_ranB[j])*xlos+(s_y_ran[i]-s_y_ranB[j])*ylos+(s_z_ran[i]-s_z_ranB[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
//        if(sqrt(xlos*xlos+ylos*ylos+zlos*zlos)==0){mu=0;}
        index_s=(int)(s/deltaS_window);
        if(index_s<N_bin)
        {
        num_eff[index_s][tid]=num_eff[index_s][tid]+1.;
        W0[index_s][tid]=W0[index_s][tid]+weight_ij/(s*s*deltaS_window);
//        W2[index_s][tid]=W2[index_s][tid]+weight_ij*Leg2(mu)*5./(s*s*deltaS_window);
//        W4[index_s][tid]=W4[index_s][tid]+weight_ij*Leg4(mu)*9./(s*s*deltaS_window);
//        W6[index_s][tid]=W6[index_s][tid]+weight_ij*Leg6(mu)*13./(s*s*deltaS_window);
//        W8[index_s][tid]=W8[index_s][tid]+weight_ij*Leg8(mu)*17./(s*s*deltaS_window);
        W2[index_s][tid]=W2[index_s][tid]+weight_ij*L2*5./(s*s*deltaS_window);
        W4[index_s][tid]=W4[index_s][tid]+weight_ij*L4*9./(s*s*deltaS_window);
        W6[index_s][tid]=W6[index_s][tid]+weight_ij*L6*13./(s*s*deltaS_window);
        W8[index_s][tid]=W8[index_s][tid]+weight_ij*L8*17./(s*s*deltaS_window);



        s_eff[index_s][tid]=s_eff[index_s][tid]+s;
        }

        }

        }

        printf("End of parallel loop.... Reassigning sectors\n");

                for(i=0;i<N_bin;i++)
                {
                        for(j=1;j<nparallel;j++)
                        {
                                W0[i][0]=W0[i][0]+W0[i][j];
                                W2[i][0]=W2[i][0]+W2[i][j];
                                W4[i][0]=W4[i][0]+W4[i][j];
                                W6[i][0]=W6[i][0]+W6[i][j];
                                W8[i][0]=W8[i][0]+W8[i][j];
                                s_eff[i][0]=s_eff[i][0]+s_eff[i][j];
                                num_eff[i][0]=num_eff[i][0]+num_eff[i][j];

                        }
               }

        printf("Complete!\n");

        normAB=W0[window_norm_bin-1][0];


        //call window
        f=fopen(name_winkAB_out,"w");
        if( strcmp(header, "yes") == 0)
        {
        fprintf(f,"#N= %e SumwtotA= %e SumwtotB= %e\n",normAB,sumw,sumwB);
        fprintf(f,"#seff \t sav \t W0 \t W2 \t W4 \t W6 \t W8\n");
        }

        for(i=1;i<N_bin;i++)
        {
        if(num_eff[i][0]==0){break;}//Once an empty bin is found break;

        fprintf(f,"%e %e %e %e %e %e %e\n",(i+0.5)*deltaS_window,(s_eff[i][0]/num_eff[i][0]),W0[i][0]/normAB,W2[i][0]/normAB,W4[i][0]/normAB,W6[i][0]/normAB,W8[i][0]/normAB);

        }
        fclose(f);

        freeTokens(W0, N_bin);
        freeTokens(W2, N_bin);
        freeTokens(W4, N_bin);
        freeTokens(W6, N_bin);
        freeTokens(W8, N_bin);
        freeTokens(s_eff,N_bin);
        freeTokens(num_eff,N_bin);
        
        
    }
    
free(s_x_ran);
free(s_y_ran);
free(s_z_ran);
free(weight_ran);
    
        if(strcmp(type_of_code, "rusticoX") == 0){
            
            free(s_x_ranB);
            free(s_y_ranB);
            free(s_z_ranB);
            free(weight_ranB);
            
        }
}
