#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include <math.h>
#include "functions.h"
#include "ps_write.h"
#include "mass_assignment.h"
//#include "structures.h"
#define indexGG(n1,n2,Ngal) (n1)*(Ngal)+(n2)
#define indexGR(ngal,nrand,Nran) (ngal)*(Nran)+(nran)

struct io_header_1//Header for the reader
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int           flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

struct particle_data
{
  float  Pos[3];
  float  Vel[3];
  float   Mass;
  int    Type;

} *P;

void load_snapshot(char *fname, int files, double *params)
{
        FILE *fd;
        char   buf[200];
        int    i,j,k,dummy,ntot_withmasses;
        int    t,n,off,pc,pc_new,pc_sph;
        int NumPart,Ngas;  
        int *Id;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

        for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
                if(files>1)
                        sprintf(buf,"%s.%d",fname,i);
                else
                        sprintf(buf,"%s",fname);

                if(!(fd=fopen(buf,"r")))
                {
                        printf("can't open file `%s`\n",buf);
                        exit(EXIT_FAILURE);
                }


                printf("Reading `%s' ...\n",buf); fflush(stdout);

                fread(&dummy, sizeof(dummy), 1, fd);
                fread(&header1, sizeof(header1), 1, fd);
                fread(&dummy, sizeof(dummy), 1, fd);

                if(files==1)
                {
                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npart[k];
                        Ngas= header1.npart[0];
                }
                else
                {
                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npartTotal[k];//Suma del numero total de particulas
                        Ngas= header1.npartTotal[0];//Suma del numero total de particulas tipo 0 (gas)
                }

                for(k=0, ntot_withmasses=0; k<5; k++)
                {
                        if(header1.mass[k]==0)
                                ntot_withmasses+= header1.npart[k];
                }

                if(i==0)
{
        printf("allocating memory...\n");

        if( P!=NULL)
        {
        free(P);
        P=NULL;
        printf("P memory freed\n");
        }

        if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
                fprintf(stderr,"failed to allocate memory.\n");
                exit(EXIT_FAILURE);
    }


        Id=malloc(NumPart*sizeof(int));
        printf("done\n");

}
         SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);//Lee las posiciones
                                pc_new++;
                        }
                }
                SKIP;

                SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);//Lee las velocidades
                                pc_new++;
                        }

                }
                SKIP;

                SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&Id[pc_new], sizeof(int), 1, fd);
                                pc_new++;
                        }
                }
                SKIP;

                fclose(fd);
    }

params[0]=NumPart;
params[1]=header1.redshift;
params[2]=header1.Omega0;
params[3]=header1.OmegaLambda;
free(Id);
}

long int count_particles_gadget(char *name_data_in)
{
        long int Ndata;
        FILE *fd;
        int dummy,k,ntot_withmasses;
         int NumPart,Ngas;

        fd=fopen(name_data_in,"r");
        if(fd==NULL){printf("Error reading %s file. Exiting now...\n",name_data_in);exit(EXIT_FAILURE);}
                fread(&dummy, sizeof(dummy), 1, fd);
                fread(&header1, sizeof(header1), 1, fd);
                fread(&dummy, sizeof(dummy), 1, fd);

                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npart[k];
                        Ngas= header1.npart[0];

                for(k=0, ntot_withmasses=0; k<5; k++)
                {
                        if(header1.mass[k]==0)
                                ntot_withmasses+= header1.npart[k];
                }


Ndata=NumPart;
return Ndata;
}

void generate_deltax_field(double *delta_k_re, double *delta_k_im, int Ngrid, int Ninterlacing, char *name_out,double Normalization)
{

  fftw_plan p;
  fftw_complex *delta_k;
  double *delta_x,d;
  long int Ntot,l,Nc2r;
  int tid;
  FILE *f;
  Nc2r=pow(Ngrid,3)/2+pow(Ngrid,2);;
  
  delta_k =  fftw_malloc(sizeof(fftw_complex)*(Nc2r));
  for(l=0;l<Nc2r;l++)
  {
    __real__ delta_k[l] = delta_k_re[l];
    __imag__ delta_k[l] = delta_k_im[l];
  }
//  free(delta_k_re);//delete them to save mem
//  free(delta_k_im);

   delta_x=(double *)calloc(Ngrid*Ngrid*Ngrid,sizeof(double));
   p =  fftw_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,delta_k,delta_x,FFTW_ESTIMATE);


   fftw_execute(p);//FFT
   fftw_destroy_plan(p);

 Ntot=(long int)(pow(Ngrid,3));
// #pragma omp parallel for private(l,tid) shared(Ngrid,Ntot,delta_x,Normalization)
f=fopen(name_out,"w");
for(l=0;l<Ntot;l++)
{
//tid=omp_get_thread_num();//thread number
d=delta_x[l]*Normalization/Ninterlacing;
fprintf(f,"%e\n",d);
}
fclose(f);
free(delta_x);

//   delta_k_re=(double *)calloc(Nc2r,sizeof(double));
//   delta_k_im=(double *)calloc(Nc2r,sizeof(double));
//  for(l=0;l<Nc2r;l++)
//  {
//    delta_k_re[l]=creal(delta_k[l]);
//    delta_k_im[l]=cimag(delta_k[l]);
//  }

fftw_free(delta_k);

}

void loop_directsum_exact_skycut(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, int ngrid, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out,char *write_kvectors, char *name_ps_kvectors)
{
int tid;
long int j,i,i1,i2;
long int xindex,yindex,zindex;
long int NGRID;//Number of k-modes between kmin and kmax
double KXX,KYY,KZZ;
double kdotx1,ckdotx1,skdotx1,musq,kampsq,kdotx2,ckdotx2,skdotx2,kam;
double SX1,SX2,SY1,SY2,SZ1,SZ2,WGHT1,WGHT2;
double SXR1,SXR2,SYR1,SYR2,SZR1,SZR2,WGHTR1,WGHTR2;
char type[200];
double Pi=(4.*atan(1.));

double **monopole_real_rr;
double **monopole_real_gr;
double **monopole_real_gg;
double **quadrupole_real_rr;
double **quadrupole_real_gr;
double **quadrupole_real_gg;
double **hexadecapole_real_rr;
double **hexadecapole_real_gr;
double **hexadecapole_real_gg;
double **dipole_real_rr;
double **dipole_real_gr;
double **dipole_real_gg;
double **octopole_real_rr;
double **octopole_real_gr;
double **octopole_real_gg;

long int ngridtot=pow(ngrid,3);
double *P0, *P2, *P4,*P1,*P3;
double *kx;

double xam12;
double mu,xam12sqrt;
double *KX;
double *KY;
double *KZ;
double SX,SY,SZ; 

kx=malloc(sizeof(double)*(ngrid));
for(i=0;i<ngrid;i++)
{
       if(i<ngrid/2+1)
       {
          kx[i]=i*1.0*(2.0*Pi/(L2-L1));
       }
       else
       {
          kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
       }
}
i=0;
for(j=0;j<ngridtot;j++)//Count the number of grid modes btween kmax and kmin
{
xindex=(long int)(j/(ngrid*ngrid*1.));
yindex=(long int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
i++;
}
}

NGRID=i;
printf("Computing aproximately %.0lf^3 (out of %d^3) number of k-modes in the range %lf < k < %lf (original range %lf < k < %lf) \n",cbrt(NGRID),ngrid,kmin,kmax, kx[0],kx[ngrid/2]);
if(NGRID>0){

KX=malloc(sizeof(double)*NGRID);
KY=malloc(sizeof(double)*NGRID);
KZ=malloc(sizeof(double)*NGRID);

i=0;
for(j=0;j<ngridtot;j++)
{
xindex=(long int)(j/(ngrid*ngrid*1.));
yindex=(long int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
KX[i]=kx[xindex];
KY[i]=kx[yindex];
KZ[i]=kx[zindex];
i++;
}
}
free(kx);

monopole_real_rr= (double **)calloc(NGRID,sizeof(double*));
monopole_real_gr= (double **)calloc(NGRID,sizeof(double*));
monopole_real_gg= (double **)calloc(NGRID,sizeof(double*));

quadrupole_real_rr= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real_gr= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real_gg= (double **)calloc(NGRID,sizeof(double*));

hexadecapole_real_rr= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_real_gr= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_real_gg= (double **)calloc(NGRID,sizeof(double*));

if(strcmp(do_odd_multipoles,"yes") == 0)
{
dipole_real_rr= (double **)calloc(NGRID,sizeof(double*));
dipole_real_gr= (double **)calloc(NGRID,sizeof(double*));
dipole_real_gg= (double **)calloc(NGRID,sizeof(double*));

octopole_real_rr= (double **)calloc(NGRID,sizeof(double*));
octopole_real_gr= (double **)calloc(NGRID,sizeof(double*));
octopole_real_gg= (double **)calloc(NGRID,sizeof(double*));

}

    for(i=0;i<NGRID;i++)
   {
    monopole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

    quadrupole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

    hexadecapole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

if(strcmp(do_odd_multipoles,"yes") == 0)
{
    dipole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    dipole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    dipole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

    octopole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    octopole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    octopole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

}

    }


#pragma omp parallel for private(SXR1,SYR1,SZR1,WGHTR1,SXR2,SYR2,SZR2,WGHTR2,SX1,SY1,SZ1,WGHT1,SX2,SY2,SZ2,WGHT2,kam,xam12sqrt,mu,j,i1,i2,tid,KXX,KYY,KZZ,kampsq,kdotx1,kdotx2,ckdotx1,ckdotx2,skdotx1,skdotx2,xam12,musq) shared(NGRID,KX,KY,KZ,Ndata,Nrand,s_x,s_y,s_z,weight,s_x_ran,s_y_ran,s_z_ran,weight_ran,do_odd_multipoles,monopole_real_gg,dipole_real_gg,quadrupole_real_gg,octopole_real_gg,hexadecapole_real_gg,monopole_real_gr,dipole_real_gr,quadrupole_real_gr,octopole_real_gr,hexadecapole_real_gr,monopole_real_rr,dipole_real_rr,quadrupole_real_rr,octopole_real_rr,hexadecapole_real_rr)
        for(j=0;j<NGRID;j++)
        {
                tid=omp_get_thread_num();//thread number
                KXX=KX[j];
                KYY=KY[j];
                KZZ=KZ[j];
                kampsq = KXX*KXX+KYY*KYY+KZZ*KZZ;
if(strcmp(do_odd_multipoles,"yes") == 0){kam=pow(kampsq,0.5);}

                for(i1=0;i1<Ndata;i1++)
                {
                     SX1=s_x[i1];SY1=s_y[i1];SZ1=s_z[i1];WGHT1=weight[i1];
                     kdotx1=KXX*SX1+KYY*SY1+KZZ*SZ1;
                     ckdotx1=cos(kdotx1)*WGHT1;
                     skdotx1=sin(kdotx1)*WGHT1;

                     for(i2=i1;i2<Ndata;i2++)
                     {
  
                       SX2=s_x[i2];SY2=s_y[i2];SZ2=s_z[i2];WGHT2=weight[i2];
                       kdotx2=KXX*SX2+KYY*SY2+KZZ*SZ2;
                       ckdotx2=cos(kdotx2)*WGHT2;
                       skdotx2=sin(kdotx2)*WGHT2;

                       xam12=0.25*(pow(SX1+SX2,2)+pow(SY1+SY2,2)+pow(SZ1+SZ2,2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);
                       

                       monopole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                       quadrupole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;                    
                          hexadecapole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;

if(strcmp(do_odd_multipoles,"yes") == 0){
                       xam12sqrt=pow(xam12,0.5);
                       mu=0.5*(kdotx1+kdotx2)/(xam12sqrt*kam);
                       dipole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu;
                       octopole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu*musq;  
}

                     }
                     for(i2=0;i2<Nrand;i2++)
                     {
                       SXR2=s_x_ran[i2];SYR2=s_y_ran[i2];SZR2=s_z_ran[i2];WGHTR2=weight_ran[i2];
                       kdotx2=KXX*SXR2+KYY*SYR2+KZZ*SZR2;
                       ckdotx2=cos(kdotx2)*WGHTR2;
                       skdotx2=sin(kdotx2)*WGHTR2;
                   
                       xam12=0.25*(pow(SX1+SXR2,2)+pow(SY1+SYR2,2)+pow(SZ1+SZR2,2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);

                       monopole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                         quadrupole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;                     
                          hexadecapole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;

if(strcmp(do_odd_multipoles,"yes") == 0){

                       xam12sqrt=pow(xam12,0.5);
                       mu=0.5*(kdotx1+kdotx2)/(xam12sqrt*kam); 
                       dipole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu;
                       octopole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu*musq;

}

                       
                }
             }
                   for(i1=0;i1<Nrand;i1++)
                   {
                     SXR1=s_x_ran[i1];SYR1=s_y_ran[i1];SZR1=s_z_ran[i1];WGHTR1=weight_ran[i1];
                     kdotx1=KXX*SXR1+KYY*SYR1+KZZ*SZR1;
                     ckdotx1=cos(kdotx1)*WGHTR1;
                     skdotx1=sin(kdotx1)*WGHTR1;
                     for(i2=i1;i2<Nrand;i2++)
                     {
                       SXR2=s_x_ran[i2];SYR2=s_y_ran[i2];SZR2=s_z_ran[i2];WGHTR2=weight_ran[i2];
                       kdotx2=KXX*SXR2+KYY*SYR2+KZZ*SZR2;
                       ckdotx2=cos(kdotx2)*WGHTR2;
                       skdotx2=sin(kdotx2)*WGHTR2;
     
                       xam12=0.25*(pow(SXR1+SXR2,2)+pow(SYR1+SYR2,2)+pow(SZR1+SZR2,2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);

                       monopole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                       quadrupole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;                       
                         hexadecapole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;

                       if(strcmp(do_odd_multipoles,"yes") == 0){
                       xam12sqrt=pow(xam12,0.5);
                       mu=0.5*(kdotx1+kdotx2)/(xam12sqrt*kam);
                       dipole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu;
                       octopole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*mu*musq;
}


                       
                      }
                   }


            }


        for(i=0;i<NGRID;i++)
        {
               for(j=1;j<n_lines_parallel;j++)
                {
                        monopole_real_gg[i][0]+=monopole_real_gg[i][j];
                        quadrupole_real_gg[i][0]+=quadrupole_real_gg[i][j];                     
                        hexadecapole_real_gg[i][0]+=hexadecapole_real_gg[i][j];
                        

                        monopole_real_gr[i][0]+=monopole_real_gr[i][j];
                        quadrupole_real_gr[i][0]+=quadrupole_real_gr[i][j];
                        hexadecapole_real_gr[i][0]+=hexadecapole_real_gr[i][j];
                        

                        monopole_real_rr[i][0]+=monopole_real_rr[i][j];
                        quadrupole_real_rr[i][0]+=quadrupole_real_rr[i][j];                        
                        hexadecapole_real_rr[i][0]+=hexadecapole_real_rr[i][j];

if(strcmp(do_odd_multipoles,"yes") == 0){

                        dipole_real_gg[i][0]+=dipole_real_gg[i][j];
                        octopole_real_gg[i][0]+=octopole_real_gg[i][j];

                        dipole_real_gr[i][0]+=dipole_real_gr[i][j];
                        octopole_real_gr[i][0]+=octopole_real_gr[i][j];


                        dipole_real_rr[i][0]+=dipole_real_rr[i][j];
                        octopole_real_rr[i][0]+=octopole_real_rr[i][j];
}
                        

                }

                       monopole_real_gg[i][0]+=(alpha*alpha*monopole_real_rr[i][0]-alpha*monopole_real_gr[i][0]);
                       quadrupole_real_gg[i][0]+=(alpha*alpha*quadrupole_real_rr[i][0]-alpha*quadrupole_real_gr[i][0]);                       
                       hexadecapole_real_gg[i][0]+=(alpha*alpha*hexadecapole_real_rr[i][0]-alpha*hexadecapole_real_gr[i][0]);
 if(strcmp(do_odd_multipoles,"yes") == 0){
                       dipole_real_gg[i][0]+=(alpha*alpha*dipole_real_rr[i][0]-alpha*dipole_real_gr[i][0]);
                       octopole_real_gg[i][0]+=(alpha*alpha*octopole_real_rr[i][0]-alpha*octopole_real_gr[i][0]);
}                      


        }

freeTokens(monopole_real_rr,NGRID);
freeTokens(monopole_real_gr,NGRID);

freeTokens(quadrupole_real_rr,NGRID);
freeTokens(quadrupole_real_gr,NGRID);

freeTokens(hexadecapole_real_rr,NGRID);
freeTokens(hexadecapole_real_gr,NGRID);

 if(strcmp(do_odd_multipoles,"yes") == 0){
freeTokens(dipole_real_rr,NGRID);
freeTokens(dipole_real_gr,NGRID);

freeTokens(octopole_real_rr,NGRID);
freeTokens(octopole_real_gr,NGRID);
}

P0= (double *)calloc(NGRID,sizeof(double));
P2= (double *)calloc(NGRID,sizeof(double));
P4= (double *)calloc(NGRID,sizeof(double));
 if(strcmp(do_odd_multipoles,"yes") == 0){
P1= (double *)calloc(NGRID,sizeof(double));
P3= (double *)calloc(NGRID,sizeof(double));
}

#pragma omp parallel for private(i) shared(NGRID,P0,P1,P2,P3,P4,monopole_real_gg, quadrupole_real_gg,hexadecapole_real_gg,dipole_real_gg, octopole_real_gg)
        for(i=0;i<NGRID;i++)
        {
             P0[i]=monopole_real_gg[i][0];
             P2[i]=quadrupole_real_gg[i][0];             
             P4[i]=hexadecapole_real_gg[i][0];
          if(strcmp(do_odd_multipoles,"yes") == 0){
             P1[i]=quadrupole_real_gg[i][0];
             P3[i]=hexadecapole_real_gg[i][0];

}    

        }

freeTokens(monopole_real_gg,NGRID);
freeTokens(quadrupole_real_gg,NGRID);
freeTokens(hexadecapole_real_gg,NGRID);

          if(strcmp(do_odd_multipoles,"yes") == 0){
freeTokens(dipole_real_gg,NGRID);
freeTokens(octopole_real_gg,NGRID);
}

printf("Writing Power Spectrum output %s...",name_ps_out);

sprintf(type,"DSE");

    write_power_spectrum_skyscuts(kmin,kmax,KX, KY,KZ, P0, NULL,P1, NULL, P2, NULL, P3, NULL, P4, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, Deltak, ngrid, NGRID, L1, L2, I22,0, name_ps_out,NULL,NULL, P_shot_noise,0, binning_type, do_anisotropy, do_odd_multipoles, Quadrupole_type, Octopole_type, Hexadecapole_type, type,1,NULL, write_kvectors, name_ps_kvectors);


printf("Ok!\n");

free(P0);
free(P2);
free(P4);
if(strcmp(do_odd_multipoles,"yes") == 0){
free(P1);
free(P3);}

free(KX);
free(KY);
free(KZ);

}//NGRID=0
else{
free(kx);

}



}




void loop_directsum_yamamoto_skycut(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double  L1, double L2, int ngrid, double P_shot_noise, double  Deltak, double  I22, double  alpha, int  n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out,char *write_kvectors, char *name_ps_kvectors)
{
int tid;
long int j,i;
long int xindex,yindex,zindex;
long int NGRID;//Number of k-modes between kmin and kmax
double KXX,KYY,KZZ,XAM,kamp,mu,XAM05;
double SXR1,SYR1,SZR1,SXR2,SYR2,SZR2,WGHTR1;
double SX1,SY1,SZ1,SX2,SY2,SZ2,WGHT1;
double kdotx,ckdotx,skdotx,musq,kampsq;
double Pi=(4.*atan(1.));
char type[200];

double **monopole_real;
double **monopole_imag;
double **quadrupole_real;
double **quadrupole_imag;
double **hexadecapole_real;
double **hexadecapole_imag;

double **dipole_real;
double **dipole_imag;
double **octopole_real;
double **octopole_imag;


double **monopole_real_rand;
double **monopole_imag_rand;
double **quadrupole_real_rand;
double **quadrupole_imag_rand;
double **hexadecapole_real_rand;
double **hexadecapole_imag_rand;

double **dipole_real_rand;
double **dipole_imag_rand;
double **octopole_real_rand;
double **octopole_imag_rand;

long int ngridtot=pow(ngrid,3);

double *deltak_re0, *deltak_im0,*deltak_re1, *deltak_im1, *deltak_re2, *deltak_im2, *deltak_re3, *deltak_im3,*deltak_re4, *deltak_im4;
double *kx;
double *xampsq;
double *xampsq_ran;
double *xamp;
double *xamp_ran;
double *KX;
double *KY;
double *KZ;
int delta1_sw,delta2_sw,delta3_sw,delta4_sw;
double SX,SY,SZ,VEC;
delta1_sw=0;
delta2_sw=0;
delta3_sw=0;
delta4_sw=0;

if(strcmp(Hexadecapole_type, "L2L2") == 0){delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L0L4") == 0){delta4_sw=1;delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L1L3") == 0){delta1_sw=1;delta3_sw=1;}

if(strcmp(Quadrupole_type, "L0L2") == 0){delta2_sw=1;}
if(strcmp(Quadrupole_type, "L1L1") == 0){delta1_sw=1;}

if(strcmp(do_odd_multipoles,"yes") == 0){
if(strcmp(Octopole_type, "L0L3") == 0){delta3_sw=1;delta1_sw=1;}
if(strcmp(Octopole_type, "L1L2") == 0){delta1_sw=1;delta2_sw=1;}
}

if(strcmp(do_anisotropy,"no") == 0){delta1_sw=0;delta2_sw=0;delta3_sw=0;delta4_sw=0;}

if(strcmp(do_anisotropy,"yes") == 0)
{

if(delta2_sw == 1 || delta4_sw == 1){
        xampsq=(double*)malloc(sizeof(double)*Ndata);
        xampsq_ran=(double*)malloc(sizeof(double)*Nrand);}

if(delta1_sw == 1 || delta3_sw == 1){
        xamp=(double*)malloc(sizeof(double)*Ndata);
        xamp_ran=(double*)malloc(sizeof(double)*Nrand);}


for(i=0;i<Ndata;i++){
SX=s_x[i]*s_x[i];
SY=s_y[i]*s_y[i];
SZ=s_z[i]*s_z[i];
if(delta2_sw == 1 || delta4_sw == 1){xampsq[i]=SX+SY+SZ;}
if(delta1_sw == 1 || delta3_sw == 1){xamp[i]=sqrt(SX+SY+SZ);}
}
for(i=0;i<Nrand;i++){
SX=s_x_ran[i]*s_x_ran[i];
SY=s_y_ran[i]*s_y_ran[i];
SZ=s_z_ran[i]*s_z_ran[i];
if(delta2_sw == 1 || delta4_sw == 1){xampsq_ran[i]=SX+SY+SZ;}
if(delta1_sw == 1 || delta3_sw == 1){xamp_ran[i]=sqrt(SX+SY+SZ);}
}

}

kx=malloc(sizeof(double)*(ngrid));
for(i=0;i<ngrid;i++)
{
       if(i<ngrid/2+1)
       {
          kx[i]=i*1.0*(2.0*Pi/(L2-L1));
       }
       else
       {
          kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
       }
}
i=0;
for(j=0;j<ngridtot;j++)//Count the number of grid modes btween kmax and kmin
{
xindex=(long int)(j/(ngrid*ngrid*1.));
yindex=(long int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
VEC=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(VEC<=kmax*kmax && VEC>kmin*kmin && SZ>=0){//open-closed binning definition here: (kmin< k <= kmax )
//if(VEC<kmax*kmax && VEC>=kmin*kmin && SZ>=0){//open-closed binning definition here: (kmin <= k < kmax )
i++;//printf("Ngrid:%d, %lf, %lf\n",i,SX,SZ);
}
}

NGRID=i;
printf("Computing aproximately %.0lf^3 = %ld (out of %d^3) number of k-modes in the range %lf < k < %lf (original range %lf < k < %lf) \n",cbrt(NGRID),NGRID,ngrid,kmin,kmax, kx[0],kx[ngrid/2]);

if(NGRID>0){
KX=malloc(sizeof(double)*NGRID);
KY=malloc(sizeof(double)*NGRID);
KZ=malloc(sizeof(double)*NGRID);


i=0;
for(j=0;j<ngridtot;j++)
{
xindex=(long int)(j/(ngrid*ngrid*1.));
yindex=(long int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
VEC=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(VEC<=kmax*kmax && VEC>kmin*kmin && SZ>=0){//open-closed binning definition here: (kmin< k <= kmax )
//if(VEC<kmax*kmax && VEC>=kmin*kmin && SZ>=0){//open-closed binning definition here: (kmin<= k < kmax )
KX[i]=kx[xindex];
KY[i]=kx[yindex];
KZ[i]=kx[zindex];//printf("%lf, %lf, %lf\n",KX[i],KY[i],KZ[i]);
i++;//printf("i:%d, %lf, %lf\n",i,SX,SZ);
}
}
free(kx);
//if(i>0){exit(0);}
if(NGRID != i){printf("Error selecting k-vectors: %ld, %ld. Check conditions in both loops above this point. Exiting now...\n",NGRID,i);exit(0);}

monopole_real= (double **)calloc(NGRID,sizeof(double*));
monopole_imag= (double **)calloc(NGRID,sizeof(double*));

if(strcmp(do_anisotropy,"yes") == 0)
{

if(delta2_sw == 1){
quadrupole_real= (double **)calloc(NGRID,sizeof(double*));
quadrupole_imag= (double **)calloc(NGRID,sizeof(double*));}

if(delta4_sw == 1){
hexadecapole_real= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_imag= (double **)calloc(NGRID,sizeof(double*));}

if(delta1_sw == 1){
dipole_real= (double **)calloc(NGRID,sizeof(double*));
dipole_imag= (double **)calloc(NGRID,sizeof(double*));}

if(delta3_sw == 1){
octopole_real= (double **)calloc(NGRID,sizeof(double*));
octopole_imag= (double **)calloc(NGRID,sizeof(double*));}

}

monopole_real_rand= (double **)calloc(NGRID,sizeof(double*));
monopole_imag_rand= (double **)calloc(NGRID,sizeof(double*));

if(strcmp(do_anisotropy,"yes") == 0)
{

if(delta2_sw == 1){
quadrupole_real_rand= (double **)calloc(NGRID,sizeof(double*));
quadrupole_imag_rand= (double **)calloc(NGRID,sizeof(double*));}

if(delta4_sw == 1){
hexadecapole_real_rand= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_imag_rand= (double **)calloc(NGRID,sizeof(double*));}

if(delta1_sw == 1){
dipole_real_rand= (double **)calloc(NGRID,sizeof(double*));
dipole_imag_rand= (double **)calloc(NGRID,sizeof(double*));}

if(delta3_sw == 1){
octopole_real_rand= (double **)calloc(NGRID,sizeof(double*));
octopole_imag_rand= (double **)calloc(NGRID,sizeof(double*));}

}
    for(i=0;i<NGRID;i++)
   {
    monopole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));
 
if(strcmp(do_anisotropy,"yes") == 0)
{

if(delta2_sw == 1){   
    quadrupole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta4_sw == 1){
    hexadecapole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta1_sw == 1){
    dipole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    dipole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta3_sw == 1){
    octopole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    octopole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));}
}

    monopole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));

if(strcmp(do_anisotropy,"yes") == 0)
{

if(delta2_sw == 1){
    quadrupole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta4_sw == 1){
    hexadecapole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta1_sw == 1){
    dipole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    dipole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));}

if(delta3_sw == 1){
    octopole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    octopole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));}
}

  }

        #pragma omp parallel for private(SXR1,SYR1,SZR1,WGHTR1,SX1,SY1,SZ1,WGHT1,j,tid,KXX,KYY,KZZ,kampsq,kamp,i,XAM,XAM05,kdotx,musq,mu,ckdotx,skdotx) shared(NGRID,KX,KY,KZ,Ndata,Nrand,xampsq,xamp,s_x,s_y,s_z,monopole_real,monopole_imag,quadrupole_real,quadrupole_imag,hexadecapole_real,hexadecapole_imag,dipole_real,dipole_imag,octopole_real,octopole_imag,weight,xampsq_ran,xamp_ran,s_x_ran,s_y_ran,s_z_ran,weight_ran, monopole_real_rand,monopole_imag_rand, dipole_real_rand,dipole_imag_rand,quadrupole_real_rand,quadrupole_imag_rand,octopole_real_rand,octopole_imag_rand,hexadecapole_real_rand,hexadecapole_imag_rand,do_anisotropy,do_odd_multipoles, delta1_sw,delta2_sw,delta3_sw,delta4_sw)
        for(j=0;j<NGRID;j++)
        {
                tid=omp_get_thread_num();//thread number
                KXX=KX[j];
                KYY=KY[j];
                KZZ=KZ[j];
if(strcmp(do_anisotropy,"yes") == 0)
{
if(delta2_sw == 1 || delta4_sw ==1){kampsq = KXX*KXX+KYY*KYY+KZZ*KZZ;}
if(delta1_sw == 1 || delta3_sw ==1){kamp = sqrt(KXX*KXX+KYY*KYY+KZZ*KZZ);}
}

        for(i=0;i<Ndata;i++)
        {
                        SX1=s_x[i];SY1=s_y[i];SZ1=s_z[i];WGHT1=weight[i];
                        kdotx=KXX*SX1+KYY*SY1+KZZ*SZ1;
                        ckdotx=cos(kdotx)*WGHT1;
                        skdotx=sin(kdotx)*WGHT1;
                         
                        if(strcmp(do_anisotropy,"yes") == 0){
                        if(delta2_sw == 1 || delta4_sw == 1 )
                        {
                        XAM=xampsq[i];                        
                        musq = kdotx*kdotx/(kampsq*XAM);}

                        if(delta1_sw == 1 || delta3_sw == 1){
                        XAM05=xamp[i];//////////////////////////
                        mu = kdotx/(kamp*XAM05);/////////////////////
                              }
                        }

                       monopole_real[j][tid]+=ckdotx;
                       monopole_imag[j][tid]+=skdotx;

                       if(strcmp(do_anisotropy,"yes") == 0)
                       {
                       if(delta2_sw == 1){
                       quadrupole_real[j][tid]+=ckdotx*musq;
                       quadrupole_imag[j][tid]+=skdotx*musq;}
                       
                       if(delta4_sw == 1){
                       hexadecapole_real[j][tid]+=ckdotx*musq*musq;
                       hexadecapole_imag[j][tid]+=skdotx*musq*musq;}

                       if(delta1_sw == 1){
                       dipole_real[j][tid]+=ckdotx*mu;
                       dipole_imag[j][tid]+=skdotx*mu;}

                       if(delta3_sw == 1){
                       octopole_real[j][tid]+=ckdotx*mu*musq;
                       octopole_imag[j][tid]+=skdotx*mu*musq;}


                       }
                        SXR1=s_x_ran[i];SYR1=s_y_ran[i];SZR1=s_z_ran[i];WGHTR1=weight_ran[i];
                        kdotx=KXX*SXR1+KYY*SYR1+KZZ*SZR1;

                      if(strcmp(do_anisotropy,"yes") == 0){

                        if(delta2_sw ==1 || delta4_sw ==1){
                        XAM=xampsq_ran[i];
                        musq = kdotx*kdotx/(kampsq*XAM);}

                       if(delta1_sw ==1 || delta3_sw ==1){
                        XAM05=xamp_ran[i];//////////////////////////
                        mu = kdotx/(kamp*XAM05);}/////////////////////

                                          }      
                      
                        ckdotx=cos(kdotx)*WGHTR1;
                        skdotx=sin(kdotx)*WGHTR1;
                      
                        monopole_real_rand[j][tid]+=ckdotx;
                        monopole_imag_rand[j][tid]+=skdotx;

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
                        quadrupole_real_rand[j][tid]+=ckdotx*musq;
                        quadrupole_imag_rand[j][tid]+=skdotx*musq;}

                        if(delta4_sw == 1){
                        hexadecapole_real_rand[j][tid]+=ckdotx*musq*musq;
                        hexadecapole_imag_rand[j][tid]+=skdotx*musq*musq;}

                       if(delta1_sw == 1){
                        dipole_real_rand[j][tid]+=ckdotx*mu;
                        dipole_imag_rand[j][tid]+=skdotx*mu;}

                       if(delta3_sw == 1){
                        octopole_real_rand[j][tid]+=ckdotx*mu*musq;
                        octopole_imag_rand[j][tid]+=skdotx*mu*musq;}

}

         }
         for(i=Ndata;i<Nrand;i++)
         {

                        SXR1=s_x_ran[i];SYR1=s_y_ran[i];SZR1=s_z_ran[i];WGHTR1=weight_ran[i];
                        kdotx=KXX*SXR1+KYY*SYR1+KZZ*SZR1;

     if(strcmp(do_anisotropy,"yes") == 0){

     if(delta2_sw ==1 || delta4_sw ==1){
                        XAM=xampsq_ran[i];
                        musq = kdotx*kdotx/(kampsq*XAM);}

     if(delta1_sw ==1 || delta3_sw ==1){
                        XAM05=xamp_ran[i];//////////////////////////                        
                        mu = kdotx/(kamp*XAM05);}/////////////////////
     }                        
                        ckdotx=cos(kdotx)*WGHTR1;
                        skdotx=sin(kdotx)*WGHTR1;

                        monopole_real_rand[j][tid]+=ckdotx;
                        monopole_imag_rand[j][tid]+=skdotx;

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
                        quadrupole_real_rand[j][tid]+=ckdotx*musq;
                        quadrupole_imag_rand[j][tid]+=skdotx*musq;}

                       if(delta4_sw == 1){
                        hexadecapole_real_rand[j][tid]+=ckdotx*musq*musq;
                        hexadecapole_imag_rand[j][tid]+=skdotx*musq*musq;}

                       if(delta1_sw == 1){
                        dipole_real_rand[j][tid]+=ckdotx*mu;
                        dipole_imag_rand[j][tid]+=skdotx*mu;}

                       if(delta3_sw == 1){
                        octopole_real_rand[j][tid]+=ckdotx*mu*musq;
                        octopole_imag_rand[j][tid]+=skdotx*mu*musq;}

}


         }

        }


     if(strcmp(do_anisotropy,"yes") == 0){

if(delta2_sw == 1 || delta4_sw == 1){
free(xampsq);
free(xampsq_ran);}

if(delta1_sw == 1 || delta3_sw == 1){
free(xamp);
free(xamp_ran);}

}
alpha=-alpha;
 for(i=0;i<NGRID;i++)
        {
                for(j=1;j<n_lines_parallel;j++)
                {

                        monopole_real[i][0]+=monopole_real[i][j];
                        monopole_imag[i][0]+=monopole_imag[i][j];

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
                        quadrupole_real[i][0]+=quadrupole_real[i][j];
                        quadrupole_imag[i][0]+=quadrupole_imag[i][j];}
                       
                       if(delta4_sw == 1){
                        hexadecapole_real[i][0]+=hexadecapole_real[i][j];
                        hexadecapole_imag[i][0]+=hexadecapole_imag[i][j];}

                       if(delta1_sw == 1){
                        dipole_real[i][0]+=dipole_real[i][j];
                        dipole_imag[i][0]+=dipole_imag[i][j];}

                       if(delta3_sw == 1){
                        octopole_real[i][0]+=octopole_real[i][j];
                        octopole_imag[i][0]+=octopole_imag[i][j];}
}

                        monopole_real_rand[i][0]+=monopole_real_rand[i][j];
                        monopole_imag_rand[i][0]+=monopole_imag_rand[i][j];

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
                        quadrupole_real_rand[i][0]+=quadrupole_real_rand[i][j];
                        quadrupole_imag_rand[i][0]+=quadrupole_imag_rand[i][j];}
                       
                       if(delta4_sw == 1){
                        hexadecapole_real_rand[i][0]+=hexadecapole_real_rand[i][j];
                        hexadecapole_imag_rand[i][0]+=hexadecapole_imag_rand[i][j];}

                       if(delta1_sw == 1){
                        dipole_real_rand[i][0]+=dipole_real_rand[i][j];
                        dipole_imag_rand[i][0]+=dipole_imag_rand[i][j];}

                       if(delta3_sw == 1){
                        octopole_real_rand[i][0]+=octopole_real_rand[i][j];
                        octopole_imag_rand[i][0]+=octopole_imag_rand[i][j];}
}



                }
                       monopole_real[i][0]+=alpha*monopole_real_rand[i][0];
                       monopole_imag[i][0]+=alpha*monopole_imag_rand[i][0];

     if(strcmp(do_anisotropy,"yes") == 0){
                       if(delta2_sw == 1){
                       quadrupole_real[i][0]+=alpha*quadrupole_real_rand[i][0];
                       quadrupole_imag[i][0]+=alpha*quadrupole_imag_rand[i][0];}

                       if(delta4_sw == 1){
                       hexadecapole_real[i][0]+=alpha*hexadecapole_real_rand[i][0];
                       hexadecapole_imag[i][0]+=alpha*hexadecapole_imag_rand[i][0];}

                       if(delta1_sw == 1){
                       dipole_real[i][0]+=alpha*dipole_real_rand[i][0];
                       dipole_imag[i][0]+=alpha*dipole_imag_rand[i][0];}

                       if(delta3_sw == 1){
                       octopole_real[i][0]+=alpha*octopole_real_rand[i][0];
                       octopole_imag[i][0]+=alpha*octopole_imag_rand[i][0];}
}


        }
freeTokens(monopole_real_rand,NGRID);
freeTokens(monopole_imag_rand,NGRID);

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
freeTokens(quadrupole_real_rand,NGRID);
freeTokens(quadrupole_imag_rand,NGRID);}

                       if(delta4_sw == 1){
freeTokens(hexadecapole_real_rand,NGRID);
freeTokens(hexadecapole_imag_rand,NGRID);}

                       if(delta1_sw == 1){
freeTokens(dipole_real_rand,NGRID);
freeTokens(dipole_imag_rand,NGRID);}

                       if(delta3_sw == 1){
freeTokens(octopole_real_rand,NGRID);
freeTokens(octopole_imag_rand,NGRID);}
}

deltak_re0= (double *)calloc(NGRID,sizeof(double));
deltak_im0= (double *)calloc(NGRID,sizeof(double));

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
deltak_re2= (double *)calloc(NGRID,sizeof(double));
deltak_im2= (double *)calloc(NGRID,sizeof(double));}

                       if(delta4_sw == 1){
deltak_re4= (double *)calloc(NGRID,sizeof(double));
deltak_im4= (double *)calloc(NGRID,sizeof(double));}

                       if(delta1_sw == 1){
deltak_re1= (double *)calloc(NGRID,sizeof(double));
deltak_im1= (double *)calloc(NGRID,sizeof(double));}

                       if(delta3_sw == 1){
deltak_re3= (double *)calloc(NGRID,sizeof(double));
deltak_im3= (double *)calloc(NGRID,sizeof(double));}

}

#pragma omp parallel for private(i) shared(NGRID,deltak_re0,deltak_re2,deltak_re4,deltak_re1,deltak_re3, deltak_im0,deltak_im2,deltak_im4,deltak_im1,deltak_im3,monopole_real,monopole_imag,quadrupole_real, quadrupole_imag,hexadecapole_real,hexadecapole_imag,dipole_real, dipole_imag,octopole_real,octopole_imag,delta1_sw,delta2_sw,delta3_sw,delta4_sw,do_anisotropy)
 for(i=0;i<NGRID;i++)
        {
             deltak_re0[i]=monopole_real[i][0];
             deltak_im0[i]=monopole_imag[i][0];
         
     if(strcmp(do_anisotropy,"yes") == 0){    
                       if(delta2_sw == 1){
             deltak_re2[i]=quadrupole_real[i][0];
             deltak_im2[i]=quadrupole_imag[i][0];}
                                    if(delta4_sw == 1){
             deltak_re4[i]=hexadecapole_real[i][0];
             deltak_im4[i]=hexadecapole_imag[i][0];}
                       if(delta1_sw == 1){
             deltak_re1[i]=dipole_real[i][0];
             deltak_im1[i]=dipole_imag[i][0];}
                       if(delta3_sw == 1){
             deltak_re3[i]=octopole_real[i][0];
             deltak_im3[i]=octopole_imag[i][0];}
             }

        }
freeTokens(monopole_real,NGRID);
freeTokens(monopole_imag,NGRID);

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
freeTokens(quadrupole_real,NGRID);
freeTokens(quadrupole_imag,NGRID);}

                       if(delta4_sw == 1){
freeTokens(hexadecapole_real,NGRID);
freeTokens(hexadecapole_imag,NGRID);}

                       if(delta1_sw == 1){
freeTokens(dipole_real,NGRID);
freeTokens(dipole_imag,NGRID);}

                       if(delta3_sw == 1){
freeTokens(octopole_real,NGRID);
freeTokens(octopole_imag,NGRID);}

}

printf("Writing Power Spectrum output,  %s...",name_ps_out);

sprintf(type,"DSY");
    write_power_spectrum_skyscuts(kmin,kmax,KX, KY,KZ, deltak_re0, deltak_im0,deltak_re1, deltak_im1, deltak_re2, deltak_im2, deltak_re3, deltak_im3, deltak_re4, deltak_im4,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, Deltak, ngrid, NGRID, L1, L2, I22,0, name_ps_out,NULL,NULL, P_shot_noise,0, binning_type, do_anisotropy, do_odd_multipoles, Quadrupole_type, Octopole_type, Hexadecapole_type, type,1,NULL,write_kvectors, name_ps_kvectors);


free(deltak_re0);
free(deltak_im0);

     if(strcmp(do_anisotropy,"yes") == 0){

                       if(delta2_sw == 1){
free(deltak_re2);
free(deltak_im2);}

                       if(delta4_sw == 1){
free(deltak_re4);
free(deltak_im4);}

                       if(delta1_sw == 1){
free(deltak_re1);
free(deltak_im1);}

                       if(delta3_sw == 1){
free(deltak_re3);
free(deltak_im3);}

}


printf("Ok!\n");

free(KX);
free(KY);
free(KZ);


}//NGRID=0
else{
free(kx);

     if(strcmp(do_anisotropy,"yes") == 0){

if(delta2_sw == 1 || delta4_sw == 1){
free(xampsq);
free(xampsq_ran);}

if(delta1_sw == 1 || delta3_sw == 1){
free(xamp);
free(xamp_ran);}

}


}



}


void loop_directsum_yamamoto_skycut_caller(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type,char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out, char *type_of_computation,char *write_kvectors, char *name_ps_kvectors)
{

double *kmin_i;
double *kmax_i;
int i;
int Nbins;
int ngrid;
double Pi=(4.*atan(1.));

if(strcmp(binning_type, "linear") == 0){Nbins=(int)((kmax-kmin)/Deltak)+1;}
if(strcmp(binning_type, "log10") == 0){Nbins=(int)((log10(kmax)-log10(kmin))/Deltak)+1;}

kmin_i= (double*)calloc(Nbins,sizeof(double));
kmax_i= (double*)calloc(Nbins,sizeof(double));

for(i=0;i<Nbins;i++)
{
//linear
if(strcmp(binning_type, "linear") == 0)
{
kmin_i[i]=kmin+Deltak*i;
kmax_i[i]=kmin+Deltak*(i+1);
}
//log
if(strcmp(binning_type, "log10") == 0)
{
kmin_i[i]=pow(10,log10(kmin)+Deltak*i);
kmax_i[i]=pow(10,log10(kmin)+Deltak*(i+1));
}

//Determine the appropiate value of ngrid
ngrid=1;
do
{

ngrid=ngrid*2;

}while(kmax_i[i]>Pi*ngrid/(L2-L1));

if(strcmp(type_of_computation, "DSY") == 0)
{
loop_directsum_yamamoto_skycut(kmin_i[i],kmax_i[i], s_x, s_y, s_z, weight, Ndata, s_x_ran, s_y_ran, s_z_ran, weight_ran, Nrand, L1, L2, ngrid, P_shot_noise, Deltak, I22, alpha, n_lines_parallel, binning_type, Quadrupole_type, Octopole_type, Hexadecapole_type, do_odd_multipoles, do_anisotropy, name_ps_out,write_kvectors, name_ps_kvectors);
}
if(strcmp(type_of_computation, "DSE") == 0)
{
loop_directsum_exact_skycut(kmin_i[i],kmax_i[i], s_x, s_y, s_z, weight, Ndata, s_x_ran, s_y_ran, s_z_ran, weight_ran, Nrand, L1, L2, ngrid, P_shot_noise, Deltak, I22, alpha, n_lines_parallel, binning_type,  Quadrupole_type, Octopole_type, Hexadecapole_type, do_odd_multipoles, do_anisotropy, name_ps_out, write_kvectors, name_ps_kvectors);
}

}

free(s_x);
free(s_y);
free(s_z);
free(weight);
free(s_x_ran);
free(s_y_ran);
free(s_z_ran);
free(weight_ran);
free(kmin_i);
free(kmax_i);

}

void modify_input_for_yamamoto(int mode_yamamoto, int mode_yamamoto_prev, double in[], int ngrid, double L1, double L2)
{
long int i,j,k;
long int c;
double kmod2;
double prev_factor,curr_factor;
double r_x,r_y,r_z;
long int ngridtot=pow(ngrid,3);
        
        #pragma omp parallel for private(c,i,j,k,kmod2,prev_factor,curr_factor,r_x,r_y,r_z) shared(ngrid,ngridtot,L1,L2,in,mode_yamamoto,mode_yamamoto_prev)
		for(c=0;c<ngridtot;c++)
		{
			i=(long int)(c/(ngrid*ngrid*1.));
			j=(long int)( (c-i*ngrid*ngrid)/(ngrid*1.));
			k=c-i*ngrid*ngrid-j*ngrid;
                        
            kmod2=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)+(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)+(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.);//kmod2 should be always different than 0
            r_x=(L1+i*(L2-L1)/ngrid*1.);
            r_y=(L1+j*(L2-L1)/ngrid*1.);
            r_z=(L1+k*(L2-L1)/ngrid*1.);

        if(mode_yamamoto==0)
        {   curr_factor=1;//in[c]*=r_x*r_x/(kmod2);     
        }
        if(mode_yamamoto_prev==0)
        {   prev_factor=1;//in[c]*=r_x*r_x/(kmod2);     
        }
        if(mode_yamamoto==1)
        {   curr_factor=r_x*r_x/kmod2;//in[c]*=r_x*r_x/(kmod2);	    
        }
        if(mode_yamamoto_prev==1)
        {   prev_factor=r_x*r_x/kmod2;//in[c]*=r_x*r_x/(kmod2);     
        }
        if(mode_yamamoto==2)
        {   curr_factor=r_x*r_y/kmod2;//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/(L1+i*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto_prev==2)
        {   prev_factor=r_x*r_y/kmod2;//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/(L1+i*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto==3)
        {   curr_factor=r_x*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/(L1+j*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto_prev==3)
        {   prev_factor=r_x*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/(L1+j*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto==4)
        {   curr_factor=r_y*r_y/kmod2;//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==4)
        {   prev_factor=r_y*r_y/kmod2;//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==5)
        {   curr_factor=r_y*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==5)
        {   prev_factor=r_y*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==6)
        {   curr_factor=r_z*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==6)
        {   prev_factor=r_z*r_z/kmod2;//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==7)
        {   curr_factor=r_x*r_x*r_x*r_x/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/(kmod2*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==7)
        {   prev_factor=r_x*r_x*r_x*r_x/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/(kmod2*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==8)
        {   curr_factor=r_y*r_y*r_y*r_y/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==8)
        {   prev_factor=r_y*r_y*r_y*r_y/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==9)
        {   curr_factor=r_z*r_z*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==9)
        {   prev_factor=r_z*r_z*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==10)
        {   curr_factor=r_x*r_x*r_x*r_y/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==10)
        {   prev_factor=r_x*r_x*r_x*r_y/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==11)
        {   curr_factor=r_x*r_x*r_x*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==11)
        {   prev_factor=r_x*r_x*r_x*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==12)
        {   curr_factor=r_y*r_y*r_y*r_x/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==12)
        {   prev_factor=r_y*r_y*r_y*r_x/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==13)
        {   curr_factor=r_y*r_y*r_y*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==13)
        {   prev_factor=r_y*r_y*r_y*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==14)
        {   curr_factor=r_z*r_z*r_z*r_x/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==14)
        {   prev_factor=r_z*r_z*r_z*r_x/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==15)
        {   curr_factor=r_z*r_z*r_z*r_y/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==15)
        {   prev_factor=r_z*r_z*r_z*r_y/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==16) 
        {   curr_factor=r_x*r_x*r_y*r_y/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==16)
        {   prev_factor=r_x*r_x*r_y*r_y/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==17)
        {   curr_factor=r_x*r_x*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==17)
        {   prev_factor=r_x*r_x*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==18)
        {   curr_factor=r_y*r_y*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==18)
        {   prev_factor=r_y*r_y*r_z*r_z/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==19)
        {   curr_factor=r_x*r_x*r_y*r_z/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==19)
        {   prev_factor=r_x*r_x*r_y*r_z/(kmod2*kmod2);//in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==20)
        {   curr_factor=r_y*r_y*r_x*r_z/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==20)
        {   prev_factor=r_y*r_y*r_x*r_z/(kmod2*kmod2);//in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==21)
        {   curr_factor=r_z*r_z*r_x*r_y/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto_prev==21)
        {   prev_factor=r_z*r_z*r_x*r_y/(kmod2*kmod2);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }

        if(mode_yamamoto==22)
        {   curr_factor=r_x*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==22)
        {   prev_factor=r_x*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }

        if(mode_yamamoto==23)
        {   curr_factor=r_y*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==23)
        {   prev_factor=r_y*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==24)
        {   curr_factor=r_z*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==24)
        {   prev_factor=r_z*pow(kmod2,-0.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==25)
        {   curr_factor=r_x*r_x*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==25)
        {   prev_factor=r_x*r_x*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==26)
        {   curr_factor=r_y*r_y*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==26)
        {   prev_factor=r_y*r_y*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==27)
        {   curr_factor=r_z*r_z*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==27)
        {   prev_factor=r_z*r_z*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==28)
        {   curr_factor=r_x*r_x*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==28)
        {   prev_factor=r_x*r_x*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==29)
        {   curr_factor=r_x*r_x*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==29)
        {   prev_factor=r_x*r_x*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==30)
        {   curr_factor=r_y*r_y*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==30)
        {   prev_factor=r_y*r_y*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==31)
        {   curr_factor=r_y*r_y*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==31)
        {   prev_factor=r_y*r_y*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==32)
        {   curr_factor=r_z*r_z*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==32)
        {   prev_factor=r_z*r_z*r_x*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==33)
        {   curr_factor=r_z*r_z*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==33)
        {   prev_factor=r_z*r_z*r_y*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==34)
        {   curr_factor=r_x*r_y*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto_prev==34)
        {   prev_factor=r_x*r_y*r_z*pow(kmod2,-1.5);//in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }

in[c]*=curr_factor/prev_factor;

              }

}


void fftw_yamamoto_skycut(int mode_yamamoto, double in[], double deltak_re[], double deltak_im[], int ngrid, double L1, double L2, int mode_mass_ass)
{
  fftw_complex *out;
 // double *out2;
 //double a1,a2;
  fftw_plan p;
  long int i,j,k;
 // double ii,jj,kk;
  long int c;
  double cx,cy,cz;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
  //printf("%ld\n",ngridtotr2c);
  long int index2;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
			if(i<ngrid/2+1)
			{
				kx[i]=i*1.0*(2.0*Pi/(L2-L1));
			}
			else
			{
				kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
			}
	}

    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ngridtotr2c));
//    if(out==NULL){printf("Error, grid-array (out fftw_yamamoto) could not be created. Exiting now...\n");exit(0);}
    fftw_plan_with_nthreads(omp_get_max_threads());
    p =  fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid,in,out,FFTW_ESTIMATE);
//   if(p==NULL){printf("Error, FFTW-plan (within fftw_yamamoto) could not be created. Exiting now...\n");exit(0);}
    fftw_execute(p);//FFT
    fftw_destroy_plan(p);

   #pragma omp parallel for private(c,i,j,k,cx,cy,cz,index2) shared(kx,deltak_re,deltak_im,out,ngrid,ngridtotr2c,mode_mass_ass,mode_yamamoto,Pi)
  for(index2=0;index2<ngridtotr2c;index2++)
  {

/*
ii=index2*1./(1.*(ngrid*ngrid/2.+ngrid));
i=(long int)(ii);
jj=(index2-i*(ngrid*ngrid/2.+ngrid))/(1.*(ngrid/2.+1));
j=(long int)(jj);
kk=index2-i*(ngrid*ngrid/2.+ngrid)-j*(ngrid/2.+1);
k=(long int)(kk);
*/
i=(long int)(index2*1./(1.*ngrid*ngrid/2.+ngrid*1.));
j=(long int)( (index2-i*(1.*ngrid*ngrid/2.+ngrid*1.))/(ngrid/2.+1.) );
//k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
k=index2-(i*(ngrid*ngrid+2*ngrid)+j*(ngrid+2))/2;

        cx=sin( kx[i]*Pi/(2.*kx[ngrid/2]) )/( kx[i]*Pi/(2.*kx[ngrid/2]) );
        cy=sin( kx[j]*Pi/(2.*kx[ngrid/2]) )/( kx[j]*Pi/(2.*kx[ngrid/2]) );
        cz=sin( kx[k]*Pi/(2.*kx[ngrid/2]) )/( kx[k]*Pi/(2.*kx[ngrid/2]) );
        if(kx[i]==0 || mode_mass_ass==0){cx=1.;}
        if(kx[j]==0 || mode_mass_ass==0){cy=1.;}
        if(kx[k]==0 || mode_mass_ass==0){cz=1.;}


if(mode_yamamoto==0){
          deltak_re[index2]=creal(out[index2])*pow(cx*cy*cz,-mode_mass_ass*1.);
          deltak_im[index2]=cimag(out[index2])*pow(cx*cy*cz,-mode_mass_ass*1.);           
}
if(mode_yamamoto==1){
           deltak_re[index2]=(creal(out[index2])*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
           deltak_im[index2]=(cimag(out[index2])*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==2){
deltak_re[index2]+=(creal(out[index2])*kx[i]*kx[j]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[i]*kx[j]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==3){
deltak_re[index2]+=(creal(out[index2])*kx[i]*kx[k]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[i]*kx[k]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==4){
deltak_re[index2]+=(creal(out[index2])*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==5){
deltak_re[index2]+=2.*(creal(out[index2])*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=2.*(cimag(out[index2])*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==6){
deltak_re[index2]+=(creal(out[index2])*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==7){
deltak_re[index2]=(creal(out[index2])*pow(kx[i],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]=(cimag(out[index2])*pow(kx[i],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
} 
if(mode_yamamoto==8){
deltak_re[index2]+=(creal(out[index2])*pow(kx[j],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=(cimag(out[index2])*pow(kx[j],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==9){
deltak_re[index2]+=(creal(out[index2])*pow(kx[k],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=(cimag(out[index2])*pow(kx[k],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==10){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[i],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[i],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==11){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[i],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[i],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==12){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[j],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[j],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==13){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[j],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[j],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==14){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[k],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[k],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==15){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[k],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[k],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==16){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[i]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[i]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==17){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[i]*kx[k],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[i]*kx[k],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==18){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[k]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[k]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==19){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[i],2))*kx[j]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[i],2))*kx[j]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==20){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[j],2))*kx[i]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[j],2))*kx[i]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==21){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[k],2))*kx[j]*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[k],2))*kx[j]*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==22){
deltak_re[index2]=(creal(out[index2])*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
deltak_im[index2]=(cimag(out[index2])*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
}
if(mode_yamamoto==23){
deltak_re[index2]+=(creal(out[index2])*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
deltak_im[index2]+=(cimag(out[index2])*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
}
if(mode_yamamoto==24){
deltak_re[index2]+=(creal(out[index2])*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
deltak_im[index2]+=(cimag(out[index2])*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-0.5);
}
if(mode_yamamoto==25){
deltak_re[index2]=(creal(out[index2])*kx[i]*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]=(cimag(out[index2])*kx[i]*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==26){
deltak_re[index2]+=(creal(out[index2])*kx[j]*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=(cimag(out[index2])*kx[j]*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==27){
deltak_re[index2]+=(creal(out[index2])*kx[k]*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=(cimag(out[index2])*kx[k]*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==28){
deltak_re[index2]+=3.*(creal(out[index2])*kx[i]*kx[i]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[i]*kx[i]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==29){
deltak_re[index2]+=3.*(creal(out[index2])*kx[i]*kx[i]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[i]*kx[i]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==30){
deltak_re[index2]+=3.*(creal(out[index2])*kx[j]*kx[j]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[j]*kx[j]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==31){
deltak_re[index2]+=3.*(creal(out[index2])*kx[j]*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[j]*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==32){
deltak_re[index2]+=3.*(creal(out[index2])*kx[k]*kx[k]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[k]*kx[k]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==33){
deltak_re[index2]+=3.*(creal(out[index2])*kx[k]*kx[k]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=3.*(cimag(out[index2])*kx[k]*kx[k]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}
if(mode_yamamoto==34){
deltak_re[index2]+=6.*(creal(out[index2])*kx[i]*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
deltak_im[index2]+=6.*(cimag(out[index2])*kx[i]*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1.5);
}




}

   
    fftw_free(out);
    free(kx);
fftw_cleanup();
//in not modified
//deltak_re filled
//deltak_im filled
//out is internally created and destroyed

} 




void loop_interlacing_skycut_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand, double L1, double L2, int ngrid, double alpha, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re0, double *deltak_im0, double *deltak_re2, double *deltak_im2,char *do_bispectrum2,char *name_density, char *type_of_input,char *type_of_survey)
{
FILE *f;
long int i,j,k,i_inter;
long int c;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b;
int i_yama;
double Sk;
long int index2;
double* delta_data;
double* delta_rand;
double* deltak_re0_b;
double* deltak_im0_b;
double* deltak_re2_b;
double* deltak_im2_b;

long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
double L2b,L1b;
double *kx;
double Pi=(4.*atan(1.));

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


  alpha=-alpha;

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

    delta_data = (double*) calloc(ngridtot, sizeof(double));

if( strcmp(type_of_input,"particles") == 0){

    printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);
    for(c=0; c<Ndata; c++)
    {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
    }

        delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
        }
    printf("Ok!\n");

    printf("Performing FFTs ...");
 #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);


}else{
     printf("Reading densities ...");
f=fopen(name_density,"r");if(f==NULL){printf("Error reading density file %s. Exiting now.\n",name_density);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&delta_data[c]);
        delta_data[c]=delta_data[c]*pow((L2-L1)/ngrid,3);
     }
fclose(f);

}


        if(i_inter==1){
       fftw_yamamoto_skycut(0, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);

if(strcmp(do_bispectrum2, "yes") == 0 && strcmp(type_of_survey,"cutsky") == 0)
{
   for(i_yama=1;i_yama<=6;i_yama++)
   {
         modify_input_for_yamamoto(i_yama,i_yama-1, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2, deltak_im2, ngrid,L1b,L2b, mode_correction);
   }
}

  }
  else{
        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);

if(strcmp(do_bispectrum2, "yes") == 0  && strcmp(type_of_survey,"cutsky") == 0)
{

       deltak_re2_b= (double*) calloc(ngridtotr2c,sizeof(double));
       deltak_im2_b= (double*) calloc(ngridtotr2c,sizeof(double));

  for(i_yama=1;i_yama<=6;i_yama++)
  {
         modify_input_for_yamamoto(i_yama,i_yama-1, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2_b, deltak_im2_b, ngrid,L1b,L2b, mode_correction);
  }



}

#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re0_b,deltak_im0_b,deltak_re2,deltak_im2,deltak_re2_b,deltak_im2_b,kx,L2,L1,i_inter,Ninterlacing,do_bispectrum2,type_of_survey)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(long int)(index2/(ngrid*ngrid/2+ngrid));
j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
Sk=kx[i]+ kx[j]+ kx[k];

phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;

if(strcmp(do_bispectrum2, "yes") == 0  && strcmp(type_of_survey,"cutsky") == 0)
{
new_deltak_re2_b=deltak_re2_b[index2]*phase_cos-deltak_im2_b[index2]*phase_sin;
new_deltak_im2_b=deltak_re2_b[index2]*phase_sin+deltak_im2_b[index2]*phase_cos;

deltak_re2[index2]+=new_deltak_re2_b;
deltak_im2[index2]+=new_deltak_im2_b;

}

}

        free(deltak_re0_b);
        free(deltak_im0_b);
if(strcmp(do_bispectrum2, "yes") == 0  && strcmp(type_of_survey,"cutsky") == 0)
{
        free(deltak_re2_b);
        free(deltak_im2_b);
}
  }

        free(delta_data);
printf("Ok!\n");


}//loop interlacing
free(kx);
//deltak_re and deltak_im ready for bispectrum computation

}


void loop_interlacing_skycut(double kmin,double kmax,int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand,double *pos_xB, double *pos_yB, double *pos_zB, double *weightB, long int NdataB, double *pos_x_randB, double *pos_y_randB, double *pos_z_randB, double *weight_randB,long int NrandB, double L1, double L2, int ngrid,  double P_shot_noise,double P_shot_noiseB, double bin_ps, double I22,double I22B, double alpha,double alphaB, int mode_correction, int n_lines_parallel, char *binning_type, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out,char *name_psAB_out,char *name_psBB_out, char *type_of_mass_assigment, char *do_bispectrum,char *type_of_code, char *output_density, char *name_density,char *name_densityB, char *type_of_input,char *type_of_survey,char *file_for_mu,int Nmu,char *write_kvectors, char *name_ps_kvectors)
{
FILE *f;
char type[200];
long int i,j,k,i_inter;
long int c;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b;
double new_deltak_re0_bB,new_deltak_im0_bB,new_deltak_re2_bB,new_deltak_im2_bB,new_deltak_re4_bB,new_deltak_im4_bB;

double new_deltak_re1_b,new_deltak_im1_b,new_deltak_re3_b,new_deltak_im3_b;
double new_deltak_re1_bB,new_deltak_im1_bB,new_deltak_re3_bB,new_deltak_im3_bB;
int i_yama,i_yama2;
double Sk;
long int index2;
  //delta(k) pointers
  double* delta_data;
  double* delta_rand;
    
  double* delta_dataB;
  double* delta_randB;
    
  double* deltak_re0;
  double* deltak_im0;
  double* deltak_re1;
  double* deltak_im1;
  double* deltak_re2;
  double* deltak_im2;
  double* deltak_re3;
  double* deltak_im3;
  double* deltak_re4;
  double* deltak_im4;
    
  double* deltak_re0B;
  double* deltak_im0B;
  double* deltak_re1B;
  double* deltak_im1B;
  double* deltak_re2B;
  double* deltak_im2B;
  double* deltak_re3B;
  double* deltak_im3B;
  double* deltak_re4B;
  double* deltak_im4B;
    
  double* deltak_re0_b;
  double* deltak_im0_b;
  double* deltak_re1_b;
  double* deltak_im1_b;
  double* deltak_re2_b;
  double* deltak_im2_b;
  double* deltak_re3_b;
  double* deltak_im3_b;
  double* deltak_re4_b;
  double* deltak_im4_b;
    
    double* deltak_re0_bB;
    double* deltak_im0_bB;
    double* deltak_re1_bB;
    double* deltak_im1_bB;
    double* deltak_re2_bB;
    double* deltak_im2_bB;
    double* deltak_re3_bB;
    double* deltak_im3_bB;
    double* deltak_re4_bB;
    double* deltak_im4_bB;
    
long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)

int delta1_sw,delta2_sw,delta3_sw,delta4_sw;

  double L2b,L1b, d;

  double *kx;
  double Pi=(4.*atan(1.));

delta1_sw=0;
delta2_sw=0;
delta3_sw=0;
delta4_sw=0;

if(strcmp(Hexadecapole_type, "L2L2") == 0){delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L0L4") == 0){delta4_sw=1;delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L1L3") == 0){delta1_sw=1;delta3_sw=1;}

if(strcmp(Quadrupole_type, "L0L2") == 0){delta2_sw=1;}
if(strcmp(Quadrupole_type, "L1L1") == 0){delta1_sw=1;}

if(strcmp(do_odd_multipoles,"yes") == 0){
if(strcmp(Octopole_type, "L0L3") == 0){delta3_sw=1;delta1_sw=1;}
if(strcmp(Octopole_type, "L1L2") == 0){delta1_sw=1;delta2_sw=1;}
}


if(strcmp(do_anisotropy,"no") == 0 || strcmp(type_of_survey,"periodicFKP") == 0){delta1_sw=0;delta2_sw=0;delta3_sw=0;delta4_sw=0;}

//printf("%ld %ld %ld %ld\n",delta1_sw,delta2_sw,delta3_sw,delta4_sw);

alpha=-alpha;
if(strcmp(type_of_code, "rusticoX") == 0){alphaB=-alphaB;}

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));   

if( strcmp(type_of_input,"particles") == 0){

    printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);


    for(c=0; c<Ndata; c++)
    {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
    }

       delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
        }
    printf("Ok!\n");

    printf("Performing FFTs ...");
 #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);

if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_density,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_data[c]*pow((L2-L1)/ngrid,-3);
        fprintf(f,"%e\n",d); // this is (Ng(x)-alpha*Nr(x) ) /Vcell.
     }
fclose(f);
}

}else{
     printf("Reading densities ...");
f=fopen(name_density,"r");if(f==NULL){printf("Error reading density file %s. Exiting now.\n",name_density);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&d);
        delta_data[c]=d*pow((L2-L1)/ngrid,3);
     }
fclose(f);

}

    if(strcmp(type_of_code, "rusticoX") == 0)
    {
          delta_dataB = (double*) calloc(ngridtot, sizeof(double));

           if( strcmp(type_of_input,"particles") == 0){

           printf("Assigning particles-B to the grid (Iteration %ld) ...", i_inter);


           for(c=0; c<NdataB; c++)
           {
              if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
           }

              delta_randB = (double*) calloc(ngridtot, sizeof(double));
              for(c=0; c<NrandB; c++)
               {
                       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
                       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
                       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
                       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
                       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
                       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c], L2b, L1b, ngrid);}
               }
           printf("Ok!\n");

           printf("Performing FFTs ...");
        #pragma omp parallel for private(c) shared(ngrid,ngridtot,alphaB,delta_dataB,delta_randB)
               for(c=0;c<ngridtot;c++)
               {
                                        delta_dataB[c]+=alphaB*delta_randB[c];
               }
               free(delta_randB);

if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_densityB,"w");
     for(c=0;c<ngridtot;c++)
     {
        fprintf(f,"%e\n",delta_dataB[c]*pow((L2-L1)/ngrid,-3));
     }
fclose(f);
}

}else{
     printf("Reading densities ...");
     f=fopen(name_densityB,"r");if(f==NULL){printf("Error reading density file %s. Exiting now.\n",name_densityB);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&d);
        delta_dataB[c]=d*pow((L2-L1)/ngrid,3);
     }
fclose(f);

}
        
    }
    
    
  if(i_inter==1){
        deltak_re0= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0=(double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);
      if(strcmp(type_of_code, "rusticoX") == 0)
      {
          deltak_re0B= (double*) calloc(ngridtotr2c,sizeof(double));
          deltak_im0B=(double*) calloc(ngridtotr2c,sizeof(double));
          fftw_yamamoto_skycut(0, delta_dataB, deltak_re0B, deltak_im0B, ngrid,L1b,L2b, mode_correction);
    
      }

if(strcmp(do_anisotropy,"yes") ==0 && strcmp(type_of_survey,"cutsky") == 0 ){

if(delta2_sw==1){
         deltak_re2=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im2=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re2B=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2B=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
}

if(delta1_sw==1){
         deltak_re1=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im1=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re1B=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im1B=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
    
}

if(delta3_sw == 1){
         deltak_re3=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im3=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re3B=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im3B=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
    
}

if(delta4_sw==1){

         deltak_re4= (double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im4= (double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re4B= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4B= (double*) calloc(ngridtotr2c,sizeof(double));
    }
    
    
}

i_yama2=0;
if(delta2_sw==1){
for(i_yama=1;i_yama<=6;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2, deltak_im2, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re2B, deltak_im2B, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}

if(delta4_sw==1){

for(i_yama=7;i_yama<=21;i_yama++)
{

         modify_input_for_yamamoto(i_yama,i_yama2, delta_data, ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4, deltak_im4, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB, ngrid,L1,L2);
    fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re4B, deltak_im4B, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}

}

if(delta1_sw==1){
for(i_yama=22;i_yama<=24;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re1, deltak_im1, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re1B, deltak_im1B, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}

if(delta3_sw==1){
for(i_yama=25;i_yama<=34;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re3, deltak_im3, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re3B, deltak_im3B, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}


}//do_aniso


  }//iteration 1
  else{  

        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);
      if(strcmp(type_of_code, "rusticoX") == 0)
      {
          deltak_re0_bB= (double*) calloc(ngridtotr2c,sizeof(double));
          deltak_im0_bB= (double*) calloc(ngridtotr2c,sizeof(double));
          fftw_yamamoto_skycut(0, delta_dataB, deltak_re0_bB, deltak_im0_bB, ngrid,L1b,L2b, mode_correction);
      }

if(strcmp(do_anisotropy,"yes") ==0 && strcmp(type_of_survey,"cutsky") == 0){

if(delta2_sw==1){
         deltak_re2_b=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im2_b=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re2_bB=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2_bB=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
}

if(delta1_sw==1){
         deltak_re1_b=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im1_b=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re1_bB=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im1_bB=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
}

if(delta3_sw==1){
         deltak_re3_b=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im3_b=(double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re3_bB=(double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im3_bB=(double*) calloc(ngridtotr2c,sizeof(double));
    }
    
}

if(delta4_sw==1){

         deltak_re4_b= (double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im4_b= (double*) calloc(ngridtotr2c,sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        deltak_re4_bB= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4_bB= (double*) calloc(ngridtotr2c,sizeof(double));
    }
    
}


i_yama2=0;
if(delta2_sw==1){
for(i_yama=1;i_yama<=6;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2_b, deltak_im2_b, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re2_bB, deltak_im2_bB, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}

if(delta4_sw==1){

for(i_yama=7;i_yama<=21;i_yama++)
{

         modify_input_for_yamamoto(i_yama,i_yama2, delta_data, ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4_b, deltak_im4_b, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB, ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re4_bB, deltak_im4_bB, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}

}

if(delta1_sw==1){
for(i_yama=22;i_yama<=24;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re1_b, deltak_im1_b, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re1_bB, deltak_im1_bB, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}

if(delta3_sw==1){
for(i_yama=25;i_yama<=34;i_yama++)
{
         modify_input_for_yamamoto(i_yama,i_yama2, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re3_b, deltak_im3_b, ngrid,L1b,L2b, mode_correction);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        modify_input_for_yamamoto(i_yama,i_yama2, delta_dataB,ngrid,L1,L2);
        fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re3_bB, deltak_im3_bB, ngrid,L1b,L2b, mode_correction);
    }
i_yama2=i_yama;
}
}

}//do-aniso


#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b,new_deltak_re1_b,new_deltak_im1_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re3_b,new_deltak_im3_b,new_deltak_re4_b,new_deltak_im4_b,new_deltak_re0_bB,new_deltak_im0_bB,new_deltak_re1_bB,new_deltak_im1_bB,new_deltak_re2_bB,new_deltak_im2_bB,new_deltak_re3_bB,new_deltak_im3_bB,new_deltak_re4_bB,new_deltak_im4_bB) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re1,deltak_im1,deltak_re2,deltak_im2,deltak_re3,deltak_im3,deltak_re4,deltak_im4,deltak_re0_b,deltak_im0_b,deltak_re1_b,deltak_im1_b,deltak_re2_b,deltak_im2_b,deltak_re3_b,deltak_im3_b,deltak_re4_b,deltak_im4_b,kx,L2,L1,i_inter,Ninterlacing,delta1_sw,delta2_sw,delta3_sw,delta4_sw,do_anisotropy,type_of_code,deltak_re0B,deltak_im0B,deltak_re1B,deltak_im1B,deltak_re2B,deltak_im2B,deltak_re3B,deltak_im3B,deltak_re4B,deltak_im4B,deltak_re0_bB,deltak_im0_bB,deltak_re1_bB,deltak_im1_bB,deltak_re2_bB,deltak_im2_bB,deltak_re3_bB,deltak_im3_bB,deltak_re4_bB,deltak_im4_bB,type_of_survey)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(long int)(index2/(ngrid*ngrid/2+ngrid));
j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

Sk=kx[i]+ kx[j]+ kx[k];

phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        new_deltak_re0_bB=deltak_re0_bB[index2]*phase_cos-deltak_im0_bB[index2]*phase_sin;
        new_deltak_im0_bB=deltak_re0_bB[index2]*phase_sin+deltak_im0_bB[index2]*phase_cos;

        deltak_re0B[index2]+=new_deltak_re0_bB;
        deltak_im0B[index2]+=new_deltak_im0_bB;
    }

if(strcmp(do_anisotropy,"yes") ==0 && strcmp(type_of_survey,"cutsky") == 0){

if(delta1_sw==1){

new_deltak_re1_b=deltak_re1_b[index2]*phase_cos-deltak_im1_b[index2]*phase_sin;
new_deltak_im1_b=deltak_re1_b[index2]*phase_sin+deltak_im1_b[index2]*phase_cos;

deltak_re1[index2]+=new_deltak_re1_b;
deltak_im1[index2]+=new_deltak_im1_b;
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        new_deltak_re1_bB=deltak_re1_bB[index2]*phase_cos-deltak_im1_bB[index2]*phase_sin;
        new_deltak_im1_bB=deltak_re1_bB[index2]*phase_sin+deltak_im1_bB[index2]*phase_cos;

        deltak_re1B[index2]+=new_deltak_re1_bB;
        deltak_im1B[index2]+=new_deltak_im1_bB;
    }

}

if(delta2_sw==1){

new_deltak_re2_b=deltak_re2_b[index2]*phase_cos-deltak_im2_b[index2]*phase_sin;
new_deltak_im2_b=deltak_re2_b[index2]*phase_sin+deltak_im2_b[index2]*phase_cos;

deltak_re2[index2]+=new_deltak_re2_b;
deltak_im2[index2]+=new_deltak_im2_b;
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        new_deltak_re2_bB=deltak_re2_bB[index2]*phase_cos-deltak_im2_bB[index2]*phase_sin;
        new_deltak_im2_bB=deltak_re2_bB[index2]*phase_sin+deltak_im2_bB[index2]*phase_cos;

        deltak_re2B[index2]+=new_deltak_re2_bB;
        deltak_im2B[index2]+=new_deltak_im2_bB;
    }

}

if(delta3_sw==1){

new_deltak_re3_b=deltak_re3_b[index2]*phase_cos-deltak_im3_b[index2]*phase_sin;
new_deltak_im3_b=deltak_re3_b[index2]*phase_sin+deltak_im3_b[index2]*phase_cos;

deltak_re3[index2]+=new_deltak_re3_b;
deltak_im3[index2]+=new_deltak_im3_b;
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        new_deltak_re3_bB=deltak_re3_bB[index2]*phase_cos-deltak_im3_bB[index2]*phase_sin;
        new_deltak_im3_bB=deltak_re3_bB[index2]*phase_sin+deltak_im3_bB[index2]*phase_cos;

        deltak_re3B[index2]+=new_deltak_re3_bB;
        deltak_im3B[index2]+=new_deltak_im3_bB;
    }

}

if(delta4_sw==1){

new_deltak_re4_b=deltak_re4_b[index2]*phase_cos-deltak_im4_b[index2]*phase_sin;
new_deltak_im4_b=deltak_re4_b[index2]*phase_sin+deltak_im4_b[index2]*phase_cos;

deltak_re4[index2]+=new_deltak_re4_b;
deltak_im4[index2]+=new_deltak_im4_b;
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        new_deltak_re4_bB=deltak_re4_bB[index2]*phase_cos-deltak_im4_bB[index2]*phase_sin;
        new_deltak_im4_bB=deltak_re4_bB[index2]*phase_sin+deltak_im4_bB[index2]*phase_cos;

        deltak_re4B[index2]+=new_deltak_re4_bB;
        deltak_im4B[index2]+=new_deltak_im4_bB;
    }

}

}//do aniso


  }//for index2

        free(deltak_re0_b);
        free(deltak_im0_b);
      if(strcmp(type_of_code, "rusticoX") == 0)
      {
          free(deltak_re0_bB);
          free(deltak_im0_bB);
      }

if(strcmp(do_anisotropy,"yes") ==0 && strcmp(type_of_survey,"periodicFKP") == 0){

if(delta1_sw==1){
        free(deltak_re1_b);
        free(deltak_im1_b);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re1_bB);
        free(deltak_im1_bB);
    }
    
}

if(delta2_sw==1){
        free(deltak_re2_b);
        free(deltak_im2_b);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re2_bB);
        free(deltak_im2_bB);
    }
    
}

if(delta3_sw==1){
        free(deltak_re3_b);
        free(deltak_im3_b);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re3_bB);
        free(deltak_im3_bB);
    }
    
}

if(delta4_sw==1){
         free(deltak_re4_b);
         free(deltak_im4_b);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re4_bB);
        free(deltak_im4_bB);
    }
    
}

}


}//else-iteration1
        free(delta_data);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(delta_dataB);
    }
printf("Ok!\n");

if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0)
{
if( strcmp(type_of_input,"particles") == 0){

printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);
free(pos_x_rand);
free(pos_y_rand);
free(pos_z_rand);
free(weight_rand);}

    if(strcmp(type_of_code, "rusticoX") == 0)
    {
 
if( strcmp(type_of_input,"particles") == 0){
        free(pos_xB);
        free(pos_yB);
        free(pos_zB);
        free(weightB);
        free(pos_x_randB);
        free(pos_y_randB);
        free(pos_z_randB);
        free(weight_randB);}
    }
}


}//loop interlacing

free(kx);


if( strcmp( output_density,"corrected") == 0 )//change write_PS do not erase them
{
d=pow((L2-L1),-3);
generate_deltax_field(deltak_re0, deltak_im0, ngrid, Ninterlacing, name_density, d);

if(strcmp(type_of_code, "rusticoX") == 0){generate_deltax_field(deltak_re0B, deltak_im0B, ngrid,  Ninterlacing, name_densityB, d);}

}

    if(strcmp(type_of_code, "rustico") == 0){
        printf("Writing Power Spectrum output %s...",name_ps_out);
        
    }
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        printf("Writing Power Spectrum output %s %s %s...",name_ps_out,name_psAB_out,name_psBB_out);

    }

sprintf(type,"FFT");

if(strcmp(type_of_survey,"cutsky") == 0){
    write_power_spectrum_skyscuts(kmin,kmax,NULL, NULL, NULL, deltak_re0, deltak_im0,deltak_re1, deltak_im1, deltak_re2, deltak_im2, deltak_re3, deltak_im3, deltak_re4, deltak_im4,deltak_re0B, deltak_im0B,deltak_re1B, deltak_im1B, deltak_re2B, deltak_im2B, deltak_re3B, deltak_im3B, deltak_re4B, deltak_im4B, bin_ps, ngrid, 0, L1, L2, I22,I22B, name_ps_out,name_psAB_out,name_psBB_out, P_shot_noise,P_shot_noiseB, binning_type, do_anisotropy, do_odd_multipoles,  Quadrupole_type, Octopole_type, Hexadecapole_type, type ,Ninterlacing,type_of_code,write_kvectors, name_ps_kvectors);
}
if(strcmp(type_of_survey,"periodicFKP") == 0){

if(Nmu==1){
write_power_spectrum_periodic(kmin,kmax,deltak_re0, deltak_im0, deltak_re0B,deltak_im0B,bin_ps, ngrid, L1, L2, Ninterlacing, name_ps_out, name_psAB_out, name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,do_odd_multipoles,do_anisotropy,type_of_code,I22,I22B,write_kvectors, name_ps_kvectors);
}
if(Nmu>1){
write_power_spectrum_periodic2D(kmin,kmax,deltak_re0,deltak_im0,deltak_re0B,deltak_im0B,bin_ps,Nmu, ngrid,L1,L2,Ninterlacing,name_ps_out,name_psAB_out,name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,file_for_mu,type_of_code,I22,I22B,write_kvectors, name_ps_kvectors);
}

}

//free deltas

free(deltak_re0);
free(deltak_im0);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re0B);
        free(deltak_im0B);
    }
if(strcmp(do_anisotropy,"yes") ==0 && strcmp(type_of_survey,"cutsky") == 0){

if(delta2_sw==1){
          free(deltak_re2);
          free(deltak_im2);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re2B);
        free(deltak_im2B);
    }
    
}

if(delta1_sw==1){
          free(deltak_re1);
          free(deltak_im1);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re1B);
        free(deltak_im1B);
    }
    
}

if(delta3_sw==1){
          free(deltak_re3);
          free(deltak_im3);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re3B);
        free(deltak_im3B);
    }
    
}

if(delta4_sw==1){

          free(deltak_re4);
          free(deltak_im4);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_re4B);
        free(deltak_im4B);
    }
    
}
}

printf("Ok!\n");

}

void loop_interlacing_skycut2(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand,double *pos_xB, double *pos_yB, double *pos_zB, double *weightB, long int NdataB, double *pos_x_randB, double *pos_y_randB, double *pos_z_randB, double *weight_randB, long int NrandB, double L1, double L2, int ngrid,  double P_shot_noise,double P_shot_noiseB, double bin_ps, double I22,double I22B, double alpha,double alphaB, int mode_correction, int n_lines_parallel, char *binning_type, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out,char *name_psAB_out,char *name_psBB_out, char *type_of_mass_assigment, char *do_bispectrum,char *type_of_code, char *output_density,  char *name_density,char *name_densityB, char *type_of_survey, char *file_for_mu, int Nmu,char *write_kvectors, char *name_ps_kvectors)
{
FILE *f;
char type[200];
long int i,j,k,i_inter;
long int c;
//int i_yama_max;
//i_yama_max;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b,new_deltak_re1_b,new_deltak_im1_b,new_deltak_re3_b,new_deltak_im3_b;
    double new_deltak_re0_bB,new_deltak_im0_bB,new_deltak_re2_bB,new_deltak_im2_bB,new_deltak_re4_bB,new_deltak_im4_bB,new_deltak_re1_bB,new_deltak_im1_bB,new_deltak_re3_bB,new_deltak_im3_bB;
int i_yama;

int delta1_sw=0;
int delta2_sw=0;
int delta3_sw=0;
int delta4_sw=0;

if(strcmp(Hexadecapole_type, "L2L2") == 0){delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L0L4") == 0){delta4_sw=1;delta2_sw=1;}
if(strcmp(Hexadecapole_type, "L1L3") == 0){delta1_sw=1;delta3_sw=1;}

if(strcmp(Quadrupole_type, "L0L2") == 0){delta2_sw=1;}
if(strcmp(Quadrupole_type, "L1L1") == 0){delta1_sw=1;}

if(strcmp(do_odd_multipoles,"yes") == 0){
if(strcmp(Octopole_type, "L0L3") == 0){delta3_sw=1;delta1_sw=1;}
if(strcmp(Octopole_type, "L1L2") == 0){delta1_sw=1;delta2_sw=1;}
}

if(strcmp(do_anisotropy,"no") == 0 || strcmp(type_of_survey,"periodicFKP") == 0){delta1_sw=0;delta2_sw=0;delta3_sw=0;delta4_sw=0;}


alpha=-alpha;
if(strcmp(type_of_code, "rusticoX") == 0){alphaB=-alphaB;}

  //delta(k) pointers
  double* delta_data;
  double* delta_rand;
    
  double* delta_dataB;
  double* delta_randB;
    
  double* deltak_re0;
  double* deltak_im0;
  double* deltak_re1;
  double* deltak_im1;
  double* deltak_re2;
  double* deltak_im2;
  double* deltak_re3;
  double* deltak_im3;
  double* deltak_re4;
  double* deltak_im4;
    
    double* deltak_re0B;
    double* deltak_im0B;
    double* deltak_re1B;
    double* deltak_im1B;
    double* deltak_re2B;
    double* deltak_im2B;
    double* deltak_re3B;
    double* deltak_im3B;
    double* deltak_re4B;
    double* deltak_im4B;
    
  double* deltak_re0_b;
  double* deltak_im0_b;
  double* deltak_re1_b;
  double* deltak_im1_b;
  double* deltak_re2_b;
  double* deltak_im2_b;
  double* deltak_re3_b;
  double* deltak_im3_b;
  double* deltak_re4_b;
  double* deltak_im4_b;
    
    double* deltak_re0_bB;
    double* deltak_im0_bB;
    double* deltak_re1_bB;
    double* deltak_im1_bB;
    double* deltak_re2_bB;
    double* deltak_im2_bB;
    double* deltak_re3_bB;
    double* deltak_im3_bB;
    double* deltak_re4_bB;
    double* deltak_im4_bB;
    
long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));
long int index2;
double weight_pos;

  double L2b,L1b,d;
  
  double *kx;
  double Pi=(4.*atan(1.));

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

if(i_inter==1)
{
        printf("Assigning particles to the grid and FFTing (Iteration %ld) ...", i_inter);

          
          for(i_yama=0;i_yama<=34;i_yama++)
          {

if(strcmp(do_anisotropy,"no") == 0 && i_yama==1){break;}

if(strcmp(type_of_survey,"periodicFKP") == 0 && i_yama==1){break;}

if(i_yama==1 && delta2_sw==0){i_yama=7;}
 
if(i_yama==7 && delta4_sw==0){i_yama=22;}

if(i_yama==22 && delta1_sw==0){i_yama=25;}

if(i_yama==25 && delta3_sw==0){i_yama=35;}

if(i_yama==35){break;}

            
        delta_data = (double*) calloc(ngridtot, sizeof(double));

        for(c=0; c<Ndata; c++)
        {

if(i_yama==0){weight_pos=1;}

if(i_yama==1){weight_pos=pos_x[c]*pos_x[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==2){weight_pos=pos_x[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==3){weight_pos=pos_x[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==4){weight_pos=pos_y[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==5){weight_pos=pos_y[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==6){weight_pos=pos_z[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}

if(i_yama==7){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==8){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==9){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==10){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==11){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==12){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==13){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==14){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==15){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==16){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==17){weight_pos=pow(pos_x[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==18){weight_pos=pow(pos_y[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==19){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==20){weight_pos=pow(pos_y[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==21){weight_pos=pow(pos_z[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}

if(i_yama==22){weight_pos=pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==23){weight_pos=pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==24){weight_pos=pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==25){weight_pos=pos_x[c]*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==26){weight_pos=pos_y[c]*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==27){weight_pos=pos_z[c]*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==28){weight_pos=pos_x[c]*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==29){weight_pos=pos_x[c]*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==30){weight_pos=pos_y[c]*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==31){weight_pos=pos_y[c]*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==32){weight_pos=pos_z[c]*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==33){weight_pos=pos_z[c]*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==34){weight_pos=pos_x[c]*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}

if(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]==0){weight_pos=0;}

       weight_pos=weight_pos*weight[c];

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
    } 

        delta_rand = (double*) calloc(ngridtot, sizeof(double));

      for(c=0;c<Nrand;c++)
      {
if(i_yama==0){weight_pos=1.;}
if(i_yama==1){weight_pos=pos_x_rand[c]*pos_x_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==2){weight_pos=pos_x_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==3){weight_pos=pos_x_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==4){weight_pos=pos_y_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==5){weight_pos=pos_y_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==6){weight_pos=pos_z_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==7){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==8){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==9){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==10){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==11){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==12){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==13){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==14){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==15){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==16){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==17){weight_pos=pow(pos_x_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==18){weight_pos=pow(pos_y_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==19){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==20){weight_pos=pow(pos_y_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==21){weight_pos=pow(pos_z_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==22){weight_pos=pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==23){weight_pos=pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==24){weight_pos=pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==25){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==26){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==27){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==28){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==29){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==30){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==31){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==32){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==33){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==34){weight_pos=pos_x_rand[c]*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}


if(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]==0){weight_pos=0;}
      weight_pos=weight_pos*weight_rand[c];

  
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}

        }

   

         
        #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++) 
        { 
                                 delta_data[c]+=alpha*delta_rand[c];
        } 
        free(delta_rand);

if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_density,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_data[c]*pow((L2-L1)/ngrid,-3);
        fprintf(f,"%e\n",d);
     }
fclose(f);
}


if(i_yama==0){
        deltak_re0= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==1)
{
        deltak_re2= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==7)
{
        deltak_re4= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4= (double*) calloc(ngridtotr2c,sizeof(double));
}

if(i_yama==22)
{
        deltak_re1= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im1= (double*) calloc(ngridtotr2c,sizeof(double));
}

if(i_yama==25)
{
        deltak_re3= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im3= (double*) calloc(ngridtotr2c,sizeof(double));

}


if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2, deltak_im2, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=7 && i_yama<=21){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4, deltak_im4, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=22 && i_yama<=24){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re1, deltak_im1, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=25 && i_yama<=34){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re3, deltak_im3, ngrid,L1b,L2b, mode_correction);}


        free(delta_data);
              
              if(strcmp(type_of_code, "rusticoX") == 0){
                  
                          delta_dataB = (double*) calloc(ngridtot, sizeof(double));

                          for(c=0; c<NdataB; c++)
                          {

                  if(i_yama==0){weight_pos=1;}

                  if(i_yama==1){weight_pos=pos_xB[c]*pos_xB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
                  if(i_yama==2){weight_pos=pos_xB[c]*pos_yB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
                  if(i_yama==3){weight_pos=pos_xB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
                  if(i_yama==4){weight_pos=pos_yB[c]*pos_yB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
                  if(i_yama==5){weight_pos=pos_yB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
                  if(i_yama==6){weight_pos=pos_zB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}

                  if(i_yama==7){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==8){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==9){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==10){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==11){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==12){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==13){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==14){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==15){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==16){weight_pos=pow(pos_xB[c],2)*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==17){weight_pos=pow(pos_xB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==18){weight_pos=pow(pos_yB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==19){weight_pos=pow(pos_xB[c],2)*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==20){weight_pos=pow(pos_yB[c],2)*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
                  if(i_yama==21){weight_pos=pow(pos_zB[c],2)*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}

                  if(i_yama==22){weight_pos=pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
                  if(i_yama==23){weight_pos=pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
                  if(i_yama==24){weight_pos=pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
                  if(i_yama==25){weight_pos=pos_xB[c]*pos_xB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==26){weight_pos=pos_yB[c]*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==27){weight_pos=pos_zB[c]*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==28){weight_pos=pos_xB[c]*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==29){weight_pos=pos_xB[c]*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==30){weight_pos=pos_yB[c]*pos_yB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==31){weight_pos=pos_yB[c]*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==32){weight_pos=pos_zB[c]*pos_zB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==33){weight_pos=pos_zB[c]*pos_zB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
                  if(i_yama==34){weight_pos=pos_xB[c]*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}

                  if(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]==0){weight_pos=0;}

                         weight_pos=weight_pos*weightB[c];

                         if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                         if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                         if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                         if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                         if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                         if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight_pos, L2b, L1b, ngrid);}
                      }

                          delta_randB = (double*) calloc(ngridtot, sizeof(double));

                        for(c=0;c<NrandB;c++)
                        {
                  if(i_yama==0){weight_pos=1.;}
                  if(i_yama==1){weight_pos=pos_x_randB[c]*pos_x_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
                  if(i_yama==2){weight_pos=pos_x_randB[c]*pos_y_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
                  if(i_yama==3){weight_pos=pos_x_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
                  if(i_yama==4){weight_pos=pos_y_randB[c]*pos_y_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
                  if(i_yama==5){weight_pos=pos_y_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
                  if(i_yama==6){weight_pos=pos_z_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
            if(i_yama==7){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);} if(i_yama==8){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                            if(i_yama==9){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                            if(i_yama==10){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==11){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==12){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==13){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==14){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==15){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==16){weight_pos=pow(pos_x_randB[c],2)*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==17){weight_pos=pow(pos_x_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==18){weight_pos=pow(pos_y_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==19){weight_pos=pow(pos_x_randB[c],2)*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==20){weight_pos=pow(pos_y_randB[c],2)*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==21){weight_pos=pow(pos_z_randB[c],2)*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
                                   if(i_yama==22){weight_pos=pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
                                   if(i_yama==23){weight_pos=pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
                                   if(i_yama==24){weight_pos=pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
                                   if(i_yama==25){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==26){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==27){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==28){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==29){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==30){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==31){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==32){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==33){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
                                   if(i_yama==34){weight_pos=pos_x_randB[c]*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}


                                   if(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]==0){weight_pos=0;}
                                         weight_pos=weight_pos*weight_randB[c];

                                     
                                                   if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}
                                                   if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}
                                                   if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}
                                                   if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}
                                                   if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}
                                                   if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_pos, L2b, L1b, ngrid);}

                                           }

                     

                           
                          #pragma omp parallel for private(c) shared(ngrid,ngridtot,alphaB,delta_dataB,delta_randB)
                          for(c=0;c<ngridtot;c++)
                          {
                                                   delta_dataB[c]+=alphaB*delta_randB[c];
                          }
                          free(delta_randB);

if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_densityB,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_dataB[c]*pow((L2-L1)/ngrid,-3);
        fprintf(f,"%e\n",d);
     }
fclose(f);
}


                  if(i_yama==0){
                          deltak_re0B= (double*) calloc(ngridtotr2c,sizeof(double));
                          deltak_im0B= (double*) calloc(ngridtotr2c,sizeof(double));
                  }
                  if(i_yama==1)
                  {
                          deltak_re2B= (double*) calloc(ngridtotr2c,sizeof(double));
                          deltak_im2B= (double*) calloc(ngridtotr2c,sizeof(double));
                  }
                  if(i_yama==7)
                  {
                          deltak_re4B= (double*) calloc(ngridtotr2c,sizeof(double));
                          deltak_im4B= (double*) calloc(ngridtotr2c,sizeof(double));
                  }

                  if(i_yama==22)
                  {
                          deltak_re1B= (double*) calloc(ngridtotr2c,sizeof(double));
                          deltak_im1B= (double*) calloc(ngridtotr2c,sizeof(double));
                  }

                  if(i_yama==25)
                  {
                          deltak_re3B= (double*) calloc(ngridtotr2c,sizeof(double));
                          deltak_im3B= (double*) calloc(ngridtotr2c,sizeof(double));

                  }


                  if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re0B, deltak_im0B, ngrid,L1b,L2b, mode_correction);}
                  if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re2B, deltak_im2B, ngrid,L1b,L2b, mode_correction);}
                  if(i_yama>=7 && i_yama<=21){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re4B, deltak_im4B, ngrid,L1b,L2b, mode_correction);}
                  if(i_yama>=22 && i_yama<=24){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re1B, deltak_im1B, ngrid,L1b,L2b, mode_correction);}
                  if(i_yama>=25 && i_yama<=34){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re3B, deltak_im3B, ngrid,L1b,L2b, mode_correction);}


                          free(delta_dataB);
                  
                  
              }
            
}

    printf("Ok!\n");

  }
  else{  
        printf("Assigning particles to the grid and FFTing (Iteration %ld) ...", i_inter);

       
        for(i_yama=0;i_yama<=34;i_yama++)
          {

if(strcmp(do_anisotropy,"no") == 0 && i_yama==1){break;}

if(strcmp(type_of_survey,"periodicFKP") == 0 && i_yama==1){break;}

if(i_yama==1 && delta2_sw==0){i_yama=7;}

if(i_yama==7 && delta4_sw==0){i_yama=22;}

if(i_yama==22 && delta1_sw==0){i_yama=25;}

if(i_yama==25 && delta3_sw==0){i_yama=35;}

if(i_yama==35){break;}



        delta_data = (double*) calloc(ngridtot, sizeof(double));
        for(c=0; c<Ndata; c++)
        {

if(i_yama==0){weight_pos=1;}
if(i_yama==1){weight_pos=pos_x[c]*pos_x[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==2){weight_pos=pos_x[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==3){weight_pos=pos_x[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==4){weight_pos=pos_y[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==5){weight_pos=pos_y[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}

if(i_yama==6){weight_pos=pos_z[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==7){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==8){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==9){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==10){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==11){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==12){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==13){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==14){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==15){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==16){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==17){weight_pos=pow(pos_x[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==18){weight_pos=pow(pos_y[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==19){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==20){weight_pos=pow(pos_y[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==21){weight_pos=pow(pos_z[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}

if(i_yama==22){weight_pos=pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==23){weight_pos=pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==24){weight_pos=pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],0.5);}
if(i_yama==25){weight_pos=pos_x[c]*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==26){weight_pos=pos_y[c]*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==27){weight_pos=pos_z[c]*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==28){weight_pos=pos_x[c]*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==29){weight_pos=pos_x[c]*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==30){weight_pos=pos_y[c]*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==31){weight_pos=pos_y[c]*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==32){weight_pos=pos_z[c]*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==33){weight_pos=pos_z[c]*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}
if(i_yama==34){weight_pos=pos_x[c]*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],1.5);}


if(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]==0){weight_pos=0;}
      

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
    }

        delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {

if(i_yama==0){weight_pos=1;}
if(i_yama==1){weight_pos=pos_x_rand[c]*pos_x_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==2){weight_pos=pos_x_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==3){weight_pos=pos_x_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==4){weight_pos=pos_y_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==5){weight_pos=pos_y_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==6){weight_pos=pos_z_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==7){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==8){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==9){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==10){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==11){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==12){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==13){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==14){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==15){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==16){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==17){weight_pos=pow(pos_x_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==18){weight_pos=pow(pos_y_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==19){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==20){weight_pos=pow(pos_y_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==21){weight_pos=pow(pos_z_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==22){weight_pos=pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==23){weight_pos=pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==24){weight_pos=pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],0.5);}
if(i_yama==25){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==26){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==27){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==28){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==29){weight_pos=pos_x_rand[c]*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==30){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==31){weight_pos=pos_y_rand[c]*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==32){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==33){weight_pos=pos_z_rand[c]*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}
if(i_yama==34){weight_pos=pos_x_rand[c]*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],1.5);}

if(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]==0){weight_pos=0;}

              if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
        }


        #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);

if(i_yama==0){
        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==1)
{
        deltak_re2_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==7)
{
        deltak_re4_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==22)
{
        deltak_re1_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im1_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==25)
{
        deltak_re3_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im3_b= (double*) calloc(ngridtotr2c,sizeof(double));
}



if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2_b, deltak_im2_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=7 && i_yama<=21){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4_b, deltak_im4_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=22 && i_yama<=24){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re1_b, deltak_im1_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=25 && i_yama<=34){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re3_b, deltak_im3_b, ngrid,L1b,L2b, mode_correction);}

free(delta_data);
              
  //B
              
if(strcmp(type_of_code, "rusticoX") == 0){
                      delta_dataB = (double*) calloc(ngridtot, sizeof(double));
                      for(c=0; c<NdataB; c++)
                      {

              if(i_yama==0){weight_pos=1;}
              if(i_yama==1){weight_pos=pos_xB[c]*pos_xB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
              if(i_yama==2){weight_pos=pos_xB[c]*pos_yB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
              if(i_yama==3){weight_pos=pos_xB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
              if(i_yama==4){weight_pos=pos_yB[c]*pos_yB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
              if(i_yama==5){weight_pos=pos_yB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}

              if(i_yama==6){weight_pos=pos_zB[c]*pos_zB[c]/(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]);}
              if(i_yama==7){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==8){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==9){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==10){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==11){weight_pos=pow(pos_xB[c],2)*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==12){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==13){weight_pos=pow(pos_yB[c],2)*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==14){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==15){weight_pos=pow(pos_zB[c],2)*pos_zB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==16){weight_pos=pow(pos_xB[c],2)*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==17){weight_pos=pow(pos_xB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==18){weight_pos=pow(pos_yB[c],2)*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==19){weight_pos=pow(pos_xB[c],2)*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==20){weight_pos=pow(pos_yB[c],2)*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}
              if(i_yama==21){weight_pos=pow(pos_zB[c],2)*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],2);}

              if(i_yama==22){weight_pos=pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
              if(i_yama==23){weight_pos=pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
              if(i_yama==24){weight_pos=pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],0.5);}
              if(i_yama==25){weight_pos=pos_xB[c]*pos_xB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==26){weight_pos=pos_yB[c]*pos_yB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==27){weight_pos=pos_zB[c]*pos_zB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==28){weight_pos=pos_xB[c]*pos_xB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==29){weight_pos=pos_xB[c]*pos_xB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==30){weight_pos=pos_yB[c]*pos_yB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==31){weight_pos=pos_yB[c]*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==32){weight_pos=pos_zB[c]*pos_zB[c]*pos_xB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==33){weight_pos=pos_zB[c]*pos_zB[c]*pos_yB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}
              if(i_yama==34){weight_pos=pos_xB[c]*pos_yB[c]*pos_zB[c]/pow(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c],1.5);}


              if(pos_xB[c]*pos_xB[c]+pos_yB[c]*pos_yB[c]+pos_zB[c]*pos_zB[c]==0){weight_pos=0;}
                    

                     if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                     if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                     if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                     if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                     if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                     if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c]*weight_pos, L2b, L1b, ngrid);}
                  }

                      delta_randB = (double*) calloc(ngridtot, sizeof(double));
                     for(c=0; c<NrandB; c++)
                      {

              if(i_yama==0){weight_pos=1;}
              if(i_yama==1){weight_pos=pos_x_randB[c]*pos_x_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==2){weight_pos=pos_x_randB[c]*pos_y_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==3){weight_pos=pos_x_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==4){weight_pos=pos_y_randB[c]*pos_y_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==5){weight_pos=pos_y_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==6){weight_pos=pos_z_randB[c]*pos_z_randB[c]/(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]);}
              if(i_yama==7){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==8){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==9){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==10){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==11){weight_pos=pow(pos_x_randB[c],2)*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==12){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==13){weight_pos=pow(pos_y_randB[c],2)*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==14){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==15){weight_pos=pow(pos_z_randB[c],2)*pos_z_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==16){weight_pos=pow(pos_x_randB[c],2)*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==17){weight_pos=pow(pos_x_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==18){weight_pos=pow(pos_y_randB[c],2)*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==19){weight_pos=pow(pos_x_randB[c],2)*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==20){weight_pos=pow(pos_y_randB[c],2)*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==21){weight_pos=pow(pos_z_randB[c],2)*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],2);}
              if(i_yama==22){weight_pos=pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
              if(i_yama==23){weight_pos=pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
              if(i_yama==24){weight_pos=pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],0.5);}
              if(i_yama==25){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==26){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==27){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==28){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==29){weight_pos=pos_x_randB[c]*pos_x_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==30){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==31){weight_pos=pos_y_randB[c]*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==32){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_x_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==33){weight_pos=pos_z_randB[c]*pos_z_randB[c]*pos_y_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}
              if(i_yama==34){weight_pos=pos_x_randB[c]*pos_y_randB[c]*pos_z_randB[c]/pow(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c],1.5);}

              if(pos_x_randB[c]*pos_x_randB[c]+pos_y_randB[c]*pos_y_randB[c]+pos_z_randB[c]*pos_z_randB[c]==0){weight_pos=0;}

                            if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                            if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                            if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                            if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                            if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                            if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_randB, pos_x_randB[c], pos_y_randB[c], pos_z_randB[c], weight_randB[c]*weight_pos, L2b, L1b, ngrid);}
                      }


                      #pragma omp parallel for private(c) shared(ngrid,ngridtot,alphaB,delta_dataB,delta_randB)
                      for(c=0;c<ngridtot;c++)
                      {
                                               delta_dataB[c]+=alphaB*delta_randB[c];
                      }
                      free(delta_randB);

              if(i_yama==0){
                      deltak_re0_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                      deltak_im0_bB= (double*) calloc(ngridtotr2c,sizeof(double));
              }
              if(i_yama==1)
              {
                      deltak_re2_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                      deltak_im2_bB= (double*) calloc(ngridtotr2c,sizeof(double));
              }
              if(i_yama==7)
              {
                      deltak_re4_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                      deltak_im4_bB= (double*) calloc(ngridtotr2c,sizeof(double));
              }
              if(i_yama==22)
              {
                      deltak_re1_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                      deltak_im1_bB= (double*) calloc(ngridtotr2c,sizeof(double));
              }
              if(i_yama==25)
              {
                      deltak_re3_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                      deltak_im3_bB= (double*) calloc(ngridtotr2c,sizeof(double));
              }



              if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re0_bB, deltak_im0_bB, ngrid,L1b,L2b, mode_correction);}
              if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re2_bB, deltak_im2_bB, ngrid,L1b,L2b, mode_correction);}
              if(i_yama>=7 && i_yama<=21){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re4_bB, deltak_im4_bB, ngrid,L1b,L2b, mode_correction);}
              if(i_yama>=22 && i_yama<=24){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re1_bB, deltak_im1_bB, ngrid,L1b,L2b, mode_correction);}
              if(i_yama>=25 && i_yama<=34){fftw_yamamoto_skycut(i_yama, delta_dataB, deltak_re3_bB, deltak_im3_bB, ngrid,L1b,L2b, mode_correction);}

              free(delta_dataB);
          }
              
}


#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b,new_deltak_re1_b,new_deltak_im1_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re3_b,new_deltak_im3_b,new_deltak_re4_b,new_deltak_im4_b,new_deltak_re0_bB,new_deltak_im0_bB,new_deltak_re1_bB,new_deltak_im1_bB,new_deltak_re2_bB,new_deltak_im2_bB,new_deltak_re3_bB,new_deltak_im3_bB,new_deltak_re4_bB,new_deltak_im4_bB) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re1,deltak_im1,deltak_re2,deltak_im2,deltak_re3,deltak_im3,deltak_re4,deltak_im4,deltak_re0_b,deltak_im0_b,deltak_re1_b,deltak_im1_b,deltak_re2_b,deltak_im2_b,deltak_re3_b,deltak_im3_b,deltak_re4_b,deltak_im4_b,deltak_re0B,deltak_im0B,deltak_re1B,deltak_im1B,deltak_re2B,deltak_im2B,deltak_re3B,deltak_im3B,deltak_re4B,deltak_im4B,deltak_re0_bB,deltak_im0_bB,deltak_re1_bB,deltak_im1_bB,deltak_re2_bB,deltak_im2_bB,deltak_re3_bB,deltak_im3_bB,deltak_re4_bB,deltak_im4_bB,kx,L2,L1,i_inter,Ninterlacing,delta1_sw,delta2_sw,delta3_sw,delta4_sw,do_anisotropy,type_of_code,type_of_survey)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(long int)(index2/(ngrid*ngrid/2+ngrid));
j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);


phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;
    if(strcmp(type_of_code, "rusticoX") == 0){
        new_deltak_re0_bB=deltak_re0_bB[index2]*phase_cos-deltak_im0_bB[index2]*phase_sin;
        new_deltak_im0_bB=deltak_re0_bB[index2]*phase_sin+deltak_im0_bB[index2]*phase_cos;

        deltak_re0B[index2]+=new_deltak_re0_bB;
        deltak_im0B[index2]+=new_deltak_im0_bB;

    }

if(strcmp(do_anisotropy,"yes") ==0 ){

if(delta1_sw==1){

new_deltak_re1_b=deltak_re1_b[index2]*phase_cos-deltak_im1_b[index2]*phase_sin;
new_deltak_im1_b=deltak_re1_b[index2]*phase_sin+deltak_im1_b[index2]*phase_cos;

deltak_re1[index2]+=new_deltak_re1_b;
deltak_im1[index2]+=new_deltak_im1_b;
    if(strcmp(type_of_code, "rusticoX") == 0){
        new_deltak_re1_bB=deltak_re1_bB[index2]*phase_cos-deltak_im1_bB[index2]*phase_sin;
        new_deltak_im1_bB=deltak_re1_bB[index2]*phase_sin+deltak_im1_bB[index2]*phase_cos;

        deltak_re1B[index2]+=new_deltak_re1_bB;
        deltak_im1B[index2]+=new_deltak_im1_bB;

    }
}

if(delta2_sw==1){

new_deltak_re2_b=deltak_re2_b[index2]*phase_cos-deltak_im2_b[index2]*phase_sin;
new_deltak_im2_b=deltak_re2_b[index2]*phase_sin+deltak_im2_b[index2]*phase_cos;

deltak_re2[index2]+=new_deltak_re2_b;
deltak_im2[index2]+=new_deltak_im2_b;
    if(strcmp(type_of_code, "rusticoX") == 0){

        new_deltak_re2_bB=deltak_re2_bB[index2]*phase_cos-deltak_im2_bB[index2]*phase_sin;
        new_deltak_im2_bB=deltak_re2_bB[index2]*phase_sin+deltak_im2_bB[index2]*phase_cos;

        deltak_re2B[index2]+=new_deltak_re2_bB;
        deltak_im2B[index2]+=new_deltak_im2_bB;
    }

}

if(delta3_sw==1){

new_deltak_re3_b=deltak_re3_b[index2]*phase_cos-deltak_im3_b[index2]*phase_sin;
new_deltak_im3_b=deltak_re3_b[index2]*phase_sin+deltak_im3_b[index2]*phase_cos;

deltak_re3[index2]+=new_deltak_re3_b;
deltak_im3[index2]+=new_deltak_im3_b;
    if(strcmp(type_of_code, "rusticoX") == 0){
        new_deltak_re3_bB=deltak_re3_bB[index2]*phase_cos-deltak_im3_bB[index2]*phase_sin;
        new_deltak_im3_bB=deltak_re3_bB[index2]*phase_sin+deltak_im3_bB[index2]*phase_cos;

        deltak_re3B[index2]+=new_deltak_re3_bB;
        deltak_im3B[index2]+=new_deltak_im3_bB;

    }
}

if(delta4_sw==1){

new_deltak_re4_b=deltak_re4_b[index2]*phase_cos-deltak_im4_b[index2]*phase_sin;
new_deltak_im4_b=deltak_re4_b[index2]*phase_sin+deltak_im4_b[index2]*phase_cos;

deltak_re4[index2]+=new_deltak_re4_b;
deltak_im4[index2]+=new_deltak_im4_b;
    if(strcmp(type_of_code, "rusticoX") == 0){
        new_deltak_re4_bB=deltak_re4_bB[index2]*phase_cos-deltak_im4_bB[index2]*phase_sin;
        new_deltak_im4_bB=deltak_re4_bB[index2]*phase_sin+deltak_im4_bB[index2]*phase_cos;

        deltak_re4B[index2]+=new_deltak_re4_bB;
        deltak_im4B[index2]+=new_deltak_im4_bB;
    }
}

}

}

        free(deltak_re0_b);
        free(deltak_im0_b);
      if(strcmp(type_of_code, "rusticoX") == 0){
          free(deltak_re0_bB);
          free(deltak_im0_bB);
      }

if(strcmp(do_anisotropy,"yes") ==0 ){
if(delta1_sw==1){
        free(deltak_re1_b);
        free(deltak_im1_b);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re1_bB);
        free(deltak_im1_bB);
    }
}
if(delta2_sw==1){
        free(deltak_re2_b);
        free(deltak_im2_b);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re2_bB);
        free(deltak_im2_bB);
    }
}
if(delta3_sw==1){
        free(deltak_re3_b);
        free(deltak_im3_b);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re3_bB);
        free(deltak_im3_bB);
    }
}
if(delta4_sw==1){
        free(deltak_re4_b);
        free(deltak_im4_b);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re4_bB);
        free(deltak_im4_bB);
    }
}
}

printf("Ok!\n");

  }


if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0 )
{
//printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);
free(pos_x_rand);
free(pos_y_rand);
free(pos_z_rand);
free(weight_rand);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(pos_xB);
        free(pos_yB);
        free(pos_zB);
        free(weightB);
        free(pos_x_randB);
        free(pos_y_randB);
        free(pos_z_randB);
        free(weight_randB);
    }

}

}
free(kx);

if( strcmp( output_density,"corrected") == 0 )//change write_PS do not erase them
{
d=pow((L2-L1),-3);
generate_deltax_field(deltak_re0, deltak_im0, ngrid, Ninterlacing, name_density,d);

if(strcmp(type_of_code, "rusticoX") == 0){generate_deltax_field(deltak_re0B, deltak_im0B, ngrid, Ninterlacing, name_densityB,d);}

}

if(strcmp(type_of_code, "rustico") == 0){
    printf("Writing Power Spectrum output %s...",name_ps_out);
    
}
if(strcmp(type_of_code, "rusticoX") == 0)
{
    printf("Writing Power Spectrum output %s %s %s...",name_ps_out,name_psAB_out,name_psBB_out);

}
sprintf(type,"FFT");

if( strcmp(type_of_survey,"cutsky") == 0){
    write_power_spectrum_skyscuts(kmin,kmax,NULL, NULL, NULL, deltak_re0, deltak_im0,deltak_re1, deltak_im1, deltak_re2, deltak_im2, deltak_re3, deltak_im3, deltak_re4, deltak_im4,deltak_re0B, deltak_im0B,deltak_re1B, deltak_im1B, deltak_re2B, deltak_im2B, deltak_re3B, deltak_im3B, deltak_re4B, deltak_im4B, bin_ps, ngrid, 0, L1, L2, I22,I22B, name_ps_out,name_psAB_out,name_psBB_out, P_shot_noise,P_shot_noiseB, binning_type, do_anisotropy, do_odd_multipoles,  Quadrupole_type, Octopole_type, Hexadecapole_type, type ,Ninterlacing,type_of_code, write_kvectors, name_ps_kvectors);
}

if( strcmp(type_of_survey,"periodicFKP") == 0){

if(Nmu == 1){

}
if(Nmu > 1){


}

}


//free deltas
free(deltak_re0);
free(deltak_im0);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re0B);
        free(deltak_im0B);
    }
if(strcmp(do_anisotropy,"yes") ==0 ){

if(delta2_sw==1){
          free(deltak_re2);
          free(deltak_im2);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re2B);
        free(deltak_im2B);
    }
    
    
}

if(delta1_sw==1){
          free(deltak_re1);
          free(deltak_im1);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re1B);
        free(deltak_im1B);
    }
    
    
}

if(delta3_sw==1){
          free(deltak_re3);
          free(deltak_im3);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re3B);
        free(deltak_im3B);
    }
    
    
}

if(delta4_sw==1){

          free(deltak_re4);
          free(deltak_im4);
    if(strcmp(type_of_code, "rusticoX") == 0){
        free(deltak_re4B);
        free(deltak_im4B);
    }
    
    
}
}


printf("Ok!\n");

}


void loop_interlacing_periodic_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata,double Ndataw, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im, char *name_density, char *type_of_input)
{
  FILE *f;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  double d;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));

    if( strcmp(type_of_input,"particles") == 0){

     printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);
     for(c=0; c<Ndata; c++)
     {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}

     }

        #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,Ndataw,delta_data)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndataw-pow(ngrid,-3);
        }



}//particles
else{//density
     printf("Reading densities ...");
f=fopen(name_density,"r");if(f==NULL){printf("Error reading density file %s. Exiting now...\n",name_density);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&d);
        delta_data[c]=d/ngridtot*1.;
     }
fclose(f);
}

     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
     {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(long int)(index2/(ngrid*ngrid/2+ngrid));
            j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;
}

        free(deltak_re_b);
        free(deltak_im_b);

        printf("Ok!\n");


}

}//end of itteration loop
free(kx);

}

void loop_interlacing_periodic(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double Ndataw,double *pos_xB, double *pos_yB, double *pos_zB, double *weightB, long int NdataB, double NdataBw, double L1, double L2, int ngrid, double P_shot_noise,double P_shot_noiseB, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out,char *name_psAB_out,char *name_psBB_out, char *type_of_mass_assigment, char *do_odd_multipoles,char *do_anisotropy, char *do_bispectrum,int Nmu, char *file_for_mu, char *type_of_code,char *output_density, char  *name_density, char *name_densityB,char *type_of_input, char *write_kvectors, char *name_ps_kvectors)
{
  FILE *f;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin; 
  double* deltak_re;
  double* deltak_im;
  
  double* deltak_reB;
  double* deltak_imB;

  double* deltak_re_b;
  double* deltak_im_b;

  double* deltak_re_bB;
  double* deltak_im_bB;
    
  double new_deltak_re_b,new_deltak_im_b;
  double new_deltak_re_bB,new_deltak_im_bB;
    
  double* delta_data;
  double* delta_dataB;

  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  double d;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{ 
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));//inizialized to 0.
//     if(delta_data==NULL){printf("Warning, grid-array (delta_data) couldn't be created. Exiting now...\n");exit(0);}
     //delta_data = (double*) malloc(sizeof(double)*ngridtot);memset(delta_data,0,sizeof(double));

    if( strcmp(type_of_input,"particles") == 0){
//printf("%lf, %lf\n",L1b,L2b);
     printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);
     for(c=0; c<Ndata; c++)
     {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}

     }

//exit(0); 
//long int dumb=0;
        #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,Ndataw,delta_data)
        for(c=0;c<ngridtot;c++)
        {
        //         if(delta_data[c]==0){dumb++;}
                                 delta_data[c]=delta_data[c]/Ndataw-pow(ngrid,-3);//print out delta_data[c]*ngrid^3 which is delta(x) 
                                 //modify delta inputed by delta/ngrid3 to match the definition from here (norm of FFTW)        
        }
//printf("%ld / %ld\n",dumb,ngridtot);
//exit(0);

if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_density,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_data[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}

}//particles
else{//density
     printf("Reading densities ...");
f=fopen(name_density,"r");if(f==NULL){printf("Error reading density file %s. Exiting now...\n",name_density);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&d);
        delta_data[c]=d/ngridtot*1.;
     }
fclose(f);
}        

    if(strcmp(type_of_code, "rusticoX") == 0){
   
        delta_dataB = (double*) calloc(ngridtot, sizeof(double));

    if( strcmp(type_of_input,"particles") == 0){

//        if(delta_dataB==NULL){printf("Warning, grid-array (delta_dataB) couldn't be created. Exiting now...\n");exit(0);}
        printf("Assigning particles-B to the grid (Iteration %ld) ...", i_inter);
        for(c=0; c<NdataB; c++)
        {
          if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
          if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
          if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
          if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
          if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
          if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}

        }


           #pragma omp parallel for private(c) shared(ngrid,ngridtot,NdataB,NdataBw,delta_dataB)
           for(c=0;c<ngridtot;c++)
           {
                                    delta_dataB[c]=delta_dataB[c]/NdataBw-pow(ngrid,-3);
           }
if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_densityB,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_dataB[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}

}//particles
else{//density
     printf("Reading densities ...");
     f=fopen(name_densityB,"r");if(f==NULL){printf("Error reading density file %s. Exiting now...\n",name_densityB);exit(0);}
     for(c=0;c<ngridtot;c++)
     {
        fscanf(f,"%lf\n",&d);
        delta_dataB[c]=d/ngridtot*1.;
     }
fclose(f);
}

        
    }

     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        deltak_re= (double*) calloc(ngridtotr2c,sizeof(double));
//        if(deltak_re==NULL){printf("Warning, grid-array (deltak_re) couldn't be created. Exiting now...\n");exit(0);}
        deltak_im= (double*) calloc(ngridtotr2c,sizeof(double));
//        if(deltak_im==NULL){printf("Warning, grid-array (deltak_im) couldn't be created. Exiting now...\n");exit(0);}
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
 
         if(strcmp(type_of_code, "rusticoX") == 0){

             deltak_reB= (double*) calloc(ngridtotr2c,sizeof(double));
             deltak_imB= (double*) calloc(ngridtotr2c,sizeof(double));
             fftw_yamamoto_skycut(0,delta_dataB, deltak_reB, deltak_imB, ngrid,L1b,L2b, mode_correction);
             free(delta_dataB);
             
         }
         
        printf("Ok!\n");

     }
     else
     {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(long int)(index2/(ngrid*ngrid/2+ngrid));
            j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;

}

        free(deltak_re_b);
        free(deltak_im_b);
         
         if(strcmp(type_of_code, "rusticoX") == 0){
             
                     deltak_re_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                     deltak_im_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                     fftw_yamamoto_skycut(0,delta_dataB, deltak_re_bB, deltak_im_bB, ngrid,L1b,L2b, mode_correction);
                     free(delta_dataB);
                   #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_bB,new_deltak_im_bB) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_reB,deltak_imB,deltak_re_bB,deltak_im_bB,i_inter,Ninterlacing,L1,L2)
                   for(index2=0;index2<ngridtotr2c;index2++)
                   {

                         i=(long int)(index2/(ngrid*ngrid/2+ngrid));
                         j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
                         k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

                          phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                          phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                          new_deltak_re_bB=deltak_re_bB[index2]*phase_cos-deltak_im_bB[index2]*phase_sin;
                          new_deltak_im_bB=deltak_re_bB[index2]*phase_sin+deltak_im_bB[index2]*phase_cos;
                          deltak_reB[index2]+=new_deltak_re_bB;
                          deltak_imB[index2]+=new_deltak_im_bB;

             }

                     free(deltak_re_bB);
                     free(deltak_im_bB);
             
         }
        
        printf("Ok!\n");

if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0 )//only free them in the last interation and if no-bispectrum is computed
{

    if( strcmp(type_of_input,"particles") == 0){
printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);}
    
    if(strcmp(type_of_code, "rusticoX") == 0){
        
    if( strcmp(type_of_input,"particles") == 0){
        free(pos_xB);
        free(pos_yB);
        free(pos_zB);
        free(weightB);}
    }
    
}


}

}//end of itteration loop

free(kx);


if( strcmp( output_density,"corrected") == 0 )//change write_PS do not erase them
{
generate_deltax_field(deltak_re, deltak_im, ngrid, Ninterlacing, name_density, 1.0);

if(strcmp(type_of_code, "rusticoX") == 0){generate_deltax_field(deltak_reB, deltak_imB, ngrid, Ninterlacing, name_densityB, 1.0);}

}


    if(strcmp(type_of_code, "rustico") == 0){printf("Writing Power Spectrum output %s...",name_ps_out);}
    if(strcmp(type_of_code, "rusticoX") == 0){printf("Writing Power Spectrum output %s %s %s...",name_ps_out,name_psAB_out,name_psBB_out);}



if(Nmu == 1){

    write_power_spectrum_periodic(kmin,kmax,deltak_re,deltak_im,deltak_reB,deltak_imB,bin_ps, ngrid,L1,L2,Ninterlacing,name_ps_out,name_psAB_out,name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,do_odd_multipoles,do_anisotropy,type_of_code,0,0,write_kvectors, name_ps_kvectors);

}
else{
    write_power_spectrum_periodic2D(kmin,kmax,deltak_re,deltak_im,deltak_reB,deltak_imB,bin_ps,Nmu, ngrid,L1,L2,Ninterlacing,name_ps_out,name_psAB_out,name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,file_for_mu,type_of_code,0,0,write_kvectors, name_ps_kvectors);

}


printf("Ok!\n");

}


void loop_interlacing_periodic_gadget(double kmin,double kmax, int Ninterlacing, char *name_data_in,char *name_dataB_in ,int gadget_files,int gadget_filesB, double L1, double L2, int ngrid, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out,char *name_psAB_out,char *name_psBB_out, char *type_of_mass_assigment,double Shot_noise_factor,char *grid_correction_string, char *RSD,char *RSDB, char *do_odd_multipoles,char *do_anisotropy, int Nmu, char *file_for_mu,char *type_of_code, char *output_density, char *name_density,char *name_densityB,char *type_of_input,char *write_kvectors, char *name_ps_kvectors)
{
  FILE *f;
  double *pos_x,*pos_y,*pos_z;
    double *pos_xB,*pos_yB,*pos_zB;
  double weight;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re;
  double* deltak_im;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;

    double* deltak_reB;
    double* deltak_imB;
    double* deltak_re_bB;
    double* deltak_im_bB;
    double new_deltak_re_bB,new_deltak_im_bB;
    double* delta_dataB;
    
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  int snapshot_num;
  char name_gadget_file[500];
  long int NumPart_file;
  long int Ndata;
  double P_shot_noise;
  double scale_factor;
  double Omatter,Olambda;
    
    char name_gadget_fileB[500];
    long int NumPart_fileB;
    long int NdataB;
    double P_shot_noiseB;
    double scale_factorB;
    double OmatterB,OlambdaB;
    double d;
    
  double params[4];
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

weight=1.0;
    for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
    {
         L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
         L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

         delta_data = (double*) calloc(ngridtot, sizeof(double));
         printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);

    Ndata=0;
    for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
      {

      sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);

      load_snapshot(name_gadget_file, 1, params);
    NumPart_file=(long int)(params[0]);
    scale_factor=1./(1+params[1]);
    Omatter=params[2];
    Olambda=params[3];

    Ndata=Ndata+NumPart_file;
     pos_x = (double*) calloc(NumPart_file, sizeof(double));
     pos_y = (double*) calloc(NumPart_file, sizeof(double));
     pos_z = (double*) calloc(NumPart_file, sizeof(double));

         for(c=0; c<NumPart_file; c++)
         {

         pos_x[c]=P[c].Pos[0]*0.001;
         pos_y[c]=P[c].Pos[1]*0.001;
    if(strcmp(RSD, "no") == 0){pos_z[c]=P[c].Pos[2]*0.001;}
    if(strcmp(RSD, "yes") == 0){pos_z[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factor)/(100.*scale_factor*sqrt(Omatter*pow(scale_factor,-3)+Olambda));}

//fprintf(f,"%e %e %e\n",pos_x[c],pos_y[c],pos_z[c]);

           if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
           if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
           if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
           if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
           if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
           if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}

         }
    free(pos_x);
    free(pos_y);
    free(pos_z);

    printf("Ok!\n");
    }

            #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
         for(c=0;c<ngridtot;c++)
            {
                                     delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
            }
if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_density,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_data[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}




         printf("Ok!\n");
         printf("Performing FFTs ...");
       if(i_inter==1)
         {
            deltak_re= (double*) calloc(ngridtotr2c,sizeof(double));
            deltak_im= (double*) calloc(ngridtotr2c,sizeof(double));
            fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
            free(delta_data);
            printf("Ok!\n");

         }
         else
         {
            deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
            deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
            fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
            free(delta_data);
          #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
        for(index2=0;index2<ngridtotr2c;index2++)
          {

                i=(long int)(index2/(ngrid*ngrid/2+ngrid));
                j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
                k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

                 phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                 phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                 new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
                 new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
                 deltak_re[index2]+=new_deltak_re_b;
                 deltak_im[index2]+=new_deltak_im_b;
    }

            free(deltak_re_b);
            free(deltak_im_b);

            printf("Ok!\n");


    }
        
        if(strcmp(type_of_code, "rusticoX") == 0)
        {

                 delta_dataB = (double*) calloc(ngridtot, sizeof(double));
                 printf("Assigning particles-B to the grid (Iteration %ld) ...", i_inter);

            NdataB=0;
            for(snapshot_num=0;snapshot_num<gadget_filesB;snapshot_num++)
              {

              sprintf(name_gadget_fileB, "%s.%d", name_dataB_in, snapshot_num);

              load_snapshot(name_gadget_fileB, 1, params);
            NumPart_fileB=(long int)(params[0]);
            scale_factorB=1./(1+params[1]);
            OmatterB=params[2];
            OlambdaB=params[3];

            NdataB=NdataB+NumPart_file;
             pos_xB = (double*) calloc(NumPart_fileB, sizeof(double));
             pos_yB = (double*) calloc(NumPart_fileB, sizeof(double));
             pos_zB = (double*) calloc(NumPart_fileB, sizeof(double));

                 for(c=0; c<NumPart_fileB; c++)
                 {

                 pos_xB[c]=P[c].Pos[0]*0.001;
                 pos_yB[c]=P[c].Pos[1]*0.001;
            if(strcmp(RSDB, "no") == 0){pos_zB[c]=P[c].Pos[2]*0.001;}
            if(strcmp(RSDB, "yes") == 0){pos_zB[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factorB)/(100.*scale_factorB*sqrt(OmatterB*pow(scale_factorB,-3)+OlambdaB));}


                   if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}
                   if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}
                   if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}
                   if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}
                   if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}
                   if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weight, L2b, L1b, ngrid);}

                 }
            free(pos_xB);
            free(pos_yB);
            free(pos_zB);

            printf("Ok!\n");
            }

                    #pragma omp parallel for private(c) shared(ngrid,ngridtot,NdataB,delta_dataB)
                 for(c=0;c<ngridtot;c++)
                    {
                                             delta_dataB[c]=delta_dataB[c]/NdataB-pow(ngrid,-3);
                    }
if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_densityB,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_dataB[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}



                 printf("Ok!\n");
                 printf("Performing FFTs ...");
               if(i_inter==1)
                 {
                    deltak_reB= (double*) calloc(ngridtotr2c,sizeof(double));
                    deltak_imB= (double*) calloc(ngridtotr2c,sizeof(double));
                    fftw_yamamoto_skycut(0,delta_dataB, deltak_reB, deltak_imB, ngrid,L1b,L2b, mode_correction);
                    free(delta_dataB);
                    printf("Ok!\n");

                 }
                 else
                 {
                    deltak_re_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                    deltak_im_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                    fftw_yamamoto_skycut(0,delta_dataB, deltak_re_bB, deltak_im_bB, ngrid,L1b,L2b, mode_correction);
                    free(delta_dataB);
                  #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_bB,new_deltak_im_bB) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_reB,deltak_imB,deltak_re_bB,deltak_im_bB,i_inter,Ninterlacing,L1,L2)
                for(index2=0;index2<ngridtotr2c;index2++)
                  {

                        i=(long int)(index2/(ngrid*ngrid/2+ngrid));
                        j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
                        k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

                         phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                         phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                         new_deltak_re_bB=deltak_re_bB[index2]*phase_cos-deltak_im_bB[index2]*phase_sin;
                         new_deltak_im_bB=deltak_re_bB[index2]*phase_sin+deltak_im_bB[index2]*phase_cos;
                         deltak_reB[index2]+=new_deltak_re_bB;
                         deltak_imB[index2]+=new_deltak_im_bB;
            }

                    free(deltak_re_bB);
                    free(deltak_im_bB);

                    printf("Ok!\n");


            }
            
        }

    }//end of itteration loop
    free(kx);

if( strcmp( output_density,"corrected") == 0 )//change write_PS do not erase them
{
generate_deltax_field(deltak_re, deltak_im, ngrid, Ninterlacing, name_density,1.0);

if(strcmp(type_of_code, "rusticoX") == 0){generate_deltax_field(deltak_reB, deltak_imB, ngrid, Ninterlacing, name_densityB,1.0);}

}

    if(strcmp(type_of_code, "rustico") == 0){
printf("Writing Power Spectrum output %s...",name_ps_out);
P_shot_noise=pow(L2-L1,3)/Ndata;
    }
    if(strcmp(type_of_code, "rusticoX") == 0){

    printf("Writing Power Spectrum output %s %s %s...",name_ps_out,name_psAB_out,name_psBB_out);
    P_shot_noise=pow(L2-L1,3)/Ndata;
    P_shot_noiseB=pow(L2-L1,3)/NdataB;
    }
    

    if(Nmu == 1){
        write_power_spectrum_periodic(kmin,kmax,deltak_re,deltak_im,deltak_reB,deltak_imB,bin_ps, ngrid,L1,L2,Ninterlacing,name_ps_out,name_psAB_out,name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,do_odd_multipoles,do_anisotropy,type_of_code,0,0,write_kvectors, name_ps_kvectors);
        
    }
    else{
        write_power_spectrum_periodic2D(kmin,kmax,deltak_re,deltak_im,deltak_reB,deltak_imB,bin_ps,Nmu, ngrid,L1,L2,Ninterlacing,name_ps_out,name_psAB_out,name_psBB_out,P_shot_noise,P_shot_noiseB,binning_type,file_for_mu,type_of_code,0,0,write_kvectors, name_ps_kvectors);
        
    }


printf("Ok!\n");
//fclose(f);
}


void loop_interlacing_periodic_gadget_for_bispectrum(int Ninterlacing, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im, char *name_data_in,int gadget_files, double *params_input, char *RSD)
{
  double *pos_x,*pos_y,*pos_z;
  long int Ndata;
  double weight=1.0;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  int snapshot_num;
  char name_gadget_file[500];
  int NumPart_file;
  double P_shot_noise;
  double scale_factor,Omatter,Olambda;
  double params[4];

       kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));
     printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);

Ndata=0;
for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {

  sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);

  load_snapshot(name_gadget_file, 1, params);
NumPart_file=(long int)(params[0]);
scale_factor=1./(1+params[1]);
Omatter=params[2];
Olambda=params[3];

Ndata=Ndata+NumPart_file;
 pos_x = (double*) calloc(NumPart_file, sizeof(double));
 pos_y = (double*) calloc(NumPart_file, sizeof(double));
 pos_z = (double*) calloc(NumPart_file, sizeof(double));

//printf("%d\n",NumPart_file);
     for(c=0; c<NumPart_file; c++)
     {
     pos_x[c]=P[c].Pos[0]*0.001; //conversion to Mpc/h (originally in kpc/h)
     pos_y[c]=P[c].Pos[1]*0.001;
if(strcmp(RSD, "no") == 0){pos_z[c]=P[c].Pos[2]*0.001;}
if(strcmp(RSD, "yes") == 0){pos_z[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factor)/(100.*scale_factor*sqrt(Omatter*pow(scale_factor,-3)+Olambda));}
        

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}

     }
free(pos_x);
free(pos_y);
free(pos_z);
printf("Ok!\n");
}
       #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
        }



     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
    {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(long int)(index2/(ngrid*ngrid/2+ngrid));
            j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;
}

        free(deltak_re_b);
        free(deltak_im_b);

        printf("Ok!\n");


}

}//end of itteration loop
free(kx);

params_input[0]=pow(L2-L1,3)/Ndata;

}

    void loop_interlacing_periodic_gadget_x_ascii(double kmin,double kmax,int Ninterlacing, char *name_dataA_in, int gadget_filesA, double *pos_xB, double *pos_yB, double *pos_zB, double *weightB, long int NdataB, double NdataBw, double L1, double L2, int ngrid, double P_shot_noiseB,  double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_psAA_out,char *name_psAB_out, char *name_psBB_out, char *type_of_mass_assigment,double Shot_noise_factor, char *do_bispectrum, char *RSD,char *do_odd_multipoles,char *do_anisotropy, int Nmu, char *file_for_mu,char *type_of_code,  char *output_density, char *name_density, char *name_densityB, char *write_kvectors, char *name_ps_kvectors)
    {
     
         FILE *f;
          double *pos_xA,*pos_yA,*pos_zA;

         double weight;
          long int i,j,k,i_inter;
          long int c;
          double phase_cos,phase_sin;
          double* deltak_reA;
          double* deltak_imA;
          double* deltak_re_bA;
          double* deltak_im_bA;
          double new_deltak_re_bA,new_deltak_im_bA;
          double* delta_dataA;

          double* deltak_reB;
          double* deltak_imB;
          double* deltak_re_bB;
          double* deltak_im_bB;
          double new_deltak_re_bB,new_deltak_im_bB;
          double* delta_dataB;

          double L2b,L1b;
          double *kx;
          double Pi=(4.*atan(1.));
          long int ngridtot=pow(ngrid,3);
          long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
          long int index2;
          int snapshot_num;

          char name_gadget_fileA[500];
          long int NumPart_fileA;
          long int NdataA;
          double P_shot_noiseA;
          double scale_factorA;
          double OmatterA,OlambdaA;
          double params[4];
          double d;
                kx=malloc(sizeof(double)*(ngrid));
                for(i=0;i<ngrid;i++)
                {
                                if(i<ngrid/2+1)
                                {
                                        kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                                }
                                else
                                {
                                        kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                                }
                }
        weight=1.0;
        for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
        {
             L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
             L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

             delta_dataA = (double*) calloc(ngridtot, sizeof(double));
             printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);

        NdataA=0;
        for(snapshot_num=0;snapshot_num<gadget_filesA;snapshot_num++)
          {

          sprintf(name_gadget_fileA, "%s.%d", name_dataA_in, snapshot_num);

          load_snapshot(name_gadget_fileA, 1, params);
        NumPart_fileA=(long int)(params[0]);
        scale_factorA=1./(1+params[1]);
        OmatterA=params[2];
        OlambdaA=params[3];
        NdataA=NdataA+NumPart_fileA;
         pos_xA = (double*) calloc(NumPart_fileA, sizeof(double));
         pos_yA = (double*) calloc(NumPart_fileA, sizeof(double));
         pos_zA = (double*) calloc(NumPart_fileA, sizeof(double));

             for(c=0; c<NumPart_fileA; c++)
            {

             pos_xA[c]=P[c].Pos[0]*0.001;
             pos_yA[c]=P[c].Pos[1]*0.001;
        if(strcmp(RSD, "no") == 0){pos_zA[c]=P[c].Pos[2]*0.001;}
        if(strcmp(RSD, "yes") == 0){pos_zA[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factorA)/(100.*scale_factorA*sqrt(OmatterA*pow(scale_factorA,-3)+OlambdaA));}


               if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}
               if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}
               if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}
               if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}
               if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}
               if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataA, pos_xA[c], pos_yA[c], pos_zA[c], weight, L2b, L1b, ngrid);}

             }
        free(pos_xA);
        free(pos_yA);
        free(pos_zA);

        printf("Ok!\n");
        }

                #pragma omp parallel for private(c) shared(ngrid,ngridtot,NdataA,delta_dataA)
             for(c=0;c<ngridtot;c++)
                {
                                         delta_dataA[c]=delta_dataA[c]/NdataA-pow(ngrid,-3);
                }
if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_density,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_dataA[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}

           printf("Ok!\n");
             printf("Performing FFTs ...");
           if(i_inter==1)
             {
                deltak_reA= (double*) calloc(ngridtotr2c,sizeof(double));
                deltak_imA= (double*) calloc(ngridtotr2c,sizeof(double));
                fftw_yamamoto_skycut(0,delta_dataA, deltak_reA, deltak_imA, ngrid,L1b,L2b, mode_correction);
                free(delta_dataA);
                printf("Ok!\n");

             }
           else
             {
                deltak_re_bA= (double*) calloc(ngridtotr2c,sizeof(double));
                deltak_im_bA= (double*) calloc(ngridtotr2c,sizeof(double));
                fftw_yamamoto_skycut(0,delta_dataA, deltak_re_bA, deltak_im_bA, ngrid,L1b,L2b, mode_correction);
                free(delta_dataA);
              #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_bA,new_deltak_im_bA) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_reA,deltak_imA,deltak_re_bA,deltak_im_bA,i_inter,Ninterlacing,L1,L2)
            for(index2=0;index2<ngridtotr2c;index2++)
              {

                    i=(long int)(index2/(ngrid*ngrid/2+ngrid));
                    j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
                    k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

                     phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                     phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                     new_deltak_re_bA=deltak_re_bA[index2]*phase_cos-deltak_im_bA[index2]*phase_sin;
                     new_deltak_im_bA=deltak_re_bA[index2]*phase_sin+deltak_im_bA[index2]*phase_cos;
                     deltak_reA[index2]+=new_deltak_re_bA;
                     deltak_imA[index2]+=new_deltak_im_bA;
        }

                free(deltak_re_bA);
                free(deltak_im_bA);

                printf("Ok!\n");


        }
        }
        //B

        for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
        {
             L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
             L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
            delta_dataB = (double*) calloc(ngridtot, sizeof(double));
                  printf("Assigning particles to the grid (Iteration %ld) ...", i_inter);
                  for(c=0; c<NdataB; c++)
              {
              if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_dataB, pos_xB[c], pos_yB[c], pos_zB[c], weightB[c], L2b, L1b, ngrid);}
              }
                         #pragma omp parallel for private(c) shared(ngrid,ngridtot,NdataB,NdataBw,delta_dataB)
                         for(c=0;c<ngridtot;c++)
                         {
                           delta_dataB[c]=delta_dataB[c]/NdataBw-pow(ngrid,-3);
                     }
if( strcmp( output_density,"uncorrected") == 0 )
{
f=fopen(name_densityB,"w");
     for(c=0;c<ngridtot;c++)
     {  d = delta_dataB[c]*ngridtot;
        fprintf(f,"%e\n",d);
     }
fclose(f);
}


             printf("Ok!\n");
             printf("Performing FFTs ...");
             if(i_inter==1)
             {
              deltak_reB= (double*) calloc(ngridtotr2c,sizeof(double));
              deltak_imB= (double*) calloc(ngridtotr2c,sizeof(double));
              fftw_yamamoto_skycut(0,delta_dataB, deltak_reB, deltak_imB, ngrid,L1b,L2b, mode_correction);
              free(delta_dataB);
              printf("Ok!\n");
           }
             else
             {
                deltak_re_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                deltak_im_bB= (double*) calloc(ngridtotr2c,sizeof(double));
                fftw_yamamoto_skycut(0,delta_dataB, deltak_re_bB, deltak_im_bB, ngrid,L1b,L2b, mode_correction);
                free(delta_dataB);
            #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_bB,new_deltak_im_bB) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_reB,deltak_imB,deltak_re_bB,deltak_im_bB,i_inter,Ninterlacing,L1,L2)
              for(index2=0;index2<ngridtotr2c;index2++)
              {
                   i=(long int)(index2/(ngrid*ngrid/2+ngrid));
                   j=(long int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
                   k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

                   phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                   phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
                   new_deltak_re_bB=deltak_re_bB[index2]*phase_cos-deltak_im_bB[index2]*phase_sin;
                   new_deltak_im_bB=deltak_re_bB[index2]*phase_sin+deltak_im_bB[index2]*phase_cos;
                   deltak_reB[index2]+=new_deltak_re_bB;
                   deltak_imB[index2]+=new_deltak_im_bB;
                }


               free(deltak_re_bB);
               free(deltak_im_bB);
        }

                printf("Ok!\n");
        if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0 )//only free them in the last interation and if no-bispectrum is computed
        {
        free(pos_xB);
        free(pos_yB);
        free(pos_zB);
        free(weightB);

        }


        }//end of itteration loop

        

        free(kx);

if( strcmp( output_density,"corrected") == 0 )//change write_PS do not erase them
{

generate_deltax_field(deltak_reA, deltak_imA, ngrid, Ninterlacing, name_density,1.0);//gadget

if(strcmp(type_of_code, "rusticoX") == 0){generate_deltax_field(deltak_reB, deltak_imB, ngrid, Ninterlacing, name_densityB,1.0);}//ascii

}

        printf("Writing Power Spectrum output %s %s %s...",name_psAA_out,name_psAB_out,name_psBB_out);
        P_shot_noiseA=pow(L2-L1,3)/NdataA;

        if(Nmu == 1){
            write_power_spectrum_periodic(kmin,kmax,deltak_reA,deltak_imA,deltak_reB,deltak_imB,bin_ps, ngrid,L1,L2,Ninterlacing,name_psAA_out,name_psAB_out,name_psBB_out,P_shot_noiseA,P_shot_noiseB,binning_type,do_odd_multipoles,do_anisotropy,type_of_code,0,0,write_kvectors, name_ps_kvectors);
            
        }
        else{
            write_power_spectrum_periodic2D(kmin,kmax,deltak_reA,deltak_imA,deltak_reB,deltak_imB,bin_ps,Nmu, ngrid,L1,L2,Ninterlacing,name_psAA_out,name_psAB_out,name_psBB_out,P_shot_noiseA,P_shot_noiseB,binning_type,file_for_mu,type_of_code,0,0,write_kvectors, name_ps_kvectors);
            
        }

        printf("Ok!\n");

    }



