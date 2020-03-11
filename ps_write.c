#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "functions.h"

double Leg1(double x)
{
double f=x;

return f;
}

double Leg2(double x)
{
double f=0.5*(3.*x*x-1.);

return f;
}

double Leg3(double x)
{
double f=0.5*(5.*x*x-3.)*x;

return f;
}

double Leg4(double x)
{
double f=1./8.*(35.*x*x*x*x-30.*x*x+3.);

return f;
}

void write_power_spectrum_skyscuts(double kmin,double kmax,double kx[], double ky[], double kz[], double P0r[], double P0i[],double P1r[], double P1i[], double P2r[], double P2i[], double P3r[], double P3i[], double P4r[], double P4i[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type, char *do_anisotropy, char *do_odd_multipoles, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type,char *type,int N_interlacing)
{
FILE *f;
long int i,j,k,l22;
int l,tid;
int mode_duplicate;
double **K;
double **Mono;
double **Quadru;
double **Di;
double **Hexadeca;
double **Octo;
long int **nmodes;
double Pi=(4.*atan(1.));
int bintype_sw;
double kmineff;
double keff;
double KX,KY,KZ;
int ivec,jvec,kvec;
int i2vec,j2vec,k2vec;
double *kvector;
long int index2;
long int ngrid_eff;
int Nk;

if(strcmp(type,"FFT") == 0){
  kvector=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kvector[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kvector[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }
mode_duplicate=1;
ngrid_eff=ngrid*ngrid*ngrid;
}
if(strcmp(type,"DSY") == 0 || strcmp(type,"DSE") == 0 ){
mode_duplicate=2;
ngrid_eff=NGRID;
}

if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

if(Nk<=0){printf("Error, Nk<=0 in ps_write.c, Exiting now...\n");exit(0);}

       int nthreads;
  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));

if(strcmp(do_anisotropy,"yes") == 0){
        Quadru= (double **)calloc(Nk,sizeof(double*));
        Hexadeca= (double **)calloc(Nk,sizeof(double*));
if(strcmp(do_odd_multipoles,"yes") == 0){
        Di= (double **)calloc(Nk,sizeof(double*));
        Octo= (double **)calloc(Nk,sizeof(double*));}}


        K=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));

if(strcmp(do_anisotropy,"yes") == 0){

                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                Hexadeca[l] = (double*)calloc(nthreads,sizeof(double));

if(strcmp(do_odd_multipoles,"yes") == 0){
                Di[l] = (double*)calloc(nthreads,sizeof(double));
                Octo[l] = (double*)calloc(nthreads,sizeof(double));}}

                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
        }
//for(i=0;i<10000;i++){printf("%lf %lf %lf %lf %lf %lf\n",P0r[i],P0i[i],P2r[i],P2i[i],P4r[i],P4i[i]);}

//exit(0);
#pragma omp parallel for private(i2vec,j2vec,k2vec,ivec,jvec,kvec,KX,KY,KZ,i,tid,l,kmineff,keff,index2) shared(NGRID,ngrid,K,kmin,kx,ky,kz,Mono,Di,Quadru,Octo,Hexadeca,P0r,P1r,P2r,P3r,P4r,P0i,P1i,P2i,P3i,P4i,nmodes,Deltak,Nk,bintype_sw,type,do_anisotropy,do_odd_multipoles, Quadrupole_type,Hexadecapole_type,Octopole_type,N_interlacing,kvector,mode_duplicate,ngrid_eff)
        for(i=0;i<ngrid_eff;i++)
        {
                tid=omp_get_thread_num();

if(strcmp(type,"FFT") == 0){
                ivec=(int)(i/(ngrid*ngrid*1.));
                jvec=(int)( (i-ivec*ngrid*ngrid)/(ngrid*1.));
                kvec=i-ivec*ngrid*ngrid-jvec*ngrid;
KX=kvector[ivec];
KY=kvector[jvec];
KZ=kvector[kvec];
}
if(strcmp(type,"DSY") == 0 || strcmp(type,"DSE") == 0){
KX=kx[i];
KY=ky[i];
KZ=kz[i];//printf("\n i=%ld/%ld, %lf, %lf, %lf\n",i,ngrid_eff,KX,KY,KZ);
}

if(bintype_sw==0){

l=(int)((pow(KX*KX+KY*KY+KZ*KZ,0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(KX*KX+KY*KY+KZ*KZ,0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
if(KZ==0)
{
 
if(strcmp(type,"DSE") == 0 && KX*KX+KY*KY+KZ*KZ>0){

               K[l][tid]+=pow(KX*KX+KY*KY+KZ*KZ,0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=P0r[i];
               nmodes[l][tid]+=1;
               Quadru[l][tid]+=P2r[i];
               Hexadeca[l][tid]+=P4r[i];

if(strcmp(do_odd_multipoles,"yes") == 0){
               Di[l][tid]+=P1r[i];
               Octo[l][tid]+=P3r[i];}

}
if(strcmp(type,"DSY") == 0 || strcmp(type,"FFT") == 0){

if(strcmp(type,"FFT") == 0){index2=((pow(ngrid,2)*ivec+ngrid*jvec+2*kvec)/2+ivec*ngrid+jvec);}
if(strcmp(type,"DSY") == 0){index2=i;}

if(KX*KX+KY*KY+KZ*KZ>0){

               K[l][tid]+=pow(KX*KX+KY*KY+KZ*KZ,0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=(pow(P0r[index2],2)+pow(P0i[index2],2));
               nmodes[l][tid]+=1;//printf("\n(1)\n");

//if(l==0){printf("A %lf %lf %lf, %lf %lf\n",KX,KY,KZ,P0r[index2],P0i[index2]);}

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(Quadrupole_type,"L2")== 0){Quadru[l][tid]+=P0r[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0i[index2]*0.5*(3.*P2i[index2]-P0i[index2]);}
if(strcmp(Quadrupole_type,"L1L1")== 0){Quadru[l][tid]+=3./2.*(pow(P1r[index2],2)+pow(P1i[index2],2))-0.5*(pow(P0r[index2],2)+pow(P0i[index2],2));}

if(strcmp(Hexadecapole_type,"L4")== 0){Hexadeca[l][tid]+=P0r[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0i[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]);}
if(strcmp(Hexadecapole_type,"L2L2")== 0){Hexadeca[l][tid]+=35./8.*(pow(P2r[index2],2)+pow(P2i[index2],2))-30./8.*(P0r[index2]*P2r[index2]+P0i[index2]*P2i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2));}
if(strcmp(Hexadecapole_type,"L1L3")== 0){Hexadeca[l][tid]+=35./8.*(P3r[index2]*P1r[index2]+P3i[index2]*P1i[index2])-30./8.*(P1r[index2]*P1r[index2]+P1i[index2]*P1i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2));}
if(strcmp(do_odd_multipoles,"yes") == 0){
Di[l][tid]+=P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2];
if(strcmp(Octopole_type,"L3") == 0){Octo[l][tid]+=0.5*P0i[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0r[index2]*(5.*P3i[index2]-3.*P1i[index2]);}
if(strcmp(Octopole_type,"L1L2") == 0){Octo[l][tid]+=5./2.*(P2r[index2]*P1r[index2]+P2i[index2]*P1i[index2])-3./2.*(P0r[index2]*P1r[index2]+P0i[index2]*P1i[index2]);}
}
}
}//k>0
}               
}
else//KZ>0 for DS (as kz can't be negative on DS options
{

if(strcmp(type,"DSE") == 0 && KX*KX+KY*KY+KZ*KZ>0){

               K[l][tid]+=mode_duplicate*pow(KX*KX+KY*KY+KZ*KZ,0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=mode_duplicate*P0r[i];
               nmodes[l][tid]+=mode_duplicate;
               Quadru[l][tid]+=mode_duplicate*P2r[i];
               Hexadeca[l][tid]+=mode_duplicate*P4r[i];

if(strcmp(do_odd_multipoles,"yes") == 0){
               Di[l][tid]+=mode_duplicate*P1r[i];
               Octo[l][tid]+=mode_duplicate*P3r[i];}
}
if(strcmp(type,"DSY") == 0 || strcmp(type,"FFT") == 0){

if(strcmp(type,"DSY") == 0){index2=i;}
if(strcmp(type,"FFT") == 0 && KZ>=0){index2=((pow(ngrid,2)*ivec+ngrid*jvec+2*kvec)/2+ivec*ngrid+jvec);}
 if(strcmp(type,"FFT") == 0 && KZ<0){
 k2vec=ngrid-kvec;
 i2vec=ngrid-ivec;
 j2vec=ngrid-jvec;
 if(i2vec==ngrid){i2vec=0;}
 if(j2vec==ngrid){j2vec=0;}
 if(k2vec==ngrid){k2vec=0;}
 index2=((pow(ngrid,2)*i2vec+ngrid*j2vec+2*k2vec)/2+i2vec*ngrid+j2vec);}

if(KX*KX+KY*KY+KZ*KZ>0){

//if(l==0){printf("B %lf %lf %lf, %lf %lf\n",KX,KY,KZ,P0r[index2],P0i[index2]);}

               K[l][tid]+=mode_duplicate*pow(KX*KX+KY*KY+KZ*KZ,0.5);//two modes coming from the condition kz>0
               nmodes[l][tid]+=mode_duplicate;//printf("\n(2)\n");
               Mono[l][tid]+=mode_duplicate*(pow(P0r[index2],2)+pow(P0i[index2],2));

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(Quadrupole_type,"L2")== 0){Quadru[l][tid]+=mode_duplicate*(P0r[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0i[index2]*0.5*(3.*P2i[index2]-P0i[index2]));}
if(strcmp(Quadrupole_type,"L1L1")== 0){Quadru[l][tid]+=mode_duplicate*(3./2.*(pow(P1r[index2],2)+pow(P1i[index2],2))-0.5*(pow(P0r[index2],2)+pow(P0i[index2],2)));}

if(strcmp(Hexadecapole_type,"L4")== 0){Hexadeca[l][tid]+=mode_duplicate*(P0r[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0i[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]));}
if(strcmp(Hexadecapole_type,"L2L2")== 0){Hexadeca[l][tid]+=mode_duplicate*(35./8.*(pow(P2r[index2],2)+pow(P2i[index2],2))-30./8.*(P0r[index2]*P2r[index2]+P0i[index2]*P2i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2)));}
if(strcmp(Hexadecapole_type,"L1L3")== 0){Hexadeca[l][tid]+=mode_duplicate*(35./8.*(P3r[index2]*P1r[index2]+P3i[index2]*P1i[index2])-30./8.*(P1r[index2]*P1r[index2]+P1i[index2]*P1i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2)));}

if(strcmp(do_odd_multipoles,"yes") == 0){

Di[l][tid]+=mode_duplicate*(P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2]);
if(strcmp(Octopole_type,"L3") == 0){Octo[l][tid]+=mode_duplicate*(0.5*P0i[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0r[index2]*(5.*P3i[index2]-3.*P1i[index2]));}
if(strcmp(Octopole_type,"L1L2") == 0){Octo[l][tid]+=mode_duplicate*(5./2.*(P2r[index2]*P1r[index2]+P2i[index2]*P1i[index2])-3./2.*(P0r[index2]*P1r[index2]+P0i[index2]*P1i[index2]));}

}
}

}//k>0
}
}
}
        }

//exit(0);
for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];

if(strcmp(do_anisotropy,"yes") == 0){
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
if(strcmp(do_odd_multipoles,"yes") == 0){
Di[l][0]+=Di[l][tid];
Octo[l][0]+=Octo[l][tid];}
}

nmodes[l][0]+=nmodes[l][tid];
}
}

  f=fopen(name_ps_out,"a");
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}

                if(nmodes[l][0]!=0 && keff>=kmin && keff<=kmax)
                {       K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
if(strcmp(do_anisotropy,"yes") == 0){

                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
                        Hexadeca[l][0]*=9./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
if(strcmp(do_odd_multipoles,"yes") == 0){
                        Di[l][0]*=3./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
                        Octo[l][0]*=7./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);}
}

if(strcmp(do_anisotropy,"no") == 0){
                     fprintf(f,"%lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,nmodes[l][0], P_shot_noise);}

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(do_odd_multipoles,"yes") == 0){
                     fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Di[l][0],Quadru[l][0],Octo[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}

if(strcmp(do_odd_multipoles,"no") == 0){
                     fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}
}


                }
        }
        fclose(f);

freeTokens(K,Nk);
freeTokens(Mono,Nk);

if(strcmp(do_anisotropy,"yes") == 0){
freeTokens(Hexadeca,Nk);
freeTokens(Quadru,Nk);
if(strcmp(do_odd_multipoles,"yes") == 0){
freeTokens(Octo,Nk);
freeTokens(Di,Nk);}
}


if(strcmp(type,"FFT") == 0){free(kvector);}

freeTokensLInt(nmodes,Nk);
//exit(0);
}


void write_power_spectrum_periodic(double kmin, double kmax, double deltak_re[], double deltak_im[], double Deltak, int  ngrid, double L1, double L2, int Ninterlacing, char *name_ps_out, double P_shot_noise, char *binning_type, char *do_odd_multipoles,char *do_anisotropy)
{
double Pi=(4.*atan(1.));
double **K;
double **k_av;
double **Mono;
double **Quadru;
double **Hexadeca;
double **Di;
double **Octo;
long int **nmodes;
int l,tid,i,j,k;
int i2,j2,k2;
long int l2;
FILE *f;
double *kx;
double argument;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=pow(ngrid,3);
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

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


       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

       int nthreads;

  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc((Nk),sizeof(double*));
if(strcmp(do_anisotropy, "yes") == 0){
        Quadru= (double **)calloc((Nk),sizeof(double*));
        Hexadeca=(double **)calloc((Nk),sizeof(double*));}

if(strcmp(do_odd_multipoles, "yes") == 0){
        Di=(double **)calloc((Nk),sizeof(double*));
        Octo=(double **)calloc((Nk),sizeof(double*));}


        K=(double **)calloc((Nk),sizeof(double*));
        k_av=(double **)calloc((Nk),sizeof(double*));
        nmodes=(long int **)calloc((Nk),sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc((nthreads),sizeof(double));
if(strcmp(do_anisotropy, "yes") == 0){
                Quadru[l] = (double*)calloc((nthreads),sizeof(double));
                Hexadeca[l] = (double*)calloc((nthreads),sizeof(double));}
                nmodes[l] = (long int*)calloc((nthreads),sizeof(long int));
                K[l] = (double*)calloc((nthreads),sizeof(double));
                k_av[l] = (double*)calloc((nthreads),sizeof(double));
if(strcmp(do_odd_multipoles, "yes") == 0){
                Di[l] = (double*)calloc((nthreads),sizeof(double));
                Octo[l] = (double*)calloc((nthreads),sizeof(double));}


        }


#pragma  omp parallel for private(index2,l2,i,j,k,l,tid,i2,k2,j2,argument,kmineff,keff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,deltak_re,deltak_im,nmodes,Mono,Di,Quadru,Octo,Hexadeca,Ninterlacing,kmin,bintype_sw,k_av,do_odd_multipoles,do_anisotropy)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(int)(l2/(ngrid*ngrid*1.));
                j=(int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){
l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){
l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
                K[l][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);
if(bintype_sw==0){k_av[l][tid]+=(l+0.5)*Deltak+(kmin);}
if(bintype_sw==1){k_av[l][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

argument=sqrt(kx[k]*kx[k]/(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]));

if(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]==0){argument=0;}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg1(argument);
Octo[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg3(argument);}


nmodes[l][tid]+=1;
}
}
else
{
k2=ngrid-k;
i2=ngrid-i;
j2=ngrid-j;
if(i2==ngrid){i2=0;}
if(j2==ngrid){j2=0;}
if(k2==ngrid){k2=0;}
index2=((pow(ngrid,2)*i2+ngrid*j2+2*k2)/2+i2*ngrid+j2);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg1(argument);
Octo[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg3(argument);}


               nmodes[l][tid]+=1;
}
}
               }//if l<Nk
      }//loop l2





for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
k_av[l][0]+=k_av[l][tid];
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];

if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][0]+=Di[l][tid];
Octo[l][0]+=Octo[l][tid];}


nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");

        for(l=0;l<Nk;l++)
        {
                if(nmodes[l][0]!=0)
                {       k_av[l][0]*=1./nmodes[l][0]*1.;
                        K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=pow(L2-L1,3)/(nmodes[l][0]*1.);
if(strcmp(do_anisotropy, "yes") == 0){
                        Quadru[l][0]*=pow(L2-L1,3)*5./(nmodes[l][0]*1.);
                        Hexadeca[l][0]*=pow(L2-L1,3)*9/(nmodes[l][0]);}
if(strcmp(do_odd_multipoles, "yes") == 0){
                        Di[l][0]*=pow(L2-L1,3)*3./(nmodes[l][0]*1.);
                        Octo[l][0]*=pow(L2-L1,3)*7/(nmodes[l][0]);}


if(k_av[l][0]<=kmax && k_av[l][0]>=kmin){

if(strcmp(do_anisotropy, "yes") == 0){
if(strcmp(do_odd_multipoles, "no") == 0){fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}
if(strcmp(do_odd_multipoles, "yes") == 0){fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Di[l][0],Quadru[l][0],Octo[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}

}
else{
fprintf(f,"%lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,nmodes[l][0], P_shot_noise);}
}

}


                }

        
        fclose(f);



freeTokens(K,Nk);
freeTokens(k_av,Nk);
freeTokens(Mono,Nk);
if(strcmp(do_anisotropy, "yes") == 0){
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);}
if(strcmp(do_odd_multipoles, "yes") == 0){
freeTokens(Di,Nk);
freeTokens(Octo,Nk);}

freeTokensLInt(nmodes,Nk);
free(deltak_re);
free(deltak_im);
free(kx);

}


void write_power_spectrum_periodic2D(double kmin, double kmax, double deltak_re[], double deltak_im[], double Deltak, int mubin, int  ngrid, double L1, double L2, int Ninterlacing, char *name_ps_out, double P_shot_noise, char *binning_type,char *file_for_mu)
{
char name_ps_out2[2000];
double Pi=(4.*atan(1.));
double ***K;
double ***k_av;
double ***Mu;
double ***mu_av;
double ***Pk;
long int ***nmodes;
int l,tid,i,j,k,lmu;
int i2,j2,k2;
long int l2;
FILE *f;
double *kx;
double argument;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=pow(ngrid,3);
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

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

       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

       int nthreads;

  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Pk= (double ***)calloc((Nk),sizeof(double**));
        K=(double ***)calloc((Nk),sizeof(double**));
        k_av=(double ***)calloc((Nk),sizeof(double**));
        nmodes=(long int ***)calloc((Nk),sizeof(long int**));
        Mu=(double ***)calloc((Nk),sizeof(double**));
        mu_av=(double ***)calloc((Nk),sizeof(double**));


        for(l=0;l<Nk;l++)
        {
                Pk[l] = (double**)calloc((mubin),sizeof(double*));
                nmodes[l] = (long int**)calloc((mubin),sizeof(long int*));
                K[l] = (double**)calloc((mubin),sizeof(double*));
                k_av[l] = (double**)calloc((mubin),sizeof(double*));
                Mu[l] = (double**)calloc((mubin),sizeof(double*));
                mu_av[l] = (double**)calloc((mubin),sizeof(double*));

                for(lmu=0;lmu<mubin;lmu++)
                {

                Pk[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                nmodes[l][lmu] = (long int*)calloc((nthreads),sizeof(long int));
                K[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                k_av[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                Mu[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                mu_av[l][lmu] = (double*)calloc((nthreads),sizeof(double));

                }

        }


#pragma  omp parallel for private(index2,l2,i,j,k,l,lmu,tid,i2,k2,j2,argument,kmineff,keff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,Mu,deltak_re,deltak_im,nmodes,Pk,Ninterlacing,kmin,bintype_sw,mubin,k_av,mu_av)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(int)(l2/(ngrid*ngrid*1.));
                j=(int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){
l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){
l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {

argument=sqrt(kx[k]*kx[k]/(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]));
if( kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k] == 0){argument=0;}

lmu=(int)(argument*mubin*1.);
if(lmu==mubin){lmu=lmu-1;}//mu_low<=mu<mu_high execpt last bin which is mu_low<=mu<=mu_high
//if(lmu<0 || lmu>mubin-1 || argument-argument!= 0){printf("Error, mu=%lf, lmu=%d. Exiting now...\n",argument,lmu);exit(0);}


                K[l][lmu][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);
if(bintype_sw==0){k_av[l][lmu][tid]+=(l+0.5)*Deltak+(kmin);}
if(bintype_sw==1){k_av[l][lmu][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

                Mu[l][lmu][tid]+=argument;
                mu_av[l][lmu][tid]+=(lmu+0.5)*1/mubin;


if(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]==0){argument=0;}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);
if(index2>0)
{
Pk[l][lmu][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
nmodes[l][lmu][tid]+=1;
}
}
else
{
k2=ngrid-k;
i2=ngrid-i;
j2=ngrid-j;
if(i2==ngrid){i2=0;}
if(j2==ngrid){j2=0;}
if(k2==ngrid){k2=0;}
index2=((pow(ngrid,2)*i2+ngrid*j2+2*k2)/2+i2*ngrid+j2);
if(index2>0)
{
Pk[l][lmu][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
               nmodes[l][lmu][tid]+=1;
}
}
               }//l<Nk if
      }//main for loop



for(l=0;l<Nk;l++)
{
for(lmu=0;lmu<mubin;lmu++)
{
for(tid=1;tid<nthreads;tid++)
{
k_av[l][lmu][0]+=k_av[l][lmu][tid];
K[l][lmu][0]+=K[l][lmu][tid];
mu_av[l][lmu][0]+=mu_av[l][lmu][tid];
Mu[l][lmu][0]+=Mu[l][lmu][tid];
Pk[l][lmu][0]+=Pk[l][lmu][tid];
nmodes[l][lmu][0]+=nmodes[l][lmu][tid];
}
}
}



if( strcmp(file_for_mu,"no") == 0 ){  f=fopen(name_ps_out,"a");}

        for(lmu=0;lmu<mubin;lmu++)
        {

if( strcmp(file_for_mu,"yes") == 0 ){sprintf(name_ps_out2,"%s_%d.txt",name_ps_out,lmu);f=fopen(name_ps_out2,"a");}

        for(l=0;l<Nk;l++)
        {

                if(nmodes[l][lmu][0]>0)
                {       
                        k_av[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;
                        K[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;

                        mu_av[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;
                        Mu[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;

                        Pk[l][lmu][0]*=pow(L2-L1,3)/(nmodes[l][lmu][0]*1.);

                 if(k_av[l][lmu][0]<=kmax && k_av[l][lmu][0]>=kmin)
                 {
                   fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][lmu][0],K[l][lmu][0],mu_av[l][lmu][0],Mu[l][lmu][0],Pk[l][lmu][0]-P_shot_noise,nmodes[l][lmu][0], P_shot_noise);
                 }


                }
        }
if(strcmp(file_for_mu,"yes") == 0){fclose(f);}

        }
if(strcmp(file_for_mu,"no") == 0){fclose(f);}



freeTokens3(K,Nk,mubin);
freeTokens3(k_av,Nk,mubin);
freeTokens3(Mu,Nk,mubin);
freeTokens3(mu_av,Nk,mubin);
freeTokens3(Pk,Nk,mubin);
freeTokensLInt3(nmodes,Nk,mubin);
free(deltak_re);
free(deltak_im);
free(kx);

}
