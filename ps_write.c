#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "functions.h"
//#include "structures.h"
#define MAX_VEC (1e7)

void write_power_spectrum_skyscuts(double kmin,double kmax,double kx[], double ky[], double kz[], double P0r[], double P0i[],double P1r[], double P1i[], double P2r[], double P2i[], double P3r[], double P3i[], double P4r[], double P4i[],double P0rB[], double P0iB[],double P1rB[], double P1iB[], double P2rB[], double P2iB[], double P3rB[], double P3iB[], double P4rB[], double P4iB[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22,double I22B, char *name_ps_out,char *name_psAB_out,char *name_psBB_out, double P_shot_noise,double P_shot_noiseB, char *binning_type, char *do_anisotropy, char *do_odd_multipoles, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type,char *type,int N_interlacing, char *type_of_code, char *write_kvectors, char *name_ps_kvectors)
{
FILE *f,*fAB,*fBB;
long int i,j,k,l22;
int l,tid;
int mode_duplicate;
double **K;
double **Mono;
double **Quadru;
double **Di;
double **Hexadeca;
double **Octo;
    double **MonoBB;
    double **QuadruBB;
    double **DiBB;
    double **HexadecaBB;
    double **OctoBB;
    
    double **MonoAB;
    double **QuadruAB;
    double **DiAB;
    double **HexadecaAB;
    double **OctoAB;
    
    double **QuadruBA;
    double **DiBA;
    double **HexadecaBA;
    double **OctoBA;

double ***KXvec,***KYvec,***KZvec;//   KXvec[vector][l-bin][tid]
long int **num_vec;
long int *printed_modes;
long int max_num_vectors=MAX_VEC;//among all threads
long int i_vectors;
char name_ps_bin[2000];
long int **nmodes;

double Pi=(4.*atan(1.));
int bintype_sw;
double kmineff;
double keff;
double KX,KY,KZ;
long int ivec,jvec,kvec;
long int i2vec,j2vec,k2vec;
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
ngrid_eff=(long int)(ngrid*1.);
ngrid_eff=ngrid_eff*ngrid;
ngrid_eff=ngrid_eff*ngrid;
}
if(strcmp(type,"DSY") == 0 || strcmp(type,"DSE") == 0 ){
mode_duplicate=2;
ngrid_eff=(long int)(NGRID);
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

    if(strcmp(type_of_code, "rusticoX") == 0)
    {
                MonoBB= (double **)calloc(Nk,sizeof(double*));

        if(strcmp(do_anisotropy,"yes") == 0){
                QuadruBB= (double **)calloc(Nk,sizeof(double*));
                HexadecaBB= (double **)calloc(Nk,sizeof(double*));
        if(strcmp(do_odd_multipoles,"yes") == 0){
                DiBB= (double **)calloc(Nk,sizeof(double*));
                OctoBB= (double **)calloc(Nk,sizeof(double*));}}
        
                MonoAB= (double **)calloc(Nk,sizeof(double*));

        if(strcmp(do_anisotropy,"yes") == 0){
                QuadruAB= (double **)calloc(Nk,sizeof(double*));
                HexadecaAB= (double **)calloc(Nk,sizeof(double*));
                QuadruBA= (double **)calloc(Nk,sizeof(double*));
                HexadecaBA= (double **)calloc(Nk,sizeof(double*));

        if(strcmp(do_odd_multipoles,"yes") == 0){
                DiAB= (double **)calloc(Nk,sizeof(double*));
                OctoAB= (double **)calloc(Nk,sizeof(double*));
                DiBA= (double **)calloc(Nk,sizeof(double*));
                OctoBA= (double **)calloc(Nk,sizeof(double*));
        }
            
        }
        
    }

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
  
            if(strcmp(type_of_code, "rusticoX") == 0)
            {
                                MonoBB[l] = (double*)calloc(nthreads,sizeof(double));

                if(strcmp(do_anisotropy,"yes") == 0){

                                QuadruBB[l] = (double*)calloc(nthreads,sizeof(double));
                                HexadecaBB[l] = (double*)calloc(nthreads,sizeof(double));

                if(strcmp(do_odd_multipoles,"yes") == 0){
                                DiBB[l] = (double*)calloc(nthreads,sizeof(double));
                                OctoBB[l] = (double*)calloc(nthreads,sizeof(double));}}
                
                                MonoAB[l] = (double*)calloc(nthreads,sizeof(double));

                if(strcmp(do_anisotropy,"yes") == 0){

                                QuadruAB[l] = (double*)calloc(nthreads,sizeof(double));
                                HexadecaAB[l] = (double*)calloc(nthreads,sizeof(double));
                                QuadruBA[l] = (double*)calloc(nthreads,sizeof(double));
                                HexadecaBA[l] = (double*)calloc(nthreads,sizeof(double));
                    
                if(strcmp(do_odd_multipoles,"yes") == 0){
                                DiAB[l] = (double*)calloc(nthreads,sizeof(double));
                                OctoAB[l] = (double*)calloc(nthreads,sizeof(double));
                                DiBA[l] = (double*)calloc(nthreads,sizeof(double));
                                OctoBA[l] = (double*)calloc(nthreads,sizeof(double));
                    
                }
                    
                }
                
                
            }
            

                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
        }

if(strcmp(write_kvectors,"yes") == 0){//K[vector][l][thread]

printed_modes = (long int*)calloc(Nk,sizeof(long int));//inizialized to 0

max_num_vectors=(long int)(max_num_vectors*1./nthreads*1.);//split equally among threads.

KXvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));//not inizialized
KYvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));
KZvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));

num_vec = (long int**)calloc(Nk,sizeof(long int*));//inizialized to 0

      for(l=0;l<Nk;l++)
      {
         num_vec[l] = (long int*) calloc(nthreads,sizeof(long int));
      }

        for(i_vectors=0;i_vectors<max_num_vectors;i_vectors++)
        {
                KXvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));
                KYvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));
                KZvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));

                for(l=0;l<Nk;l++)
                {
                    KXvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                    KYvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                    KZvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                }

        }
}

    
#pragma omp parallel for private(i2vec,j2vec,k2vec,ivec,jvec,kvec,KX,KY,KZ,i,tid,l,kmineff,keff,index2) shared(NGRID,ngrid,K,kmin,kx,ky,kz,Mono,Di,Quadru,Octo,Hexadeca,MonoBB,DiBB,QuadruBB,OctoBB,HexadecaBB,MonoAB,DiAB,QuadruAB,OctoAB,HexadecaAB,DiBA,QuadruBA,OctoBA,HexadecaBA,P0r,P1r,P2r,P3r,P4r,P0i,P1i,P2i,P3i,P4i,nmodes,Deltak,Nk,bintype_sw,type,do_anisotropy,do_odd_multipoles, Quadrupole_type,Hexadecapole_type,Octopole_type,N_interlacing,kvector,mode_duplicate,ngrid_eff,type_of_code,write_kvectors,KXvec,KYvec,KZvec,num_vec,max_num_vectors)
        for(i=0;i<ngrid_eff;i++)
        {
                tid=omp_get_thread_num();

if(strcmp(type,"FFT") == 0){
                ivec=(long int)(i/(ngrid*ngrid*1.));
                jvec=(long int)( (i-ivec*ngrid*ngrid)/(ngrid*1.));
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
//if(l<0){l=0;}
if( (pow(KX*KX+KY*KY+KZ*KZ,0.5)-kmin)/Deltak<0  ){l=-1;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(KX*KX+KY*KY+KZ*KZ,0.5))-log10(kmin))/Deltak);//-1.;
//if(l<0){l=0;}
if( ( log10(pow(KX*KX+KY*KY+KZ*KZ,0.5))-log10(kmin))/Deltak<0  ){l=-1;}
}

                if(l<Nk && l>=0 && KX*KX+KY*KY+KZ*KZ>0)
                {
if(KZ==0)
{
 
if(strcmp(type,"DSE") == 0 && KX*KX+KY*KY+KZ*KZ>0){

               K[l][tid]+=pow(KX*KX+KY*KY+KZ*KZ,0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=P0r[i];
               nmodes[l][tid]+=1;

               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors ){
               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=KZ; 
               num_vec[l][tid]++;}

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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
                       MonoBB[l][tid]+=(pow(P0rB[index2],2)+pow(P0iB[index2],2));
                       MonoAB[l][tid]+=(P0r[index2]*P0rB[index2])+(P0i[index2]*P0iB[index2]);//im part non-zero
        
    }
               nmodes[l][tid]+=1;//printf("\n(1)\n");

               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors ){
               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=KZ;
               num_vec[l][tid]++;}


//if(l==0){printf("A %lf %lf %lf, %lf %lf\n",KX,KY,KZ,P0r[index2],P0i[index2]);}

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(Quadrupole_type,"L0L2")== 0){
    Quadru[l][tid]+=(P0r[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0i[index2]*0.5*(3.*P2i[index2]-P0i[index2]));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
     QuadruBB[l][tid]+=P0rB[index2]*0.5*(3.*P2rB[index2]-P0rB[index2])+P0iB[index2]*0.5*(3.*P2iB[index2]-P0iB[index2]);
     QuadruAB[l][tid]+=P0r[index2]*0.5*(3.*P2rB[index2]-P0rB[index2])+P0i[index2]*0.5*(3.*P2iB[index2]-P0iB[index2]);
     QuadruBA[l][tid]+=P0rB[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0iB[index2]*0.5*(3.*P2i[index2]-P0i[index2]);
    }
    
}
if(strcmp(Quadrupole_type,"L1L1")== 0){
    Quadru[l][tid]+=3./2.*(pow(P1r[index2],2)+pow(P1i[index2],2))-0.5*(pow(P0r[index2],2)+pow(P0i[index2],2));
    
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    QuadruBB[l][tid]+=3./2.*(pow(P1rB[index2],2)+pow(P1iB[index2],2))-0.5*(pow(P0rB[index2],2)+pow(P0iB[index2],2));
    QuadruAB[l][tid]+=3./2.*(P1r[index2]*P1rB[index2]+P1i[index2]*P1iB[index2])-0.5*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]);
    QuadruBA[l][tid]+=3./2.*(P1r[index2]*P1rB[index2]+P1i[index2]*P1iB[index2])-0.5*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]);//(identical)
    }
}

if(strcmp(Hexadecapole_type,"L0L4")== 0){
    Hexadeca[l][tid]+=P0r[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0i[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=P0rB[index2]/8.*(35.*P4rB[index2]-30.*P2rB[index2]+3.*P0rB[index2])+P0iB[index2]/8.*(35.*P4iB[index2]-30.*P2iB[index2]+3.*P0iB[index2]);
    HexadecaAB[l][tid]+=P0r[index2]/8.*(35.*P4rB[index2]-30.*P2rB[index2]+3.*P0rB[index2])+P0i[index2]/8.*(35.*P4iB[index2]-30.*P2iB[index2]+3.*P0iB[index2]);
    HexadecaBA[l][tid]+=P0rB[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0iB[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]);
    }
}
if(strcmp(Hexadecapole_type,"L2L2")== 0){
    Hexadeca[l][tid]+=35./8.*(pow(P2r[index2],2)+pow(P2i[index2],2))-30./8.*(P0r[index2]*P2r[index2]+P0i[index2]*P2i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=35./8.*(pow(P2rB[index2],2)+pow(P2iB[index2],2))-30./8.*(P0rB[index2]*P2rB[index2]+P0iB[index2]*P2iB[index2])+3./8.*(pow(P0rB[index2],2)+pow(P0iB[index2],2));
    HexadecaAB[l][tid]+=35./8.*(P2r[index2]*P2rB[index2]+P2i[index2]*P2iB[index2])-30./8.*(P0r[index2]*P2rB[index2]+P0i[index2]*P2iB[index2])+3./8.*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]);//(is non-symmetric)
   HexadecaBA[l][tid]+=35./8.*(P2rB[index2]*P2r[index2]+P2iB[index2]*P2i[index2])-30./8.*(P0rB[index2]*P2r[index2]+P0iB[index2]*P2i[index2])+3./8.*(P0rB[index2]*P0r[index2]+P0iB[index2]*P0i[index2]);
    }
    
}
if(strcmp(Hexadecapole_type,"L1L3")== 0){
    Hexadeca[l][tid]+=35./8.*(P3r[index2]*P1r[index2]+P3i[index2]*P1i[index2])-30./8.*(P1r[index2]*P1r[index2]+P1i[index2]*P1i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=35./8.*(P3rB[index2]*P1rB[index2]+P3iB[index2]*P1iB[index2])-30./8.*(P1rB[index2]*P1rB[index2]+P1iB[index2]*P1iB[index2])+3./8.*(pow(P0rB[index2],2)+pow(P0iB[index2],2));
    HexadecaAB[l][tid]+=35./8.*(P3rB[index2]*P1r[index2]+P3iB[index2]*P1i[index2])-30./8.*(P1rB[index2]*P1r[index2]+P1iB[index2]*P1i[index2])+3./8.*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]);
    HexadecaBA[l][tid]+=35./8.*(P3r[index2]*P1rB[index2]+P3i[index2]*P1iB[index2])-30./8.*(P1r[index2]*P1rB[index2]+P1i[index2]*P1iB[index2])+3./8.*(P0rB[index2]*P0r[index2]+P0iB[index2]*P0i[index2]);

    }
    
}
if(strcmp(do_odd_multipoles,"yes") == 0){
Di[l][tid]+=P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    DiBB[l][tid]+=P0iB[index2]*P1rB[index2]-P0rB[index2]*P1iB[index2];
    DiAB[l][tid]+=P0i[index2]*P1rB[index2]-P0r[index2]*P1iB[index2];
    DiBA[l][tid]+=P0iB[index2]*P1r[index2]-P0rB[index2]*P1i[index2];


    }
if(strcmp(Octopole_type,"L0L3") == 0){
    Octo[l][tid]+=0.5*P0i[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0r[index2]*(5.*P3i[index2]-3.*P1i[index2]);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    OctoBB[l][tid]+=0.5*P0iB[index2]*(5.*P3rB[index2]-3.*P1rB[index2])-0.5*P0rB[index2]*(5.*P3iB[index2]-3.*P1iB[index2]);
    OctoAB[l][tid]+=0.5*P0i[index2]*(5.*P3rB[index2]-3.*P1rB[index2])-0.5*P0r[index2]*(5.*P3iB[index2]-3.*P1iB[index2]);
    OctoBA[l][tid]+=0.5*P0iB[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0rB[index2]*(5.*P3i[index2]-3.*P1i[index2]);
    }
    
}
if(strcmp(Octopole_type,"L1L2") == 0){
    Octo[l][tid]+=5./2.*(P2i[index2]*P1r[index2]-P2r[index2]*P1i[index2])-3./2.*(P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2]);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    OctoBB[l][tid]+=5./2.*(P2iB[index2]*P1rB[index2]-P2rB[index2]*P1iB[index2])-3./2.*(P0iB[index2]*P1rB[index2]-P0rB[index2]*P1iB[index2]);
    OctoAB[l][tid]+=5./2.*(P2iB[index2]*P1r[index2]-P2rB[index2]*P1i[index2])-3./2.*(P0i[index2]*P1rB[index2]-P0r[index2]*P1iB[index2]);
    OctoBA[l][tid]+=5./2.*(P2i[index2]*P1rB[index2]-P2r[index2]*P1iB[index2])-3./2.*(P0iB[index2]*P1r[index2]-P0rB[index2]*P1i[index2]);
    }
    
}
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
              

               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors-mode_duplicate+1 ){

               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=KZ;
               num_vec[l][tid]++;
               if(mode_duplicate==2)
               {
               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=-KZ;//check this
               num_vec[l][tid]++;
               }

               }

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
               
               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors-mode_duplicate+1 ){
               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=KZ;
               num_vec[l][tid]++;

               if(mode_duplicate==2)
               {
               KXvec[num_vec[l][tid]][l][tid]=KX;
               KYvec[num_vec[l][tid]][l][tid]=KY;
               KZvec[num_vec[l][tid]][l][tid]=-KZ;//check this
               num_vec[l][tid]++;
               }

               }

               Mono[l][tid]+=mode_duplicate*(pow(P0r[index2],2)+pow(P0i[index2],2));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    MonoBB[l][tid]+=mode_duplicate*(pow(P0rB[index2],2)+pow(P0iB[index2],2));
    MonoAB[l][tid]+=mode_duplicate*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]);

    }

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(Quadrupole_type,"L0L2")== 0){
    Quadru[l][tid]+=mode_duplicate*(P0r[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0i[index2]*0.5*(3.*P2i[index2]-P0i[index2]));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    QuadruBB[l][tid]+=mode_duplicate*(P0rB[index2]*0.5*(3.*P2rB[index2]-P0rB[index2])+P0iB[index2]*0.5*(3.*P2iB[index2]-P0iB[index2]));
    QuadruAB[l][tid]+=mode_duplicate*(P0r[index2]*0.5*(3.*P2rB[index2]-P0rB[index2])+P0i[index2]*0.5*(3.*P2iB[index2]-P0iB[index2]));
    QuadruBA[l][tid]+=mode_duplicate*(P0rB[index2]*0.5*(3.*P2r[index2]-P0r[index2])+P0iB[index2]*0.5*(3.*P2i[index2]-P0i[index2]));
    }
    
}
if(strcmp(Quadrupole_type,"L1L1")== 0){
    Quadru[l][tid]+=mode_duplicate*(3./2.*(pow(P1r[index2],2)+pow(P1i[index2],2))-0.5*(pow(P0r[index2],2)+pow(P0i[index2],2)));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    QuadruBB[l][tid]+=mode_duplicate*(3./2.*(pow(P1rB[index2],2)+pow(P1iB[index2],2))-0.5*(pow(P0rB[index2],2)+pow(P0iB[index2],2)));
    QuadruAB[l][tid]+=mode_duplicate*(3./2.*(P1r[index2]*P1rB[index2]+P1i[index2]*P1iB[index2])-0.5*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]));
    QuadruBA[l][tid]+=mode_duplicate*(3./2.*(P1rB[index2]*P1r[index2]+P1iB[index2]*P1i[index2])-0.5*(P0rB[index2]*P0r[index2]+P0iB[index2]*P0i[index2]));

    }
    
}

if(strcmp(Hexadecapole_type,"L0L4")== 0){
    Hexadeca[l][tid]+=mode_duplicate*(P0r[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0i[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=mode_duplicate*(P0rB[index2]/8.*(35.*P4rB[index2]-30.*P2rB[index2]+3.*P0rB[index2])+P0iB[index2]/8.*(35.*P4iB[index2]-30.*P2iB[index2]+3.*P0iB[index2]));
    HexadecaAB[l][tid]+=mode_duplicate*(P0r[index2]/8.*(35.*P4rB[index2]-30.*P2rB[index2]+3.*P0rB[index2])+P0i[index2]/8.*(35.*P4iB[index2]-30.*P2iB[index2]+3.*P0iB[index2]));
    HexadecaBA[l][tid]+=mode_duplicate*(P0rB[index2]/8.*(35.*P4r[index2]-30.*P2r[index2]+3.*P0r[index2])+P0iB[index2]/8.*(35.*P4i[index2]-30.*P2i[index2]+3.*P0i[index2]));

    }
    
}
if(strcmp(Hexadecapole_type,"L2L2")== 0){
    Hexadeca[l][tid]+=mode_duplicate*(35./8.*(pow(P2r[index2],2)+pow(P2i[index2],2))-30./8.*(P0r[index2]*P2r[index2]+P0i[index2]*P2i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2)));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=mode_duplicate*(35./8.*(pow(P2rB[index2],2)+pow(P2iB[index2],2))-30./8.*(P0rB[index2]*P2rB[index2]+P0iB[index2]*P2iB[index2])+3./8.*(pow(P0rB[index2],2)+pow(P0iB[index2],2)));
    HexadecaAB[l][tid]+=mode_duplicate*(35./8.*(P2rB[index2]*P2r[index2]+P2iB[index2]*P2i[index2])-30./8.*(P0rB[index2]*P2r[index2]+P0iB[index2]*P2i[index2])+3./8.*(P0rB[index2]*P0r[index2]+P0iB[index2]*P0i[index2]));
    HexadecaBA[l][tid]+=mode_duplicate*(35./8.*(P2r[index2]*P2rB[index2]+P2i[index2]*P2iB[index2])-30./8.*(P0r[index2]*P2rB[index2]+P0i[index2]*P2iB[index2])+3./8.*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]));

    }
    
}
if(strcmp(Hexadecapole_type,"L1L3")== 0){
    Hexadeca[l][tid]+=mode_duplicate*(35./8.*(P3r[index2]*P1r[index2]+P3i[index2]*P1i[index2])-30./8.*(P1r[index2]*P1r[index2]+P1i[index2]*P1i[index2])+3./8.*(pow(P0r[index2],2)+pow(P0i[index2],2)));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    HexadecaBB[l][tid]+=mode_duplicate*(35./8.*(P3rB[index2]*P1rB[index2]+P3iB[index2]*P1iB[index2])-30./8.*(P1rB[index2]*P1rB[index2]+P1iB[index2]*P1iB[index2])+3./8.*(P0rB[index2]*P0rB[index2]+P0iB[index2]*P0iB[index2]));
    HexadecaAB[l][tid]+=mode_duplicate*(35./8.*(P3rB[index2]*P1r[index2]+P3iB[index2]*P1i[index2])-30./8.*(P1rB[index2]*P1r[index2]+P1iB[index2]*P1i[index2])+3./8.*(P0rB[index2]*P0r[index2]+P0iB[index2]*P0i[index2]));
    HexadecaBA[l][tid]+=mode_duplicate*(35./8.*(P3r[index2]*P1rB[index2]+P3i[index2]*P1iB[index2])-30./8.*(P1r[index2]*P1rB[index2]+P1i[index2]*P1iB[index2])+3./8.*(P0r[index2]*P0rB[index2]+P0i[index2]*P0iB[index2]));

    }
    
}

if(strcmp(do_odd_multipoles,"yes") == 0){

Di[l][tid]+=mode_duplicate*(P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2]);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    DiBB[l][tid]+=mode_duplicate*(P0iB[index2]*P1rB[index2]-P0rB[index2]*P1iB[index2]);
    DiAB[l][tid]+=mode_duplicate*(P0iB[index2]*P1r[index2]-P0rB[index2]*P1i[index2]);
    DiBA[l][tid]+=mode_duplicate*(P0i[index2]*P1rB[index2]-P0r[index2]*P1iB[index2]);

    }
if(strcmp(Octopole_type,"L0L3") == 0){
    Octo[l][tid]+=mode_duplicate*(0.5*P0i[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0r[index2]*(5.*P3i[index2]-3.*P1i[index2]));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    OctoBB[l][tid]+=mode_duplicate*(0.5*P0iB[index2]*(5.*P3rB[index2]-3.*P1rB[index2])-0.5*P0rB[index2]*(5.*P3iB[index2]-3.*P1iB[index2]));
    OctoAB[l][tid]+=mode_duplicate*(0.5*P0i[index2]*(5.*P3rB[index2]-3.*P1rB[index2])-0.5*P0r[index2]*(5.*P3iB[index2]-3.*P1iB[index2]));
    OctoBA[l][tid]+=mode_duplicate*(0.5*P0iB[index2]*(5.*P3r[index2]-3.*P1r[index2])-0.5*P0rB[index2]*(5.*P3i[index2]-3.*P1i[index2]));

    }
    
}
if(strcmp(Octopole_type,"L1L2") == 0){
    Octo[l][tid]+=mode_duplicate*(5./2.*(P2i[index2]*P1r[index2]-P2r[index2]*P1i[index2])-3./2.*(P0i[index2]*P1r[index2]-P0r[index2]*P1i[index2]));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    OctoBB[l][tid]+=mode_duplicate*(5./2.*(P2iB[index2]*P1rB[index2]-P2rB[index2]*P1iB[index2])-3./2.*(P0iB[index2]*P1rB[index2]-P0rB[index2]*P1iB[index2]));
    OctoAB[l][tid]+=mode_duplicate*(5./2.*(P2iB[index2]*P1r[index2]-P2rB[index2]*P1i[index2])-3./2.*(P0iB[index2]*P1r[index2]-P0rB[index2]*P1i[index2]));
    OctoBA[l][tid]+=mode_duplicate*(5./2.*(P2i[index2]*P1rB[index2]-P2r[index2]*P1iB[index2])-3./2.*(P0i[index2]*P1rB[index2]-P0r[index2]*P1iB[index2]));

    }
    
}

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

if( strcmp(write_kvectors,"yes") == 0 ){
printed_modes[l]=printed_modes[l]+num_vec[l][tid];//missing thread 0!
}

K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    MonoBB[l][0]+=MonoBB[l][tid];
        MonoAB[l][0]+=MonoAB[l][tid];

    }

if(strcmp(do_anisotropy,"yes") == 0){
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    QuadruBB[l][0]+=QuadruBB[l][tid];
    HexadecaBB[l][0]+=HexadecaBB[l][tid];
        
    QuadruAB[l][0]+=QuadruAB[l][tid];
    HexadecaAB[l][0]+=HexadecaAB[l][tid];
    QuadruBA[l][0]+=QuadruBA[l][tid];
    HexadecaBA[l][0]+=HexadecaBA[l][tid];
    }
if(strcmp(do_odd_multipoles,"yes") == 0){
Di[l][0]+=Di[l][tid];
Octo[l][0]+=Octo[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    DiBB[l][0]+=DiBB[l][tid];
    OctoBB[l][0]+=OctoBB[l][tid];
        DiAB[l][0]+=DiAB[l][tid];
        OctoAB[l][0]+=OctoAB[l][tid];
        DiBA[l][0]+=DiBA[l][tid];
        OctoBA[l][0]+=OctoBA[l][tid];
    }
    
}
}

nmodes[l][0]+=nmodes[l][tid];
}
}


if( strcmp(write_kvectors,"yes") == 0)
{
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}
        printed_modes[l]=printed_modes[l]+num_vec[l][0];//thread 0 was missing!
        sprintf(name_ps_bin,"%s_bin%d.txt",name_ps_kvectors,l);
        f=fopen(name_ps_bin,"w");
        fprintf(f,"#bin %d\n",l); 
        fprintf(f,"#kfundamental %e\n",2*Pi/(L2-L1));
        fprintf(f,"#binning size %e\n",Deltak);
        fprintf(f,"#binning type %s\n",binning_type);
        fprintf(f,"#K-eff= %e\n",K[l][0]*1./nmodes[l][0]*1.);
        fprintf(f,"#K-center= %e\n",keff);
        fprintf(f,"#Total number of modes %ld\n",nmodes[l][0]);
        fprintf(f,"#Printed number of modes %ld\n",printed_modes[l]);
        fprintf(f,"#kx\t ky\t kz\n");
           for(tid=0;tid<nthreads;tid++)
           {
               for(i_vectors=0;i_vectors<num_vec[l][tid];i_vectors++)
               {
                   fprintf(f,"%e\t %e\t %e\n",KXvec[i_vectors][l][tid],KYvec[i_vectors][l][tid],KZvec[i_vectors][l][tid]);
               }
           }
        
        fclose(f);
       }
free(printed_modes);
freeTokensLInt(num_vec,Nk);
freeTokens3(KXvec,max_num_vectors,Nk);
freeTokens3(KYvec,max_num_vectors,Nk);
freeTokens3(KZvec,max_num_vectors,Nk);
}


  f=fopen(name_ps_out,"a");
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    fAB=fopen(name_psAB_out,"a");
    fBB=fopen(name_psBB_out,"a");
    }
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}

                if(nmodes[l][0]!=0 && keff>=kmin && keff<=kmax)
                {       K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
                    if(strcmp(type_of_code, "rusticoX") == 0)
                    {
                    MonoBB[l][0]*=1./(nmodes[l][0]*1.*I22B*N_interlacing*N_interlacing);
                    MonoAB[l][0]*=1./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);

                    }
if(strcmp(do_anisotropy,"yes") == 0){

                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
                        Hexadeca[l][0]*=9./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    QuadruBB[l][0]*=5./(nmodes[l][0]*1.*I22B*N_interlacing*N_interlacing);
    HexadecaBB[l][0]*=9./(nmodes[l][0]*1.*I22B*N_interlacing*N_interlacing);
        
        QuadruAB[l][0]*=5./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        HexadecaAB[l][0]*=9./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        
        QuadruBA[l][0]*=5./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        HexadecaBA[l][0]*=9./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
    }
if(strcmp(do_odd_multipoles,"yes") == 0){
                        Di[l][0]*=3./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
                        Octo[l][0]*=7./(nmodes[l][0]*1.*I22*N_interlacing*N_interlacing);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    DiBB[l][0]*=3./(nmodes[l][0]*1.*I22B*N_interlacing*N_interlacing);
    OctoBB[l][0]*=7./(nmodes[l][0]*1.*I22B*N_interlacing*N_interlacing);
        
        DiAB[l][0]*=3./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        OctoAB[l][0]*=7./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        
        DiBA[l][0]*=3./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
        OctoBA[l][0]*=7./(nmodes[l][0]*1.*sqrt(I22*I22B)*N_interlacing*N_interlacing);
    }
    
}
}

if(strcmp(do_anisotropy,"no") == 0){
                     fprintf(f,"%lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,nmodes[l][0], P_shot_noise);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    fprintf(fBB,"%lf %lf %lf %ld %lf\n",keff, K[l][0],MonoBB[l][0]-P_shot_noiseB,nmodes[l][0], P_shot_noiseB);
    fprintf(fAB,"%lf %lf %lf %ld 0\n",keff, K[l][0],MonoAB[l][0],nmodes[l][0]);

    }
    
}

if(strcmp(do_anisotropy,"yes") == 0){

if(strcmp(do_odd_multipoles,"yes") == 0){
                     fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Di[l][0],Quadru[l][0],Octo[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    fprintf(fBB,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],MonoBB[l][0]-P_shot_noiseB,DiBB[l][0],QuadruBB[l][0],OctoBB[l][0],HexadecaBB[l][0],nmodes[l][0], P_shot_noiseB);
    fprintf(fAB,"%lf %lf %lf %lf %lf %lf %lf %ld 0\n",keff, K[l][0],MonoAB[l][0],0.5*DiAB[l][0]+0.5*DiBA[l][0],0.5*QuadruAB[l][0]+0.5*QuadruBA[l][0],0.5*OctoAB[l][0]+0.5*OctoBA[l][0],0.5*HexadecaAB[l][0]+0.5*HexadecaBA[l][0],nmodes[l][0]);

    }
    
}

if(strcmp(do_odd_multipoles,"no") == 0){
                     fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    fprintf(fBB,"%lf %lf %lf %lf %lf %ld %lf\n",keff, K[l][0],MonoBB[l][0]-P_shot_noiseB,QuadruBB[l][0],HexadecaBB[l][0],nmodes[l][0], P_shot_noiseB);
    fprintf(fAB,"%lf %lf %lf %lf %lf %ld 0\n",keff, K[l][0],MonoAB[l][0],0.5*QuadruAB[l][0]+0.5*QuadruBA[l][0],0.5*HexadecaAB[l][0]+0.5*HexadecaBA[l][0],nmodes[l][0]);

    }
    
}
}


                }
        }
        fclose(f);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    fclose(fBB);
    fclose(fAB);

    }

freeTokens(K,Nk);
freeTokens(Mono,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    freeTokens(MonoBB,Nk);
    freeTokens(MonoAB,Nk);
    }

if(strcmp(do_anisotropy,"yes") == 0){
freeTokens(Hexadeca,Nk);
freeTokens(Quadru,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    freeTokens(HexadecaBB,Nk);
    freeTokens(QuadruBB,Nk);
        freeTokens(HexadecaAB,Nk);
        freeTokens(QuadruAB,Nk);
        freeTokens(HexadecaBA,Nk);
        freeTokens(QuadruBA,Nk);
    }
if(strcmp(do_odd_multipoles,"yes") == 0){
freeTokens(Octo,Nk);
freeTokens(Di,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
    freeTokens(OctoBB,Nk);
    freeTokens(DiBB,Nk);
    freeTokens(OctoAB,Nk);
    freeTokens(DiAB,Nk);
    freeTokens(OctoBA,Nk);
    freeTokens(DiBA,Nk);
    }
    
}
}


if(strcmp(type,"FFT") == 0){free(kvector);}

freeTokensLInt(nmodes,Nk);
//exit(0);
}


void write_power_spectrum_periodic(double kmin,double kmax,double deltak_re[],double deltak_im[],double deltak_reB[],double deltak_imB[],double Deltak, int ngrid,double L1,double L2,int Ninterlacing,char *name_ps_out, char *name_psAB_out, char *name_psBB_out,double P_shot_noise,double P_shot_noiseB,char *binning_type,char *do_odd_multipoles,char *do_anisotropy,char *type_of_code,double I22,double I22B,char *write_kvectors, char *name_ps_kvectors)
{
double Pi=(4.*atan(1.));
double **K;
double **k_av;
double **Mono;
double **Quadru;
double **Hexadeca;
double **Di;
double **Octo;
    double **MonoBB;
    double **QuadruBB;
    double **HexadecaBB;
    double **DiBB;
    double **OctoBB;
    
    double **MonoAB;
    double **QuadruAB;
    double **HexadecaAB;
    double **DiAB;
    double **OctoAB;

double ***KXvec,***KYvec,***KZvec;//   KXvec[vector][l-bin][tid]
long int **num_vec;
long int *printed_modes;
long int max_num_vectors=MAX_VEC;//among all threads 4pik^2Dk (1+Dk^2/(12k^2)) / kf^3 for each k-bin
long int i_vectors;
char name_ps_bin[2000];
    
long int **nmodes;
long int l,tid,i,j,k;
long int i2,j2,k2;
long int l2;
FILE *f,*fAB,*fBB;
double *kx;
double argument;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=(long int)(ngrid*1.);
ngridtot=ngridtot*ngrid;
ngridtot=ngridtot*ngrid;
int FKP=1;
//Pure periodic cases
if(I22==0){I22=pow(L2-L1,-3);FKP=0;}
if(I22B==0){I22B=pow(L2-L1,-3);}


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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        MonoAB= (double **)calloc((Nk),sizeof(double*));
        MonoBB= (double **)calloc((Nk),sizeof(double*));

    }
if(strcmp(do_anisotropy, "yes") == 0){
        Quadru= (double **)calloc((Nk),sizeof(double*));
        Hexadeca=(double **)calloc((Nk),sizeof(double*));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB= (double **)calloc((Nk),sizeof(double*));
        HexadecaBB=(double **)calloc((Nk),sizeof(double*));

        QuadruAB= (double **)calloc((Nk),sizeof(double*));
        HexadecaAB=(double **)calloc((Nk),sizeof(double*));

    }
    
}

if(strcmp(do_odd_multipoles, "yes") == 0){
        Di=(double **)calloc((Nk),sizeof(double*));
        Octo=(double **)calloc((Nk),sizeof(double*));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiBB=(double **)calloc((Nk),sizeof(double*));
        OctoBB=(double **)calloc((Nk),sizeof(double*));

        DiAB=(double **)calloc((Nk),sizeof(double*));
        OctoAB=(double **)calloc((Nk),sizeof(double*));

    }
    
}


        K=(double **)calloc((Nk),sizeof(double*));
        k_av=(double **)calloc((Nk),sizeof(double*));
        nmodes=(long int **)calloc((Nk),sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc((nthreads),sizeof(double));
            if(strcmp(type_of_code, "rusticoX") == 0)
            {
                MonoAB[l] = (double*)calloc((nthreads),sizeof(double));
                MonoBB[l] = (double*)calloc((nthreads),sizeof(double));

            }
if(strcmp(do_anisotropy, "yes") == 0){
                Quadru[l] = (double*)calloc((nthreads),sizeof(double));
                Hexadeca[l] = (double*)calloc((nthreads),sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB[l] = (double*)calloc((nthreads),sizeof(double));
        HexadecaBB[l] = (double*)calloc((nthreads),sizeof(double));
        
        QuadruAB[l] = (double*)calloc((nthreads),sizeof(double));
        HexadecaAB[l] = (double*)calloc((nthreads),sizeof(double));
    }
    
    
    
}
                nmodes[l] = (long int*)calloc((nthreads),sizeof(long int));

                K[l] = (double*)calloc((nthreads),sizeof(double));
                k_av[l] = (double*)calloc((nthreads),sizeof(double));
if(strcmp(do_odd_multipoles, "yes") == 0){
                Di[l] = (double*)calloc((nthreads),sizeof(double));
                Octo[l] = (double*)calloc((nthreads),sizeof(double));
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiAB[l] = (double*)calloc((nthreads),sizeof(double));
        OctoAB[l] = (double*)calloc((nthreads),sizeof(double));

        DiBB[l] = (double*)calloc((nthreads),sizeof(double));
        OctoBB[l] = (double*)calloc((nthreads),sizeof(double));

    }
    
}


        }

if(strcmp(write_kvectors,"yes") == 0){//K[vector][l][thread]

printed_modes = (long int*)calloc(Nk,sizeof(long int));//inizialized to 0

max_num_vectors=(long int)(max_num_vectors*1./nthreads*1.);//split equally among threads.

KXvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));//not inizialized
KYvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));
KZvec = (double ***) malloc(sizeof(double**)*(max_num_vectors));

num_vec = (long int**)calloc(Nk,sizeof(long int*));//inizialized to 0

      for(l=0;l<Nk;l++)
      {
         num_vec[l] = (long int*) calloc(nthreads,sizeof(long int));
      }

        for(i_vectors=0;i_vectors<max_num_vectors;i_vectors++)
        {
                KXvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));
                KYvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));
                KZvec[i_vectors] = (double **) malloc( sizeof(double*)*(Nk));

                for(l=0;l<Nk;l++)
                {
                    KXvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                    KYvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                    KZvec[i_vectors][l] = (double *) malloc( sizeof(double)*(nthreads));
                }

        }
}

#pragma  omp parallel for private(index2,l2,i,j,k,l,tid,i2,k2,j2,argument,kmineff,keff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,deltak_re,deltak_im,nmodes,Mono,Di,Quadru,Octo,Hexadeca,MonoBB,DiBB,QuadruBB,OctoBB,HexadecaBB,MonoAB,DiAB,QuadruAB,OctoAB,HexadecaAB,Ninterlacing,kmin,bintype_sw,k_av,do_odd_multipoles,do_anisotropy,type_of_code,write_kvectors,KXvec,KYvec,KZvec,num_vec,max_num_vectors)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(long int)(l2/(ngrid*ngrid*1.));
                j=(long int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){
l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
//if(l<0){l=0;}
if( (pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak<0   ){l=-1;}
}

if(bintype_sw==1){
l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
//if(l<0){l=0;}
if( ( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak<0  ){l=-1;}
}

                if(l<Nk && l>=0 && pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)>0 )
                {
                K[l][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);
if(bintype_sw==0){k_av[l][tid]+=(l+0.5)*Deltak+(kmin);}
if(bintype_sw==1){k_av[l][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

argument=kx[k]/sqrt(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]);

if(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]==0){argument=0;}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        MonoAB[l][tid]+=(deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.);
        MonoBB[l][tid]+=pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2);

    }
if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg2(argument);
        HexadecaBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg4(argument);
        
        QuadruAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg2(argument);
        HexadecaAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg4(argument);
    }
    
}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg1(argument);
Octo[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg3(argument);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg1(argument);
        OctoBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg3(argument);

        DiAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg1(argument);
        OctoAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg3(argument);

    }
    
}


nmodes[l][tid]+=1;

               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors){
               KXvec[num_vec[l][tid]][l][tid]=kx[i];
               KYvec[num_vec[l][tid]][l][tid]=kx[j];
               KZvec[num_vec[l][tid]][l][tid]=kx[k];
               num_vec[l][tid]++;
               }

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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        MonoBB[l][tid]+=pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2);

        MonoAB[l][tid]+=(deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.);

    }
if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg2(argument);
        HexadecaBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg4(argument);
        
        QuadruAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg2(argument);
        HexadecaAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg4(argument);
    }
    
}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg1(argument);
Octo[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg3(argument);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg1(argument);
        OctoBB[l][tid]+=(pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2))*Leg3(argument);

        DiAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg1(argument);
        OctoAB[l][tid]+=((deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.))*Leg3(argument);

    }
    
}


               nmodes[l][tid]+=1;
               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][tid]<max_num_vectors){
               KXvec[num_vec[l][tid]][l][tid]=kx[i];
               KYvec[num_vec[l][tid]][l][tid]=kx[j];
               KZvec[num_vec[l][tid]][l][tid]=kx[k];
               num_vec[l][tid]++;
               }

}
}
               }//if l<Nk
      }//loop l2





for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{

if( strcmp(write_kvectors,"yes") == 0 ){
printed_modes[l]=printed_modes[l]+num_vec[l][tid];//missing thread 0!
}

k_av[l][0]+=k_av[l][tid];
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        MonoAB[l][0]+=MonoAB[l][tid];
        MonoBB[l][0]+=MonoBB[l][tid];

    }

if(strcmp(do_anisotropy, "yes") == 0){
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB[l][0]+=QuadruBB[l][tid];
        HexadecaBB[l][0]+=HexadecaBB[l][tid];

        QuadruAB[l][0]+=QuadruAB[l][tid];
        HexadecaAB[l][0]+=HexadecaAB[l][tid];

    }
    
}

if(strcmp(do_odd_multipoles, "yes") == 0){
Di[l][0]+=Di[l][tid];
Octo[l][0]+=Octo[l][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiBB[l][0]+=DiBB[l][tid];
        OctoBB[l][0]+=OctoBB[l][tid];
        
        DiAB[l][0]+=DiAB[l][tid];
        OctoAB[l][0]+=OctoAB[l][tid];
    }
    
}


nmodes[l][0]+=nmodes[l][tid];
}
}


if( strcmp(write_kvectors,"yes") == 0)
{
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}
        printed_modes[l]=printed_modes[l]+num_vec[l][0];//thread 0 was missing!
        sprintf(name_ps_bin,"%s_bin%ld.txt",name_ps_kvectors,l);
        f=fopen(name_ps_bin,"w");
        fprintf(f,"#bin %ld\n",l);
        fprintf(f,"#kfundamental %e\n",2*Pi/(L2-L1));
        fprintf(f,"#binning size %e\n",Deltak);
        fprintf(f,"#binning type %s\n",binning_type);
        fprintf(f,"#K-eff= %e\n",K[l][0]*1./nmodes[l][0]*1.);
        fprintf(f,"#K-center= %e\n",keff);
        fprintf(f,"#Total number of modes %ld\n",nmodes[l][0]);
        fprintf(f,"#Printed number of modes %ld\n",printed_modes[l]);
        fprintf(f,"#kx\t ky\t kz\n");
           for(tid=0;tid<nthreads;tid++)
           {
               for(i_vectors=0;i_vectors<num_vec[l][tid];i_vectors++)
               {
                   fprintf(f,"%e\t %e\t %e\n",KXvec[i_vectors][l][tid],KYvec[i_vectors][l][tid],KZvec[i_vectors][l][tid]);
               }
           }
        
        fclose(f);
        }
free(printed_modes);
freeTokensLInt(num_vec,Nk);
freeTokens3(KXvec,max_num_vectors,Nk);
freeTokens3(KYvec,max_num_vectors,Nk);
freeTokens3(KZvec,max_num_vectors,Nk);
}


  f=fopen(name_ps_out,"a");
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fAB=fopen(name_psAB_out,"a");
        fBB=fopen(name_psBB_out,"a");

    }

        for(l=0;l<Nk;l++)
        {
                if(nmodes[l][0]!=0)
                {       k_av[l][0]*=1./nmodes[l][0]*1.;
                        K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1.0/(nmodes[l][0]*1.*I22);
                    if(strcmp(type_of_code, "rusticoX") == 0)
                    {
                        MonoBB[l][0]*=1.0/(nmodes[l][0]*1.*I22B);
                        MonoAB[l][0]*=1.0/(nmodes[l][0]*1.*sqrt(I22*I22B));

                    }
if(strcmp(do_anisotropy, "yes") == 0){
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        Hexadeca[l][0]*=9./(nmodes[l][0]*1.*I22);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        QuadruBB[l][0]*=5./(nmodes[l][0]*1.*I22B);
        HexadecaBB[l][0]*=9./(nmodes[l][0]*1.*I22B);
        
        QuadruAB[l][0]*=5./(nmodes[l][0]*1.*sqrt(I22*I22B));
        HexadecaAB[l][0]*=9./(nmodes[l][0]*1.*sqrt(I22*I22B));
    }
    
}
if(strcmp(do_odd_multipoles, "yes") == 0){
                        Di[l][0]*=3./(nmodes[l][0]*1.*I22);
                        Octo[l][0]*=7./(nmodes[l][0]*1.*I22);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        DiBB[l][0]*=3./(nmodes[l][0]*1.*I22B);
        OctoBB[l][0]*=7./(nmodes[l][0]*1.*I22B);
        
        DiAB[l][0]*=3./(nmodes[l][0]*1.*sqrt(I22*I22B));
        OctoAB[l][0]*=7./(nmodes[l][0]*1.*sqrt(I22*I22B));
    }
    
}


if(k_av[l][0]<=kmax && k_av[l][0]>=kmin){

if(strcmp(do_anisotropy, "yes") == 0){
if(strcmp(do_odd_multipoles, "no") == 0){fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}
if(strcmp(do_odd_multipoles, "yes") == 0){fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Di[l][0],Quadru[l][0],Octo[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        if(strcmp(do_odd_multipoles, "no") == 0){fprintf(fBB,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],MonoBB[l][0]-P_shot_noiseB,QuadruBB[l][0],HexadecaBB[l][0],nmodes[l][0], P_shot_noiseB);}
        if(strcmp(do_odd_multipoles, "yes") == 0){fprintf(fBB,"%lf %lf %lf %lf %lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],MonoBB[l][0]-P_shot_noiseB,DiBB[l][0],QuadruBB[l][0],OctoBB[l][0],HexadecaBB[l][0],nmodes[l][0], P_shot_noiseB);}
        
        if(strcmp(do_odd_multipoles, "no") == 0){fprintf(fAB,"%lf %lf %lf %lf %lf %ld 0\n",k_av[l][0],K[l][0],MonoAB[l][0],QuadruAB[l][0],HexadecaAB[l][0],nmodes[l][0]);}
        if(strcmp(do_odd_multipoles, "yes") == 0){fprintf(fAB,"%lf %lf %lf %lf %lf %lf %lf %ld 0\n",k_av[l][0],K[l][0],MonoAB[l][0],DiAB[l][0],QuadruAB[l][0],OctoAB[l][0],HexadecaAB[l][0],nmodes[l][0]);}

    }

}
else{
fprintf(f,"%lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,nmodes[l][0], P_shot_noise);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fprintf(fAB,"%lf %lf %lf %ld 0\n",k_av[l][0],K[l][0],MonoAB[l][0],nmodes[l][0]);
        fprintf(fBB,"%lf %lf %lf %ld %lf\n",k_av[l][0],K[l][0],MonoBB[l][0]-P_shot_noiseB,nmodes[l][0], P_shot_noiseB);

    }
    
}
}

}


                }

        
        fclose(f);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fclose(fBB);
        fclose(fAB);

    }



freeTokens(K,Nk);
freeTokens(k_av,Nk);
freeTokens(Mono,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        freeTokens(MonoBB,Nk);
        freeTokens(MonoAB,Nk);

    }
if(strcmp(do_anisotropy, "yes") == 0){
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        freeTokens(QuadruBB,Nk);
        freeTokens(HexadecaBB,Nk);
        freeTokens(QuadruAB,Nk);
        freeTokens(HexadecaAB,Nk);
    }
    
}
if(strcmp(do_odd_multipoles, "yes") == 0){
freeTokens(Di,Nk);
freeTokens(Octo,Nk);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        freeTokens(DiBB,Nk);
        freeTokens(OctoBB,Nk);
        freeTokens(DiAB,Nk);
        freeTokens(OctoAB,Nk);
    }
    
}

freeTokensLInt(nmodes,Nk);

//free these onnly if it's periodic (but not for periodicFKP)
if(FKP==0){

free(deltak_re);
free(deltak_im);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_reB);
        free(deltak_imB);
    }
}

free(kx);

}


void write_power_spectrum_periodic2D(double kmin,double kmax,double deltak_re[],double deltak_im[],double deltak_reB[],double deltak_imB[],double Deltak,int mubin, int ngrid,double L1,double L2,int Ninterlacing,char *name_ps_out,char *name_psAB_out,char *name_psBB_out,double P_shot_noise,double P_shot_noiseB,char *binning_type,char *file_for_mu,char *type_of_code, double I22,double I22B,char *write_kvectors, char *name_ps_kvectors)
{
char name_ps_out2[2000];
char name_psAB_out2[2000];
char name_psBB_out2[2000];

double Pi=(4.*atan(1.));
double ***K;
double ***k_av;
double ***Mu;
double ***mu_av;
double ***Pk;
double ***PkAB;
double ***PkBB;

double ****KXvec,****KYvec,****KZvec;//   KXvec[vector][l-bin][tid]
long int ***num_vec;
long int **printed_modes;
long int max_num_vectors=1e7;//among all threads 4pik^2Dk (1+Dk^2/(12k^2)) / kf^3 for each k-bin
long int i_vectors;
char name_ps_bin[2000];


long int ***nmodes;
long int l,tid,i,j,k,lmu;
long int i2,j2,k2;
long int l2;
FILE *f,*fAB,*fBB;
double *kx;
double argument;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=(long int)(ngrid*1.);
ngridtot=ngridtot*ngrid;
ngridtot=ngridtot*ngrid;
int FKP=1;
if(I22==0){I22=pow(L2-L1,-3);FKP=0;}
if(I22B==0){I22B=pow(L2-L1,-3);}

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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        PkAB= (double ***)calloc((Nk),sizeof(double**));
        PkBB= (double ***)calloc((Nk),sizeof(double**));

    }
        K=(double ***)calloc((Nk),sizeof(double**));
        k_av=(double ***)calloc((Nk),sizeof(double**));
        nmodes=(long int ***)calloc((Nk),sizeof(long int**));
        Mu=(double ***)calloc((Nk),sizeof(double**));
        mu_av=(double ***)calloc((Nk),sizeof(double**));


        for(l=0;l<Nk;l++)
        {
                Pk[l] = (double**)calloc((mubin),sizeof(double*));
            if(strcmp(type_of_code, "rusticoX") == 0)
            {
                PkAB[l] = (double**)calloc((mubin),sizeof(double*));
                PkBB[l] = (double**)calloc((mubin),sizeof(double*));

            }
                nmodes[l] = (long int**)calloc((mubin),sizeof(long int*));
                K[l] = (double**)calloc((mubin),sizeof(double*));
                k_av[l] = (double**)calloc((mubin),sizeof(double*));
                Mu[l] = (double**)calloc((mubin),sizeof(double*));
                mu_av[l] = (double**)calloc((mubin),sizeof(double*));

                for(lmu=0;lmu<mubin;lmu++)
                {

                Pk[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                    if(strcmp(type_of_code, "rusticoX") == 0)
                    {
                        PkAB[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                        PkBB[l][lmu] = (double*)calloc((nthreads),sizeof(double));

                    }
                nmodes[l][lmu] = (long int*)calloc((nthreads),sizeof(long int));
                K[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                k_av[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                Mu[l][lmu] = (double*)calloc((nthreads),sizeof(double));
                mu_av[l][lmu] = (double*)calloc((nthreads),sizeof(double));

                }

        }

if(strcmp(write_kvectors,"yes") == 0){//K[vector][l][thread]

printed_modes = (long int**)calloc(Nk,sizeof(long int*));//inizialized to 0
for(l=0;l<Nk;l++)
{
printed_modes[l]=(long int*)calloc(mubin,sizeof(long int));
}

max_num_vectors=(long int)(max_num_vectors*1./(nthreads*1.*mubin));//split equally among threads.

KXvec = (double ****) malloc(sizeof(double***)*(max_num_vectors));//not inizialized
KYvec = (double ****) malloc(sizeof(double***)*(max_num_vectors));
KZvec = (double ****) malloc(sizeof(double***)*(max_num_vectors));

num_vec = (long int***)calloc(Nk,sizeof(long int**));//inizialized to 0
      for(l=0;l<Nk;l++)
      {
         num_vec[l] = (long int**) calloc(mubin,sizeof(long int*));
         for(lmu=0;lmu<mubin;lmu++)
         {
             num_vec[l][lmu] = (long int*) calloc(nthreads,sizeof(long int));
         }
      }

        for(i_vectors=0;i_vectors<max_num_vectors;i_vectors++)
        {
                KXvec[i_vectors] = (double ***) malloc( sizeof(double**)*(Nk));
                KYvec[i_vectors] = (double ***) malloc( sizeof(double**)*(Nk));
                KZvec[i_vectors] = (double ***) malloc( sizeof(double**)*(Nk));

                for(l=0;l<Nk;l++)
                {
                 
                    KXvec[i_vectors][l] = (double **) malloc( sizeof(double*)*(mubin));
                    KYvec[i_vectors][l] = (double **) malloc( sizeof(double*)*(mubin));
                    KZvec[i_vectors][l] = (double **) malloc( sizeof(double*)*(mubin));
                    for(lmu=0;lmu<mubin;lmu++)
                    {
                         KXvec[i_vectors][l][lmu] = (double *) malloc( sizeof(double)*(nthreads));
                         KYvec[i_vectors][l][lmu] = (double *) malloc( sizeof(double)*(nthreads));
                         KZvec[i_vectors][l][lmu] = (double *) malloc( sizeof(double)*(nthreads));
                    }
                }

        }
}

#pragma  omp parallel for private(index2,l2,i,j,k,l,lmu,tid,i2,k2,j2,argument,kmineff,keff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,Mu,deltak_re,deltak_im,deltak_reB,deltak_imB,nmodes,Pk,PkAB,PkBB,Ninterlacing,kmin,bintype_sw,mubin,k_av,mu_av,type_of_code,write_kvectors,KXvec,KYvec,KZvec,num_vec,max_num_vectors)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(long int)(l2/(ngrid*ngrid*1.));
                j=(long int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){
l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
//if(l<0){l=0;}

if( (pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak<0 ){l=-1;}

}

if(bintype_sw==1){
l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
//if(l<0){l=0;}
if( ( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak<0 ){l=-1;}
}

                if(l<Nk && l>=0 && pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)>0 )
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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        PkAB[l][lmu][tid]+=(deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.);
        PkBB[l][lmu][tid]+=pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2);

    }
nmodes[l][lmu][tid]+=1;
               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][lmu][tid]<max_num_vectors){
               KXvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[i];
               KYvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[j];
               KZvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[k];
               num_vec[l][lmu][tid]++;
               }


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
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        PkBB[l][lmu][tid]+=pow(deltak_reB[index2]/Ninterlacing*1.,2)+pow(deltak_imB[index2]/Ninterlacing*1.,2);
        PkAB[l][lmu][tid]+=(deltak_re[index2]/Ninterlacing*1.)*(deltak_reB[index2]/Ninterlacing*1.)+(deltak_im[index2]/Ninterlacing*1.)*(deltak_imB[index2]/Ninterlacing*1.);

    }
               nmodes[l][lmu][tid]+=1;
               if(strcmp(write_kvectors,"yes") == 0 && num_vec[l][lmu][tid]<max_num_vectors){
               KXvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[i];
               KYvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[j];
               KZvec[num_vec[l][lmu][tid]][l][lmu][tid]=kx[k];

               num_vec[l][lmu][tid]++;
               }

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

if( strcmp(write_kvectors,"yes") == 0 ){
printed_modes[l][lmu]=printed_modes[l][lmu]+num_vec[l][lmu][tid];//missing thread 0!
}

k_av[l][lmu][0]+=k_av[l][lmu][tid];
K[l][lmu][0]+=K[l][lmu][tid];
mu_av[l][lmu][0]+=mu_av[l][lmu][tid];
Mu[l][lmu][0]+=Mu[l][lmu][tid];
Pk[l][lmu][0]+=Pk[l][lmu][tid];
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        PkAB[l][lmu][0]+=PkAB[l][lmu][tid];
        PkBB[l][lmu][0]+=PkBB[l][lmu][tid];

    }
nmodes[l][lmu][0]+=nmodes[l][lmu][tid];
}
}
}


if( strcmp(write_kvectors,"yes") == 0)
{

     for(lmu=0;lmu<mubin;lmu++)
     {
        for(l=0;l<Nk;l++)
        {
//if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
//if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}
        printed_modes[l][lmu]=printed_modes[l][lmu]+num_vec[l][lmu][0];//thread 0 was missing!
        sprintf(name_ps_bin,"%s_bink%ld_binmu%ld.txt",name_ps_kvectors,l,lmu);
        f=fopen(name_ps_bin,"w");
        fprintf(f,"#bin-k %ld\n",l);
        fprintf(f,"#bin-mu %ld\n",lmu);
        fprintf(f,"#kfundamental %e\n",2*Pi/(L2-L1));
        fprintf(f,"#binning-k size %e\n",Deltak);
        fprintf(f,"#binning-mu size %e\n",1./mubin*1.);
        fprintf(f,"#binning type %s\n",binning_type);
        fprintf(f,"#K-eff= %e\n",K[l][lmu][0]*1./nmodes[l][lmu][0]*1.);
        fprintf(f,"#K-center= %e\n",k_av[l][lmu][0]*1./nmodes[l][lmu][0]*1.);
        fprintf(f,"#mu-eff= %e\n",Mu[l][lmu][0]*1./nmodes[l][lmu][0]*1.);
        fprintf(f,"#mu-center= %e\n",mu_av[l][lmu][0]*1./nmodes[l][lmu][0]*1.);
        fprintf(f,"#Total number of modes %ld\n",nmodes[l][lmu][0]);
        fprintf(f,"#Printed number of modes %ld\n",printed_modes[l][lmu]);
        fprintf(f,"#kx\t ky\t kz\n");
           for(tid=0;tid<nthreads;tid++)
           {
               for(i_vectors=0;i_vectors<num_vec[l][lmu][tid];i_vectors++)
               {
                   fprintf(f,"%e\t %e\t %e\n",KXvec[i_vectors][l][lmu][tid],KYvec[i_vectors][l][lmu][tid],KZvec[i_vectors][l][lmu][tid]);
               }
           }
        
     
        fclose(f);
        }

     }

freeTokensLInt(printed_modes,Nk);
freeTokensLInt3(num_vec,Nk,mubin);
freeTokens4(KXvec,max_num_vectors,Nk,mubin);
freeTokens4(KYvec,max_num_vectors,Nk,mubin);
freeTokens4(KZvec,max_num_vectors,Nk,mubin);

}


if( strcmp(file_for_mu,"no") == 0 ){
    f=fopen(name_ps_out,"a");
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fAB=fopen(name_psAB_out,"a");
        fBB=fopen(name_psBB_out,"a");

    }
    
}

        for(lmu=0;lmu<mubin;lmu++)
        {

if( strcmp(file_for_mu,"yes") == 0 ){sprintf(name_ps_out2,"%s_%ld.txt",name_ps_out,lmu);f=fopen(name_ps_out2,"a");}
            if(strcmp(type_of_code, "rusticoX") == 0)
            {
                if( strcmp(file_for_mu,"yes") == 0 ){sprintf(name_psAB_out2,"%s_%ld.txt",name_psAB_out,lmu);fAB=fopen(name_psAB_out2,"a");}
                if( strcmp(file_for_mu,"yes") == 0 ){sprintf(name_psBB_out2,"%s_%ld.txt",name_psBB_out,lmu);fBB=fopen(name_psBB_out2,"a");}

            }

        for(l=0;l<Nk;l++)
        {

                if(nmodes[l][lmu][0]>0)
                {       
                        k_av[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;
                        K[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;

                        mu_av[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;
                        Mu[l][lmu][0]*=1./nmodes[l][lmu][0]*1.;

                        Pk[l][lmu][0]*=1.0/(nmodes[l][lmu][0]*1.*I22);
                    if(strcmp(type_of_code, "rusticoX") == 0)
                    {
                        PkBB[l][lmu][0]*=1.0/(nmodes[l][lmu][0]*1.*I22B);
                        PkAB[l][lmu][0]*=1.0/(nmodes[l][lmu][0]*1.*sqrt(I22*I22B));

                    }

                 if(k_av[l][lmu][0]<=kmax && k_av[l][lmu][0]>=kmin)
                 {
                   fprintf(f,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][lmu][0],K[l][lmu][0],mu_av[l][lmu][0],Mu[l][lmu][0],Pk[l][lmu][0]-P_shot_noise,nmodes[l][lmu][0], P_shot_noise);
                     if(strcmp(type_of_code, "rusticoX") == 0)
                     {
                         fprintf(fAB,"%lf %lf %lf %lf %lf %ld 0\n",k_av[l][lmu][0],K[l][lmu][0],mu_av[l][lmu][0],Mu[l][lmu][0],PkAB[l][lmu][0],nmodes[l][lmu][0]);
                         fprintf(fBB,"%lf %lf %lf %lf %lf %ld %lf\n",k_av[l][lmu][0],K[l][lmu][0],mu_av[l][lmu][0],Mu[l][lmu][0],PkBB[l][lmu][0]-P_shot_noise,nmodes[l][lmu][0], P_shot_noiseB);

                     }
                 }


                }
        }
if(strcmp(file_for_mu,"yes") == 0){
    fclose(f);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fclose(fAB);
        fclose(fBB);

    }
    
}

        }
if(strcmp(file_for_mu,"no") == 0){
    fclose(f);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        fclose(fAB);
        fclose(fBB);

    }
    
}



freeTokens3(K,Nk,mubin);
freeTokens3(k_av,Nk,mubin);
freeTokens3(Mu,Nk,mubin);
freeTokens3(mu_av,Nk,mubin);
freeTokens3(Pk,Nk,mubin);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        freeTokens3(PkAB,Nk,mubin);
        freeTokens3(PkBB,Nk,mubin);

    }
freeTokensLInt3(nmodes,Nk,mubin);

if(FKP==0){

free(deltak_re);
free(deltak_im);
    if(strcmp(type_of_code, "rusticoX") == 0)
    {
        free(deltak_reB);
        free(deltak_imB);

    }
}

free(kx);

}
