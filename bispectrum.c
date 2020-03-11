#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include <string.h> /* memset */
#include <unistd.h> /* close */
#include "fftw_compute.h"
#include "order_algorithm.h"
#include "functions.h"
#include "cubature.h"
#include "ps_write.h"

double Minimum(double a, double b)
{
double min;
min=a;
if(b<a){min=b;}

return min;
}

int Sign(double x)
{
int  f;

if(x>0){f=1;}
if(x==0){f=0;}
if(x<0){f=-1;}

return f;
}

double Number_triangles_exact_equi(double k,double Dk)
{
double triangles;
double Pi=4.*atan(1.);
triangles=(1./3840.)*pow(Pi,2)*(40.*pow(Dk+2*k,3)*(5.*pow(Dk,3)-30.*pow(Dk,2)*k+132.*Dk*pow(k,2)-88.*pow(k,3))-6.*pow(Dk-6.*k,2)*(19.*pow(Dk,3)-22.*pow(Dk,2)*k-188.*Dk*pow(k,2)+376.*pow(k,3))*fabs(Dk-6.*k)-960.*pow(Dk-2.*k,3)*(3.*pow(Dk,2)-10.*Dk*k+12.*pow(k,2))*fabs(Dk- 2.*k)-2538.*pow(Dk,5)*fabs(3.*Dk-2.*k)+8460.*pow(Dk,4)*k*fabs(3.*Dk-2.*k)-5520.*pow(Dk,3)*pow(k,2)*fabs(3.*Dk-2.*k)-9120.*pow(Dk,2)*pow(k,3)*fabs(3.*Dk-2.*k)+12000.*Dk*pow(k,4)*fabs(3.*Dk-2.*k)-3648.*pow(k,5)*fabs(3.*Dk-2.*k)+29.*pow(Dk,6)*Sign(Dk-6.*k)-564.*pow(Dk,5)*k*Sign(Dk-6.*k)+2700.*pow(Dk,4)*pow(k,2)*Sign(Dk-6.*k)+12960.*pow(Dk,3)*pow(k,3)*Sign(Dk-6.*k)-162000.*pow(Dk,2)*pow(k,4)*Sign(Dk-6.*k)+513216.*Dk*pow(k,5)*Sign(Dk-6.*k)-513216.*pow(k,6)*Sign(Dk-6.*k)+2760.*pow(Dk,6)*Sign(Dk-2.*k)-32640.*pow(Dk,5)*k*Sign(Dk-2.*k)+161760.*pow(Dk,4)*pow(k,2)*Sign(Dk-2.*k)-430080.*pow(Dk,3)*pow(k,3)*Sign(Dk-2.*k)+647040.*pow(Dk,2)*pow(k,4)*Sign(Dk-2.*k)-522240.*Dk*pow(k,5)*Sign(Dk-2.*k)+176640.*pow(k,6)*Sign(Dk-2.*k)+pow(3.*Dk-2.*k,4)*(99.*pow(Dk,2)-132.*Dk*k-116.*pow(k,2))*Sign((3.*Dk)/2.-k));

return triangles;
}


double Number_triangles_exact(double k1,double k2,double k3,double Dk)
{
double triangles;
triangles=0;
return triangles;
}
void FFT_k(int ngrid, double *out_k1, double k_input, double Deltakbis, double kf, long int *L,long int *IJK, long int *Kz,double *deltak_re, double *deltak_im, double keff_out[], long int neff_out[], int mode)
{
fftw_plan p;
long int ngridr2c=(pow(ngrid,3)/2+pow(ngrid,2));
long int ngridtot=pow(ngrid,3);
fftw_complex *in_bis;
double ki,ijk;
long int line_test;
long int N_eff1,N_eff1_unique;
double  K_eff1;
//double  K_eff1z;
long int K_eff1_now,K_eff1_bef;
int i,j,k;
long int l_new;
double L_lmu;
double argument;
if(mode == 0){L_lmu=1.0;}

  in_bis= (fftw_complex*) fftw_malloc(ngridr2c*sizeof(fftw_complex));
  memset(in_bis, 0, (ngridr2c)*sizeof(fftw_complex));//set to 0 elements of in_bis1
    fftw_plan_with_nthreads(omp_get_max_threads());
  p=fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid,in_bis,out_k1,FFTW_ESTIMATE);

 
  ki=k_input;
  ijk=pow(ki/kf,2);
  line_test=(long int)(4.4132*pow(ijk,1.4956));
  if( IJK[line_test]>=ijk)
  {
       do{line_test=line_test-1;}while(IJK[line_test]>=ijk);
  }
  else
  {
     line_test=line_test-1;
     do{line_test=line_test+1;}while(IJK[line_test+1]<ijk);
  }
  line_test=line_test+1;
  N_eff1=0;
  N_eff1_unique=0;
  K_eff1=0;
  //K_eff1z=0;

  do{
i=(int)(L[line_test]/(ngrid*ngrid*1.));
j=(int)((L[line_test]-i*(ngrid*ngrid))/(ngrid*1.));
k=L[line_test]-i*ngrid*ngrid-j*ngrid;
if(k<=ngrid/2)
{

argument=Kz[line_test]*1./sqrt(IJK[line_test]*1.);
if(mode==2){L_lmu=Leg2(argument);}

      l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
      __real__ in_bis[l_new]=deltak_re[l_new]*L_lmu;
      __imag__ in_bis[l_new]=deltak_im[l_new]*L_lmu;

if(N_eff1==0)//First iteration
{     K_eff1_now=IJK[line_test];
      K_eff1_bef=-1;
}
else
{
K_eff1_bef=K_eff1_now;
K_eff1_now=IJK[line_test];

if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}
}

      K_eff1+=sqrt(IJK[line_test]*1.)*kf;
      //K_eff1z+=Kz[line_test]*1.*kf;
      N_eff1++;
}
else
{
k=ngrid-k;
i=ngrid-i;
j=ngrid-j;
if(i==ngrid){i=0;}
if(j==ngrid){j=0;}
if(k==ngrid){k=0;}
l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;

argument=Kz[line_test]*1./sqrt(IJK[line_test]*1.);
if(mode==2){L_lmu=Leg2(argument);}


      __real__ in_bis[l_new]=deltak_re[l_new]*L_lmu;//creal
      __imag__ in_bis[l_new]=deltak_im[l_new]*L_lmu;//imag
      K_eff1+=sqrt(IJK[line_test]*1.)*kf;
      //K_eff1z+=Kz[line_test]*1.*kf;

if(N_eff1==0)//First iteration
{     K_eff1_now=IJK[line_test];
      K_eff1_bef=-1;
}
else
{ 
K_eff1_bef=K_eff1_now;
K_eff1_now=IJK[line_test];

if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}
} 

      N_eff1++;
}

line_test++;

   }while( sqrt(IJK[line_test])*kf<Deltakbis+ki );
   K_eff1=K_eff1/N_eff1*1.;
   //K_eff1z=K_eff1z/N_eff1*1.;


    fftw_execute(p);//compute the fft
    fftw_destroy_plan(p);
    fftw_free(in_bis);

keff_out[0]=K_eff1;
//keff_out[1]=K_eff1z;
neff_out[0]=N_eff1;
neff_out[1]=N_eff1_unique+1;
}

void FFT_kt(int ngrid, double *out_k1,double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK, long int *Kz,double *deltak_re, double *deltak_im, double keff_out[], long int neff_out[])
{
fftw_plan p, p_NT;
long int ngridr2c=(pow(ngrid,3)/2+pow(ngrid,2));
fftw_complex *in_bis, *in_bis_NT;
double ki,ijk;
long int line_test;
long int N_eff1;
double  K_eff1;
//double  K_eff1z;
long int K_eff1_now,K_eff1_bef;
long int N_eff1_unique;
int i,j,k;
long int l_new;

  in_bis= (fftw_complex*) fftw_malloc(ngridr2c*sizeof(fftw_complex));
  memset(in_bis, 0, (ngridr2c)*sizeof(fftw_complex));//set to 0 elements of in_bis1
    fftw_plan_with_nthreads(omp_get_max_threads());
  p=fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid,in_bis,out_k1,FFTW_ESTIMATE);

   in_bis_NT= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngridr2c);
   memset(in_bis_NT, 0, ngridr2c*sizeof(fftw_complex));//set to 0 elements of in_bis
    fftw_plan_with_nthreads(omp_get_max_threads());
   p_NT=fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid,in_bis_NT,out_k1_NT,FFTW_ESTIMATE);

  ki=k_input;
  ijk=pow(ki/kf,2);
  line_test=(long int)(4.4132*pow(ijk,1.4956));
  if( IJK[line_test]>=ijk)
  {
       do{line_test=line_test-1;}while(IJK[line_test]>=ijk);
  }
  else
  {
     line_test=line_test-1;
     do{line_test=line_test+1;}while(IJK[line_test+1]<ijk);
  }
  line_test=line_test+1;
  N_eff1=0;
  K_eff1=0;
  //K_eff1z=0;
  N_eff1_unique=0;


  do{
i=(int)(L[line_test]/(ngrid*ngrid*1.));
j=(int)((L[line_test]-i*(ngrid*ngrid))/(ngrid*1.));
k=L[line_test]-i*ngrid*ngrid-j*ngrid;
if(k<=ngrid/2)
{
l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
      __real__ in_bis[l_new]=deltak_re[l_new];
      __imag__ in_bis[l_new]=deltak_im[l_new];
      __real__ in_bis_NT[l_new]=1;
      __imag__ in_bis_NT[l_new]=0;
      K_eff1+=sqrt(IJK[line_test]*1.)*kf;
    //  K_eff1z+=Kz[line_test]*1.*kf;

if(N_eff1==0)//First iteration
{     K_eff1_now=IJK[line_test];
      K_eff1_bef=-1;
}
else
{
K_eff1_bef=K_eff1_now;
K_eff1_now=IJK[line_test];

if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}
}
      N_eff1++;
}
else
{
k=ngrid-k;
i=ngrid-i;
j=ngrid-j;
if(i==ngrid){i=0;}
if(j==ngrid){j=0;}
if(k==ngrid){k=0;}
l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
      __real__ in_bis[l_new]=deltak_re[l_new];//creal
      __imag__ in_bis[l_new]=deltak_im[l_new];//imag
      __real__ in_bis_NT[l_new]=1;
      __imag__ in_bis_NT[l_new]=0;
      K_eff1+=sqrt(IJK[line_test]*1.)*kf;
      //K_eff1z+=Kz[line_test]*1.*kf;

if(N_eff1==0)//First iteration
{     K_eff1_now=IJK[line_test];
      K_eff1_bef=-1;
}
else
{
K_eff1_bef=K_eff1_now;
K_eff1_now=IJK[line_test];

if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}
} 

      N_eff1++;
}

line_test++;

   }while( sqrt(IJK[line_test])*kf<Deltakbis+ki );
   K_eff1=K_eff1/N_eff1*1.;
   //K_eff1z=K_eff1z/N_eff1*1.;

  

    fftw_execute(p);//compute the fft
    fftw_destroy_plan(p);
    fftw_free(in_bis);

    fftw_execute(p_NT);//compute the fft
    fftw_destroy_plan(p_NT);
    fftw_free(in_bis_NT);


keff_out[0]=K_eff1;
//keff_out[1]=K_eff1z;
neff_out[0]=N_eff1;
neff_out[1]=N_eff1_unique+1;
}

void FFT_t(int ngrid, double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK)
{
fftw_plan p_NT;
long int ngridr2c=(pow(ngrid,3)/2+pow(ngrid,2));
fftw_complex *in_bis_NT;
double ki,ijk;
long int line_test;
int i,j,k;
long int l_new;


   in_bis_NT= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngridr2c);
   memset(in_bis_NT, 0, ngridr2c*sizeof(fftw_complex));//set to 0 elements of in_bis
    fftw_plan_with_nthreads(omp_get_max_threads());
   p_NT=fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid,in_bis_NT,out_k1_NT,FFTW_ESTIMATE);

  ki=k_input;
  ijk=pow(ki/kf,2);
  line_test=(long int)(4.4132*pow(ijk,1.4956));
  if( IJK[line_test]>=ijk)
  {
       do{line_test=line_test-1;}while(IJK[line_test]>=ijk);
  }
  else
  {
     line_test=line_test-1;
     do{line_test=line_test+1;}while(IJK[line_test+1]<ijk);
  }
  line_test=line_test+1;
 
  do{
i=(int)(L[line_test]/(ngrid*ngrid*1.));
j=(int)((L[line_test]-i*(ngrid*ngrid))/(ngrid*1.));
k=L[line_test]-i*ngrid*ngrid-j*ngrid;
if(k<=ngrid/2)
{
      l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
      __real__ in_bis_NT[l_new]=1;
      __imag__ in_bis_NT[l_new]=0;
}
else
{
k=ngrid-k;
i=ngrid-i;
j=ngrid-j;
if(i==ngrid){i=0;}
if(j==ngrid){j=0;}
if(k==ngrid){k=0;}
l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
      __real__ in_bis_NT[l_new]=1;
      __imag__ in_bis_NT[l_new]=0;
}

line_test++;

   }while( sqrt(IJK[line_test])*kf<Deltakbis+ki );

    fftw_execute(p_NT);//compute the fft
    fftw_destroy_plan(p_NT);
    fftw_free(in_bis_NT);
}

void write_triangles_bin(int ngrid,double *Triangles_k,long int *Weight_k, double k_input, double Deltakbis, double  kf, long int *L,long int *IJK)
{
double ki,ijk;
long int line_test;
int i,j,k;
long int l_triangles;
long int K_eff1_now,K_eff1_bef;
long int N_eff1;
  ki=k_input;
  ijk=pow(ki/kf,2);
  line_test=(long int)(4.4132*pow(ijk,1.4956));
  if( IJK[line_test]>=ijk)
  {
       do{line_test=line_test-1;}while(IJK[line_test]>=ijk);
  }
  else
  {
     line_test=line_test-1;
     do{line_test=line_test+1;}while(IJK[line_test+1]<ijk);
  }
  line_test=line_test+1;

l_triangles=0;
N_eff1=0;
  do{
i=(int)(L[line_test]/(ngrid*ngrid*1.));
j=(int)((L[line_test]-i*(ngrid*ngrid))/(ngrid*1.));
k=L[line_test]-i*ngrid*ngrid-j*ngrid;
if(k<=ngrid/2)
{
 
if(N_eff1==0)//First iteration
{     K_eff1_now=IJK[line_test];
      K_eff1_bef=-1;
      Triangles_k[l_triangles]=sqrt(IJK[line_test]*1.)*kf;
      Weight_k[l_triangles]=1;
}
else
{
K_eff1_bef=K_eff1_now;
K_eff1_now=IJK[line_test];

    if(K_eff1_bef==K_eff1_now)
    {
     Weight_k[l_triangles]+=1;
    }
    else
    {
      l_triangles++;
      Triangles_k[l_triangles]=sqrt(IJK[line_test]*1.)*kf;
      Weight_k[l_triangles]=1;
    }

}

N_eff1++;

}
line_test++;

   }while( sqrt(IJK[line_test])*kf<Deltakbis+ki );

}


void bispectrum_calculator(double *K1,double **K2, double ***K3, int Nk1, int *Nk2, int **Nk3, double *deltak_re, double *deltak_im,double *deltak2_re, double *deltak2_im,int Ninterlacing, double kmin,double kmax,double L1,double L2, int ngrid, double Deltakbis, double I33,double I22, double IN, double Bsn, double Psn, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *triangles_num, char *write_triangles, char *triangles_id, long int *input_bin, char *do_bispectrum2)
{
//double K1z,K2z,K3z;
long int bin=input_bin[0];
int i1,i2,i3;
int ii1,ii2,ii3;
int i,j,k;
long int l,l_new;
FILE *f,*ftriangles;
FILE *f002,*f020,*f200;
long int ngridr2c,ngridtot;
double Pi=4.*atan(1.0);
double kf,kth;
double B,Q,Q002,Q020,Q200;
double B002,B020,B200;
int periodic;
long int lines_temp_file;
int tid;
long int *L;
long int *IJK;
long int *Kz;
char name_of_order_file[1000];
double *Power, *Number;//,*Number_theo;
double *Power1,*Power2,*Power3;
double *Power11,*Power22,*Power33;

double *Power002,*Power020,*Power200;
double K_eff1,K_eff2,K_eff3;
long int N_eff1,N_eff2,N_eff3;
long int N_eff1_unique,N_eff2_unique,N_eff3_unique;
double ijk;
double ki;
double B_shot_noise,Q_shot_noise,Q_shot_noise002,Q_shot_noise020,Q_shot_noise200;
double B_shot_noise002,B_shot_noise020,B_shot_noise200;
long int line_test;
double Number_triangles;
double *out_k1,*out_k2,*out_k3;
double *out_k1_2,*out_k2_2,*out_k3_2;
double *out_k1_NT,*out_k2_NT,*out_k3_NT;
double P1,P2,P3;
double P11,P22,P33;

double keff_out[1];
long int neff_out[2];
double *Triangles_k1,*Triangles_k2,*Triangles_k3;
long int *Weight_k1,*Weight_k2,*Weight_k3;
double k1_bin,k2_bin,k3_bin;
long int w_k1_bin,w_k2_bin,w_k3_bin,wtot;
long int Number_triangles_norm;
char name_triangles_file[1000];
//double argument1,argument2,argument3;
if(I33==0 && I22==0 && IN==0 && Bsn==0){periodic=1;/*periodic case*/}
else{periodic=0;}

ngridtot=pow(ngrid,3);
ngridr2c=(pow(ngrid,3)/2+pow(ngrid,2));

kf=2.*Pi/(L2-L1);

Power = (double*)malloc((n_lines_parallel)*sizeof(double));

if( strcmp(do_bispectrum2, "yes") == 0 )
{
Power002 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power020 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power200 = (double*)malloc((n_lines_parallel)*sizeof(double));
}


if(strcmp(triangles_num, "FFT") == 0)
{
Number= (double*)malloc((n_lines_parallel)*sizeof(double));
}
Power1 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power2 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power3 = (double*)malloc((n_lines_parallel)*sizeof(double));

if( strcmp(do_bispectrum2, "yes") == 0)
{
Power11 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power22 = (double*)malloc((n_lines_parallel)*sizeof(double));
Power33 = (double*)malloc((n_lines_parallel)*sizeof(double));
}

sprintf(name_of_order_file,"wisdom%d.dat",ngrid);


f=fopen(name_of_order_file,"r");
if(f==NULL)
{
printf("Ordering file for ngrid=%d not found. Creating one...",ngrid);
prepare_pointer_to_delta(ngrid);
printf("Ok!\n");

lines_temp_file=countlines(name_of_order_file);
f=fopen(name_of_order_file,"r");
}
else
{
printf("Ordering file for ngrid=%d found.\n",ngrid);
lines_temp_file=countlines(name_of_order_file);
}

L = (long int*) malloc(lines_temp_file*sizeof(long int));
IJK = (long int*) malloc(lines_temp_file*sizeof(long int));
Kz = (long int*) malloc(lines_temp_file*sizeof(long int));

    for(l=0;l<lines_temp_file;l++)//cut these files at IJK<(N/2)**2
    {
       fscanf(f,"%ld %ld %ld\n",&L[l],&IJK[l],&Kz[l]);
    }
fclose(f);

for(i1=0;i1<Nk1;i1++)
{
   out_k1=(double *) malloc(sizeof(double)*ngridtot);
if( strcmp(do_bispectrum2, "yes") == 0){out_k1_2=(double *) malloc(sizeof(double)*ngridtot);}

   if(strcmp(triangles_num, "FFT") == 0)
   {   
          out_k1_NT=(double*) malloc(sizeof(double)*ngridtot);
          FFT_kt(ngrid,out_k1,out_k1_NT,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out);
   K_eff1=keff_out[0];
//   K1z=keff_out[1];
   N_eff1=neff_out[0];
   N_eff1_unique=neff_out[1];
  if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 0){FFT_k(ngrid,out_k1_2,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
  if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 1){FFT_k(ngrid,out_k1_2,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}


//Alternative option (less RAM memory, more time)
//            FFT_k(ngrid,out_k1,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out);
//            K_eff1=keff_out[0];
//            N_eff1=neff_out[0];
//            N_eff1_unique=neff_out[1];
//            out_k1_NT=(double*) malloc(sizeof(double)*ngridtot);
//            FFT_t(ngrid,out_k1_NT,K1[i1],Deltakbis, kf, L,IJK);

   if(strcmp(write_triangles, "yes") == 0)
   {
          Triangles_k1=(double*) malloc(sizeof(double)*N_eff1_unique);
          Weight_k1=(long int*) malloc(sizeof(long int)*N_eff1_unique);
          write_triangles_bin(ngrid,Triangles_k1,Weight_k1,K1[i1],Deltakbis, kf, L,IJK);
   }

   }
   else
   {
          FFT_k(ngrid,out_k1,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,0);
          K_eff1=keff_out[0];
  //        K1z=keff_out[1];
          N_eff1=neff_out[0];
          N_eff1_unique=neff_out[1];
          if( strcmp(do_bispectrum2, "yes") == 0 && periodic ==0){FFT_k(ngrid,out_k1_2,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
          if( strcmp(do_bispectrum2, "yes") == 0 && periodic ==1){FFT_k(ngrid,out_k1_2,K1[i1],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}


          if(strcmp(write_triangles, "yes") == 0 || strcmp(triangles_num, "APR_SUM") == 0 ||  strcmp(triangles_num, "EXA_SUM") == 0)
          {
          Triangles_k1=(double*) malloc(sizeof(double)*N_eff1_unique);
          Weight_k1=(long int*) malloc(sizeof(long int)*N_eff1_unique);
          write_triangles_bin(ngrid,Triangles_k1,Weight_k1,K1[i1],Deltakbis, kf, L,IJK);    
          }          
   }

f=fopen(name_bs_out,"a");

if( strcmp(do_bispectrum2, "yes") == 0 )
{
f002=fopen(name_bs002_out,"a");
f020=fopen(name_bs020_out,"a");
f200=fopen(name_bs200_out,"a");
}

for(i2=0;i2<Nk2[i1];i2++)
{
   out_k2=(double *) malloc(sizeof(double)*ngridtot);
if( strcmp(do_bispectrum2, "yes") == 0){ out_k2_2=(double *) malloc(sizeof(double)*ngridtot);}

   if(strcmp(triangles_num, "FFT") == 0)
   {      
        out_k2_NT=(double*) malloc(sizeof(double)*ngridtot);
        FFT_kt(ngrid,out_k2,out_k2_NT,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out);
        K_eff2=keff_out[0];
    //    K2z=keff_out[1];
        N_eff2=neff_out[0];
        N_eff2_unique=neff_out[1];
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic ==0){FFT_k(ngrid,out_k2_2,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic ==1){FFT_k(ngrid,out_k2_2,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}




   if(strcmp(write_triangles, "yes") == 0)
   {
          Triangles_k2=(double*) malloc(sizeof(double)*N_eff2_unique);
          Weight_k2=(long int*) malloc(sizeof(long int)*N_eff2_unique);
          write_triangles_bin(ngrid,Triangles_k2,Weight_k2,K2[i1][i2],Deltakbis, kf, L,IJK);
   }

   }
   else
   {
        FFT_k(ngrid,out_k2,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,0);
        K_eff2=keff_out[0];
      //  K2z=keff_out[1];
        N_eff2=neff_out[0];
        N_eff2_unique=neff_out[1];
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 0){FFT_k(ngrid,out_k2_2,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 1){FFT_k(ngrid,out_k2_2,K2[i1][i2],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}


              if(strcmp(write_triangles, "yes") == 0 || strcmp(triangles_num, "APR_SUM") == 0 || strcmp(triangles_num, "EXA_SUM") == 0)
               {
               Triangles_k2=(double*) malloc(sizeof(double)*N_eff2_unique);
               Weight_k2=(long int*) malloc(sizeof(long int)*N_eff2_unique);
               write_triangles_bin(ngrid,Triangles_k2,Weight_k2,K2[i1][i2],Deltakbis, kf, L,IJK); 
               }

   }
  
for(i3=0;i3<Nk3[i1][i2];i3++)
{
   out_k3=(double *) malloc(sizeof(double)*ngridtot);
if( strcmp(do_bispectrum2, "yes") == 0){ out_k3_2=(double *) malloc(sizeof(double)*ngridtot);}

   if(strcmp(triangles_num, "FFT") == 0)
   { 
        out_k3_NT=(double*) malloc(sizeof(double)*ngridtot);
        FFT_kt(ngrid,out_k3,out_k3_NT,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out);
        K_eff3=keff_out[0];
        //K3z=keff_out[1];
        N_eff3=neff_out[0];
        N_eff3_unique=neff_out[1];
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 0){FFT_k(ngrid,out_k3_2,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 1){FFT_k(ngrid,out_k3_2,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}


        if(strcmp(write_triangles, "yes") == 0)
        {
          Triangles_k3=(double*) malloc(sizeof(double)*N_eff3_unique);
          Weight_k3=(long int*) malloc(sizeof(long int)*N_eff3_unique);
          write_triangles_bin(ngrid,Triangles_k3,Weight_k3,K3[i1][i2][i3],Deltakbis, kf, L,IJK);
        }

   }
   else
   {
        FFT_k(ngrid,out_k3,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,0);
        K_eff3=keff_out[0];
        //K3z=keff_out[1];
        N_eff3=neff_out[0];
        N_eff3_unique=neff_out[1];
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 0){FFT_k(ngrid,out_k3_2,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak2_re,deltak2_im, keff_out, neff_out,0);}
        if( strcmp(do_bispectrum2, "yes") == 0 && periodic == 1){FFT_k(ngrid,out_k3_2,K3[i1][i2][i3],Deltakbis, kf, L,IJK,Kz,deltak_re,deltak_im, keff_out, neff_out,2);}


        if(strcmp(write_triangles, "yes") == 0 || strcmp(triangles_num, "APR_SUM") == 0 ||  strcmp(triangles_num, "EXA_SUM") == 0)
        {
        Triangles_k3=(double*) malloc(sizeof(double)*N_eff3_unique);
        Weight_k3=(long int*) malloc(sizeof(long int)*N_eff3_unique);
        write_triangles_bin(ngrid,Triangles_k3,Weight_k3,K3[i1][i2][i3],Deltakbis, kf, L,IJK); 
        }
   }


B=0;
B002=0;
B020=0;
B200=0;

Number_triangles=0;
Number_triangles_norm=0;
P1=0;P2=0;P3=0;
P11=0;P22=0;P33=0;

/*
if( strcmp(do_bispectrum2, "yes") == 0 && periodic==1)
{
argument3=sqrt(K3z*K3z/(K_eff3*K_eff3));if(K_eff3==0){argument3=0;}
argument2=sqrt(K2z*K2z/(K_eff2*K_eff2));if(K_eff2==0){argument2=0;}
argument1=sqrt(K1z*K1z/(K_eff1*K_eff1));if(K_eff1==0){argument1=0;}
}
*/

for(l=0;l<n_lines_parallel;l++)
{

Power[l]=0;
if( strcmp(do_bispectrum2, "yes") == 0 ){
Power002[l]=0;
Power020[l]=0;
Power200[l]=0;
}

if(strcmp(triangles_num, "FFT") == 0)
{
Number[l]=0;
}
Power1[l]=0;
Power2[l]=0;
Power3[l]=0;

if( strcmp(do_bispectrum2, "yes") == 0 ){
Power11[l]=0;
Power22[l]=0;
Power33[l]=0;
}

}

if(periodic==1)//periodic case
{
I33=pow(L2-L1,-6);
I22=pow(L2-L1,-3);
IN=Psn;
Bsn=Psn*Psn;
}


#pragma omp parallel for private(l,tid,i,j,k) shared(ngridtot,Power,Power002,Power020,Power200,Number,out_k1,out_k2,out_k3,out_k1_2,out_k2_2,out_k3_2,I33,out_k1_NT,out_k2_NT,out_k3_NT,Ninterlacing,Power1,Power2,Power3,Power11,Power22,Power33,I22,L2,L1,ngrid,Pi,K_eff1,K_eff2,K_eff3,Deltakbis,kf,do_bispectrum2,periodic)
for(l=0;l<ngridtot;l++)//This needs to go parallel
{
tid=omp_get_thread_num();
Power[tid]=Power[tid]+out_k1[l]*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33);

if( strcmp(do_bispectrum2, "yes") == 0 && periodic==0){
Power002[tid]=Power002[tid]+5.*out_k1[l]*out_k2[l]*(3./2.*out_k3_2[l]-1./2.*out_k3[l])*pow(Ninterlacing*ngrid,-3)/(I33);
Power020[tid]=Power020[tid]+5.*out_k1[l]*(3./2.*out_k2_2[l]-1./2.*out_k2[l])*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33);
Power200[tid]=Power200[tid]+5.*(3./2.*out_k1_2[l]-1./2.*out_k1[l])*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33);
}
if( strcmp(do_bispectrum2, "yes") == 0 && periodic==1){
//Power002[tid]=Power002[tid]+5.*out_k1[l]*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33)*Leg2(argument3);
//Power020[tid]=Power020[tid]+5.*out_k1[l]*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33)*Leg2(argument2);
//Power200[tid]=Power200[tid]+5.*out_k1[l]*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33)*Leg2(argument1);
Power002[tid]=Power002[tid]+5.*out_k1[l]*out_k2[l]*out_k3_2[l]*pow(Ninterlacing*ngrid,-3)/(I33);
Power020[tid]=Power020[tid]+5.*out_k1[l]*out_k2_2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33);
Power200[tid]=Power200[tid]+5.*out_k1_2[l]*out_k2[l]*out_k3[l]*pow(Ninterlacing*ngrid,-3)/(I33);

}


if(strcmp(triangles_num, "FFT") == 0)
{
Number[tid]=Number[tid]+out_k1_NT[l]*out_k2_NT[l]*out_k3_NT[l]*pow(ngrid,-3);
}
Power1[tid]=Power1[tid]+pow(ngrid,-3)*out_k1[l]*out_k1[l]/(I22*Ninterlacing*Ninterlacing);
Power2[tid]=Power2[tid]+pow(ngrid,-3)*out_k2[l]*out_k2[l]/(I22*Ninterlacing*Ninterlacing);
Power3[tid]=Power3[tid]+pow(ngrid,-3)*out_k3[l]*out_k3[l]/(I22*Ninterlacing*Ninterlacing);

if( strcmp(do_bispectrum2, "yes") == 0 && periodic==1){

//Power11[tid]=Power11[tid]+pow(ngrid,-3)*5*out_k1[l]*out_k1[l]/(I22*Ninterlacing*Ninterlacing)*Leg2(argument1);
//Power22[tid]=Power22[tid]+pow(ngrid,-3)*5*out_k2[l]*out_k2[l]/(I22*Ninterlacing*Ninterlacing)*Leg2(argument2);
//Power33[tid]=Power33[tid]+pow(ngrid,-3)*5*out_k3[l]*out_k3[l]/(I22*Ninterlacing*Ninterlacing)*Leg2(argument3);
Power11[tid]=Power11[tid]+pow(ngrid,-3)*5*out_k1[l]*out_k1_2[l]/(I22*Ninterlacing*Ninterlacing);
Power22[tid]=Power22[tid]+pow(ngrid,-3)*5*out_k2[l]*out_k2_2[l]/(I22*Ninterlacing*Ninterlacing);
Power33[tid]=Power33[tid]+pow(ngrid,-3)*5*out_k3[l]*out_k3_2[l]/(I22*Ninterlacing*Ninterlacing);

}

if( strcmp(do_bispectrum2, "yes") == 0 && periodic==0){
Power11[tid]=Power11[tid]+pow(ngrid,-3)*5*(out_k1[l]*0.5*(3.*out_k1_2[l]-out_k1[l]))/(I22*Ninterlacing*Ninterlacing);
Power22[tid]=Power22[tid]+pow(ngrid,-3)*5*(out_k2[l]*0.5*(3.*out_k2_2[l]-out_k2[l]))/(I22*Ninterlacing*Ninterlacing);
Power33[tid]=Power33[tid]+pow(ngrid,-3)*5*(out_k3[l]*0.5*(3.*out_k3_2[l]-out_k3[l]))/(I22*Ninterlacing*Ninterlacing);
}


}



if(strcmp(write_triangles, "yes") == 0)
{
sprintf(name_triangles_file,"%s_bin%ld.txt",triangles_id,bin);
ftriangles=fopen(name_triangles_file,"w");
fprintf(ftriangles,"#Effective %.16lf %.16lf %.16lf\n",K_eff1,K_eff2,K_eff3);
fprintf(ftriangles,"#Centre %.16lf %.16lf %.16lf\n",K1[i1],K2[i1][i2],K3[i1][i2][i3]);
}

if(strcmp(write_triangles, "yes") == 0 || strcmp(triangles_num, "APR_SUM") == 0 ||  strcmp(triangles_num, "EXA_SUM") == 0)
{
for(ii1=0;ii1<N_eff1_unique;ii1++)
{
k1_bin=Triangles_k1[ii1];
w_k1_bin=Weight_k1[ii1];
   for(ii2=0;ii2<N_eff2_unique;ii2++)
   {
      k2_bin=Triangles_k2[ii2];
      w_k2_bin=Weight_k2[ii2];
      for(ii3=0;ii3<N_eff3_unique;ii3++)
      {
          k3_bin=Triangles_k3[ii3];
          w_k3_bin=Weight_k3[ii3];
          wtot=w_k1_bin*w_k2_bin*w_k3_bin;

if(strcmp(write_triangles, "yes") == 0){fprintf(ftriangles,"%.16lf %.16lf %.16lf %ld\n",k1_bin,k2_bin,k3_bin,wtot);}

           if(strcmp(triangles_num, "APR_SUM") == 0)
           {
               Number_triangles+=k1_bin*k2_bin*k3_bin*wtot;
               Number_triangles_norm+=wtot;
           }
           if( strcmp(triangles_num, "EXA_SUM") == 0)
           {
              
              Number_triangles+=Number_triangles_exact(k1_bin,k2_bin,k3_bin,Deltakbis)*wtot;
              if(K1[i1]==K3[i1][i2][i3] && K1[i1]==K2[i1][i2]){Number_triangles+=Number_triangles_exact_equi(k1_bin,Deltakbis)*wtot;}
               Number_triangles_norm+=wtot;

           }

      }
   }
}
if(strcmp(write_triangles, "yes") == 0){fclose(ftriangles);}

}                                    
           if(strcmp(triangles_num, "APR_SUM") == 0)
           {
               Number_triangles*=8.*Pi*Pi*pow(Deltakbis,3)*pow(kf,-6)*(1./Number_triangles_norm*1.);//Dimensionless
           }
           if(strcmp(triangles_num, "APR_EFF") == 0)
           {
               Number_triangles=8.*Pi*Pi*pow(Deltakbis,3)*pow(kf,-6)*K_eff1*K_eff2*K_eff3;//Dimensionless
           }
           if(strcmp(triangles_num, "EXA_EFF") == 0)
           {
             Number_triangles=Number_triangles_exact(K_eff1,K_eff2,K_eff3,Deltakbis)*pow(kf,-6);
             if(K1[i1]==K3[i1][i2][i3] && K1[i1]==K2[i1][i2]){Number_triangles=Number_triangles_exact_equi(K_eff1,Deltakbis)*pow(kf,-6);}
           }
           if(strcmp(triangles_num, "EXA_SUM") == 0)
           {
               Number_triangles*=pow(kf,-6)*(1./Number_triangles_norm*1.);//Dimensionless
           }




for(l=0;l<n_lines_parallel;l++)
{
B=B+Power[l];

if( strcmp(do_bispectrum2, "yes") == 0 ){
B002=B002+Power002[l];
B020=B020+Power020[l];
B200=B200+Power200[l];
}


if(strcmp(triangles_num, "FFT") == 0)
{
Number_triangles=Number_triangles+Number[l];
}
P1+=Power1[l];
P2+=Power2[l];
P3+=Power3[l];

if( strcmp(do_bispectrum2, "yes") == 0 ){

P11+=Power11[l];
P22+=Power22[l];
P33+=Power33[l];
}

}

B=B/Number_triangles;
if( strcmp(do_bispectrum2, "yes") == 0 ){
B200=B200/Number_triangles;
B020=B020/Number_triangles;
B002=B002/Number_triangles;
}
P1=P1/N_eff1-Psn;
P2=P2/N_eff2-Psn;
P3=P3/N_eff3-Psn;


B_shot_noise=IN*( P1 + P2 + P3 )+Bsn;
B=B-B_shot_noise;

P11=P11/N_eff1;
P22=P22/N_eff2;
P33=P33/N_eff3;
B_shot_noise002=IN*( P33 );
B_shot_noise020=IN*( P22 );
B_shot_noise200=IN*( P11 );

if( strcmp(do_bispectrum2, "yes") == 0 ){
B002=B002-B_shot_noise002;
B020=B020-B_shot_noise020;
B200=B200-B_shot_noise200;
}

Q_shot_noise=B_shot_noise/(P1*P2+P1*P3+P2*P3);
Q=B/(P1*P2+P1*P3+P2*P3);

Q_shot_noise002=B_shot_noise002/(P1*P2+P1*P3+P2*P3);
Q002=B002/(P1*P2+P1*P3+P2*P3);

Q_shot_noise020=B_shot_noise020/(P1*P2+P1*P3+P2*P3);
Q020=B020/(P1*P2+P1*P3+P2*P3);

Q_shot_noise200=B_shot_noise200/(P1*P2+P1*P3+P2*P3);
Q200=B200/(P1*P2+P1*P3+P2*P3);


fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",K_eff1, K1[i1], K_eff2,  K2[i1][i2], K_eff3, K3[i1][i2][i3],B, B_shot_noise,Q,Q_shot_noise, Number_triangles);

if( strcmp(do_bispectrum2, "yes") == 0 )
{
fprintf(f002,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",K_eff1, K1[i1], K_eff2,  K2[i1][i2], K_eff3, K3[i1][i2][i3],B002, B_shot_noise002,Q002,Q_shot_noise002, Number_triangles);

//extra conditions to avoid writing duplicated triangles for 020 and 200
if(K2[i1][i2]!=K3[i1][i2][i3]){
fprintf(f020,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",K_eff1, K1[i1], K_eff2,  K2[i1][i2], K_eff3, K3[i1][i2][i3],B020, B_shot_noise020,Q020,Q_shot_noise020, Number_triangles);}

if(K1[i1]!=K2[i1][i2] && K1[i1]!=K3[i1][i2][i3]){
fprintf(f200,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",K_eff1, K1[i1], K_eff2,  K2[i1][i2], K_eff3, K3[i1][i2][i3],B200, B_shot_noise200,Q200,Q_shot_noise200, Number_triangles);}
}

bin++;

//printf("aaaaA");
//if(K1[i1]==K2[i1][i2] && K1[i1]==K3[i1][i2][i3])
//{  
//printf("%lf %lf %lf %lf %lf %e %lf %lf \n", K_eff1, K_eff2, K_eff3,B,B_shot_noise, Number_triangles,P1,Psn);
//fprintf(g,"%lf %lf %lf %lf %lf %lf %lf\n", K_eff1, K_eff2, K_eff3,B,B_shot_noise, Number_triangles,P1);
//bin++;
//}

free(out_k3);
if( strcmp(do_bispectrum2, "yes") == 0){free(out_k3_2);}

if(strcmp(triangles_num, "APR_SUM") == 0 || strcmp(write_triangles, "yes") == 0)
{
free(Triangles_k3);
free(Weight_k3);
}
if(strcmp(triangles_num, "FFT") == 0){free(out_k3_NT);}
}//loop k3

free(out_k2);
if( strcmp(do_bispectrum2, "yes") == 0){free(out_k2_2);}
if(strcmp(triangles_num, "APR_SUM") == 0 || strcmp(write_triangles, "yes") == 0)
{
free(Triangles_k2);
free(Weight_k2);
}
if(strcmp(triangles_num, "FFT") == 0){free(out_k2_NT);}
}//loop k2
fclose(f);

if( strcmp(do_bispectrum2, "yes") == 0 )
{
fclose(f002);
fclose(f020);
fclose(f200);
}


free(out_k1);
if( strcmp(do_bispectrum2, "yes") == 0){free(out_k1_2);}
if(strcmp(triangles_num, "APR_SUM") == 0 || strcmp(write_triangles, "yes") == 0)
{
free(Triangles_k1);
free(Weight_k1);
}
if(strcmp(triangles_num, "FFT") == 0){free(out_k1_NT);}

}//end of k1-loop

free(L);
free(IJK);
free(Kz);
free(Power);

if(strcmp(do_bispectrum2, "yes") == 0)
{
free(Power002);
free(Power020);
free(Power200);
}

free(Power1);
free(Power2);
free(Power3);
if(strcmp(do_bispectrum2, "yes") == 0)
{
free(Power11);
free(Power22);
free(Power33);
}
if(strcmp(triangles_num, "FFT") == 0){free(Number);}

input_bin[0]=bin;
//printf("finish");
//exit(0);

}

        
long int count_triangles(double kmin,double kmax, double Deltak, char *triangle_shapes)
{
long int Nktot;
long int Nk;
long int i1,i2,i3;
long int i;
double K1,K2,K3;
double squeezed_factor=0.1;// number <<1.0
Nk=(long int)((kmax-kmin)/Deltak)+1;//linear and log?
i=0;//Count number of triangles

if(strcmp(triangle_shapes, "ALL") == 0)
{
 for(i1=0;i1<Nk;i1++)
 {
   K1=kmin+(i1+1)*Deltak;
   for(i2=i1;i2<Nk;i2++)
   {
      K2=kmin+(i2+1)*Deltak;
      for(i3=i2;i3<Nk;i3++)
      {
          K3=kmin+(i3+1)*Deltak;
          if(K3<K2+K1){i++;}
      }
   }

 }
}

if(strcmp(triangle_shapes, "EQU") == 0)
{
 for(i1=0;i1<Nk;i1++)
 {
   K1=kmin+(i1+1)*Deltak;
   K2=kmin+(i1+1)*Deltak;
   K3=kmin+(i1+1)*Deltak;
   if(K3<K2+K1){i++;}

 }

}

if(strcmp(triangle_shapes, "ISO") == 0)
{
 for(i1=0;i1<Nk;i1++)
 {
   K1=kmin+(i1+1)*Deltak;
   for(i2=i1;i2<Nk;i2++)
   {
      K2=kmin+(i2+1)*Deltak;
      K3=kmin+(i2+1)*Deltak;
      if(K3<K2+K1){i++;}
      
   }

 }
}

if(strcmp(triangle_shapes, "SQU") == 0)
{
 for(i1=0;i1<Nk;i1++)
 {
   K1=kmin+(i1+1)*Deltak;
   for(i2=i1;i2<Nk;i2++)
   {
      K2=kmin+(i2+1)*Deltak;
      for(i3=i2;i3<Nk;i3++)
      {
          K3=kmin+(i3+1)*Deltak;
          if(K3<K2+K1 && fabs(K3-K2)<=K1 && K1<=squeezed_factor*Minimum(K2,K3) ){i++;}
      }
   }

 }
}

        Nktot=i;
return Nktot;
}
void select_triangles(double kmin,double kmax, double Deltak, double L1, double L2, double *k1, double *k2, double *k3, int *which_grid, long int Nktot,int ngrid_in, char *do_multigrid, char *triangle_shapes)
{
double Pi=4.*atan(1.);
long int Nk;
long int i1,i2,i3;
long int i;
int i_ngrid;
int ngrid;
double K1,K2,K3;
double ktop;
double factor=1.00;//This factor may depend on grid interpolation and interlacing (less than 1)
double squeezed_factor=0.1;// factor for squeezed triangles <<1
int ngrid_max=round(log(ngrid_in*1.)/log(2.));

Nk=(long int)((kmax-kmin)/Deltak)+1;//linear and log?
i=0;
if(strcmp(triangle_shapes, "ALL") == 0)
{
   for(i1=0;i1<Nk;i1++)
   {
     K1=kmin+(i1+1)*Deltak;
     for(i2=i1;i2<Nk;i2++)
     {
       K2=kmin+(i2+1)*Deltak;
       for(i3=i2;i3<Nk;i3++)
       {
          K3=kmin+(i3+1)*Deltak;
          if(K3<K2+K1)
          {
             k1[i]=K1;
             k2[i]=K2;
             k3[i]=K3;
             i++;
          }
       }
     }
   }
}

if(strcmp(triangle_shapes, "EQU") == 0)
{
   for(i1=0;i1<Nk;i1++)
   {
     K1=kmin+(i1+1)*Deltak;
     K2=kmin+(i1+1)*Deltak;
     K3=kmin+(i1+1)*Deltak;
        if(K3<K2+K1)
        {
           k1[i]=K1;
           k2[i]=K2;
           k3[i]=K3;
           i++;
        }
}

}
if(strcmp(triangle_shapes, "ISO") == 0)
{
   for(i1=0;i1<Nk;i1++)
   {
     K1=kmin+(i1+1)*Deltak;
     for(i2=i1;i2<Nk;i2++)
     {
       K2=kmin+(i2+1)*Deltak;
       K3=kmin+(i2+1)*Deltak;
          if(K3<K2+K1)
          {
             k1[i]=K1;
             k2[i]=K2;
             k3[i]=K3;
             i++;
  
       }
     }
   }

}

if(strcmp(triangle_shapes, "SQU") == 0)
{

   for(i1=0;i1<Nk;i1++)
   {
     K1=kmin+(i1+1)*Deltak;
     for(i2=i1;i2<Nk;i2++)
     {
       K2=kmin+(i2+1)*Deltak;
       for(i3=i2;i3<Nk;i3++)
       {
          K3=kmin+(i3+1)*Deltak;
          if(K3<K2+K1 && fabs(K3-K2)<=K1 && K1<=squeezed_factor*Minimum(K2,K3) )
          {
             k1[i]=K1;
             k2[i]=K2;
             k3[i]=K3;
             i++;
          }
       }
     }
   }
}

if(strcmp(do_multigrid, "yes") == 0)
{
     i_ngrid=ngrid_max;
     do
     {
        ngrid=pow(2,i_ngrid);
        ktop=(2.*Pi/(L2-L1)*ngrid/3.)*factor;//2/3 of k-Nyquist times reduction factor

        for(i=0;i<Nktot;i++)
        {

            if(k3[i]+k1[i]+k2[i]+3*Deltak<3*ktop)
            {
               which_grid[i]=ngrid;
            }
         }

     i_ngrid--;
     }while(ngrid>0);

}
else
{
     ngrid=pow(2,ngrid_max);
     ktop=(2.*Pi/(L2-L1)*ngrid/3.)*factor;//2/3 of k-Nyquist times reduction factor

     for(i=0;i<Nktot;i++)
     {
            if(k3[i]+k1[i]+k2[i]+3*Deltak<3*ktop)
            {
                 which_grid[i]=ngrid;
            }
     }
}

}

void loop_bispectrum_skycut_caller(double kmin,double kmax, int Ninterlacing,  double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, int ngrid, double P_shot_noise, double Deltakbis, double I33,double I22, double IN, double Bsn,  double alpha, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *do_multigrid, char *triangle_shapes, char *do_bispectrum2)
{
long int input_bin[1];
int i,j,ngrid_i;
long int l;
double *K1,*K2,*K3;
double *K1sel,*K2sel,*K3sel;
int *grid;
long int Nktot;
long int Nktot_i;
double *deltak_re, *deltak_im;
double *deltak2_re, *deltak2_im;
long int ngridtotr2c;
int Nk1,*Nk2,*I2,**Nk3,**I3;
double *K1eff,**K2eff,***K3eff;
long int i1,i2,i3,ll;
int periodic;
        Nktot=count_triangles(kmin,kmax,Deltakbis,triangle_shapes);
        grid= (int *) calloc(Nktot,sizeof(int));
        K1= (double*) calloc(Nktot,sizeof(double));
        K2= (double*) calloc(Nktot,sizeof(double));
        K3= (double*) calloc(Nktot,sizeof(double));

select_triangles(kmin,kmax,Deltakbis,L1,L2,K1,K2,K3,grid,Nktot,ngrid,do_multigrid,triangle_shapes);

//loop over grid cells
l=0;
input_bin[0]=1;
do
{
l++;
ngrid_i=pow(2,l);
printf("Grid %d...",ngrid_i);
//Count number of triangles of specific grid-size
Nktot_i=0;
for(i=0;i<Nktot;i++)
{
if(ngrid_i==grid[i]){Nktot_i++;}
}
printf("%ld triangles\n",Nktot_i);

if(Nktot_i>0)//There are some triangles associated to that gridsize.
{

        ngridtotr2c=(pow(ngrid_i,3)/2+pow(ngrid_i,2));
        deltak_re = (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im = (double*) calloc(ngridtotr2c,sizeof(double));

//put particles in a grid cell and forier transform it to get deltak(k)
if(I33==0 && I22==0 && IN==0 && Bsn==0 && alpha==0 )//periodic
{
periodic=1;
 loop_interlacing_periodic_for_bispectrum(Ninterlacing, s_x, s_y, s_z, weight, Ndata, L1, L2, ngrid_i, n_lines_parallel, type_of_mass_assigment, mode_correction, deltak_re, deltak_im);
}
else//Cutsky
{
periodic=0;

if( strcmp(do_bispectrum2, "yes") == 0 )
{
        deltak2_re = (double*) calloc(ngridtotr2c,sizeof(double));
        deltak2_im = (double*) calloc(ngridtotr2c,sizeof(double));
}

loop_interlacing_skycut_for_bispectrum(Ninterlacing, s_x, s_y, s_z, weight, Ndata, s_x_ran, s_y_ran, s_z_ran, weight_ran, Nrand, L1, L2, ngrid_i, alpha, n_lines_parallel, type_of_mass_assigment, mode_correction, deltak_re, deltak_im, deltak2_re, deltak2_im,do_bispectrum2);
}

//Compute the bispectrum of selected triangles
K1sel=malloc(sizeof(double)*(Nktot_i));
K2sel=malloc(sizeof(double)*(Nktot_i));
K3sel=malloc(sizeof(double)*(Nktot_i));

j=0;
for(i=0;i<Nktot;i++)
{
   if(ngrid_i==grid[i])
   {
      K1sel[j]=K1[i];
      K2sel[j]=K2[i];
      K3sel[j]=K3[i];
j++;
   }
}


//Re-arrange triangles
Nk1=1;//We know there's at least one triangle
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]!=K1sel[ll+1]){Nk1++;}//Number of different values of K1 explored
}
Nk2 = (int*) calloc(Nk1,sizeof(int));
I2 = (int*) calloc(Nk1,sizeof(int));
for(i1=0;i1<Nk1;i1++){Nk2[i1]=1;I2[i1]=1;}//We know there's at least one triangle

i1=1;
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]!=K2sel[ll+1]){Nk2[i1-1]++;}//Number of different values of K2 explored
   if(K1sel[ll]!=K1sel[ll+1]){i1++;}

}

Nk3 = (int**) calloc(Nk1,sizeof(int*));
for(i1=0;i1<Nk1;i1++)
{
Nk3[i1] = (int*) calloc(Nk2[i1],sizeof(int));
}

for(i1=0;i1<Nk1;i1++)
{
  for(i2=0;i2<Nk2[i1];i2++)
  {
     Nk3[i1][i2]=1;
  }
}

i1=1;
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]==K2sel[ll+1] && K3sel[ll]!=K3sel[ll+1]){ Nk3[i1-1][I2[i1-1]-1]++; }//Number of different values of K3 explored
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]!=K2sel[ll+1]){I2[i1-1]++;}
   if(K1sel[ll]!=K1sel[ll+1]){i1++;}

}

K1eff = (double*) calloc(Nk1,sizeof(double));
K2eff = (double**) calloc(Nk1,sizeof(double*));
K3eff = (double***) calloc(Nk1,sizeof(double**));


for(i1=0;i1<Nk1;i1++)
{
   K2eff[i1] = (double*) calloc(Nk2[i1],sizeof(double));
   K3eff[i1] = (double**) calloc(Nk2[i1],sizeof(double*));
  
   for(i2=0;i2<Nk2[i1];i2++)
   {
      K3eff[i1][i2] = (double*) calloc(Nk3[i1][i2],sizeof(double));
   }

}



K1eff[0]=K1sel[0];
i1=1;
for(ll=1;ll<Nktot_i;ll++)
{
   if(K1sel[ll]!=K1sel[ll-1]){K1eff[i1]=K1sel[ll];i1++;}
}

K2eff[0][0]=K2sel[0];
for(i1=0;i1<Nk1;i1++){I2[i1]=1;}
i1=0;
for(ll=1;ll<Nktot_i;ll++)
{

if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]!=K2sel[ll-1]){K2eff[i1][I2[i1]]=K2sel[ll];I2[i1]++;}
if(K1sel[ll]!=K1sel[ll-1]){i1++;K2eff[i1][0]=K2sel[ll];}

}


I3 = (int**) calloc(Nk1,sizeof(int*));
for(i1=0;i1<Nk1;i1++)
{
I3[i1] = (int*) calloc(Nk2[i1],sizeof(int));
}


for(i1=0;i1<Nk1;i1++)
{
  for(i2=0;i2<Nk2[i1];i2++)
  {
     I3[i1][i2]=1;
  }
}

for(i1=0;i1<Nk1;i1++){I2[i1]=0;}
i1=0;
K3eff[0][0][0]=K3sel[0];
for(ll=1;ll<Nktot_i;ll++)
{
   if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]==K2sel[ll-1] && K3sel[ll]!=K3sel[ll-1]){K3eff[i1][I2[i1]][I3[i1][I2[i1]]]=K3sel[ll]; I3[i1][I2[i1]]++;}
   if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]!=K2sel[ll-1]){I2[i1]++;K3eff[i1][I2[i1]][0]=K3sel[ll];}
   if(K1sel[ll]!=K1sel[ll-1]){i1++;K3eff[i1][0][0]=K3sel[ll];}
}


free(K1sel);
free(K2sel);
free(K3sel);

//compute and write bispectrum with selected triangles.
bispectrum_calculator(K1eff, K2eff, K3eff, Nk1, Nk2, Nk3, deltak_re, deltak_im,  deltak2_re, deltak2_im,Ninterlacing, kmin, kmax, L1, L2, ngrid_i, Deltakbis, I33,I22, IN, Bsn,P_shot_noise, n_lines_parallel, binning_type, name_bs_out,name_bs002_out,name_bs020_out,name_bs200_out,triangles_num, write_triangles, triangles_id, input_bin,do_bispectrum2);


free(K1eff);
freeTokens(K2eff,Nk1);
freeTokens2(K3eff,Nk1,Nk2);

free(Nk2);
free(I2);
freeTokensInt(Nk3,Nk1);
freeTokensInt(I3,Nk1);

//free
free(deltak_re);
free(deltak_im);

if( periodic == 0 && strcmp(do_bispectrum2, "yes") == 0)
{
free(deltak2_re);
free(deltak2_im);
}

}


}while(ngrid_i<ngrid);//check if more conditions apply


free(s_x);
free(s_y);
free(s_z);
free(weight);

if(periodic==0)//no periodic case
{
free(s_x_ran);
free(s_y_ran);
free(s_z_ran);
free(weight_ran);
}
free(K1);
free(K2);
free(K3);
free(grid);
}


void loop_bispectrum_periodic_for_gadget_caller(double kmin,double kmax, int Ninterlacing, double L1, double L2, int ngrid, double Deltakbis, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *name_data_in,int gadget_files, char *do_multigrid, char *triangle_shapes,char *RSD, char *do_bispectrum2)
{
long int input_bin[1];
int i,j,ngrid_i;
long int l;
double *K1,*K2,*K3;
double *K1sel,*K2sel,*K3sel;
int *grid;
long int Nktot;
long int Nktot_i;
double *deltak_re, *deltak_im;
long int ngridtotr2c;
int Nk1,*Nk2,*I2,**Nk3,**I3;
double *K1eff,**K2eff,***K3eff;
long int i1,i2,i3,ll;
double params[1];
double Power_spectrum_shot_noise;

        Nktot=count_triangles(kmin,kmax,Deltakbis,triangle_shapes);
        grid= (int *) calloc(Nktot,sizeof(int));
        K1= (double*) calloc(Nktot,sizeof(double));
        K2= (double*) calloc(Nktot,sizeof(double));
        K3= (double*) calloc(Nktot,sizeof(double));

select_triangles(kmin,kmax,Deltakbis,L1,L2,K1,K2,K3,grid,Nktot,ngrid,do_multigrid,triangle_shapes);
l=0;
input_bin[0]=1;
do
{
l++;
ngrid_i=pow(2,l);
printf("Grid %d...",ngrid_i);
Nktot_i=0;
for(i=0;i<Nktot;i++)
{
if(ngrid_i==grid[i]){Nktot_i++;}
}
printf("%ld triangles\n",Nktot_i);

if(Nktot_i>0)//There are some triangles associated to that gridsize.
{

        ngridtotr2c=(pow(ngrid_i,3)/2+pow(ngrid_i,2));
        deltak_re = (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im = (double*) calloc(ngridtotr2c,sizeof(double));


 loop_interlacing_periodic_gadget_for_bispectrum(Ninterlacing, L1, L2, ngrid_i, n_lines_parallel, type_of_mass_assigment, mode_correction, deltak_re, deltak_im, name_data_in,gadget_files,params,RSD);

Power_spectrum_shot_noise=params[0];

K1sel=malloc(sizeof(double)*(Nktot_i));
K2sel=malloc(sizeof(double)*(Nktot_i));
K3sel=malloc(sizeof(double)*(Nktot_i));

j=0;
for(i=0;i<Nktot;i++)
{
   if(ngrid_i==grid[i])
   {
      K1sel[j]=K1[i];
      K2sel[j]=K2[i];
      K3sel[j]=K3[i];
j++;
   }
}
Nk1=1;//We know there's at least one triangle
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]!=K1sel[ll+1]){Nk1++;}//Number of different values of K1 explored
}
Nk2 = (int*) calloc(Nk1,sizeof(int));
I2 = (int*) calloc(Nk1,sizeof(int));
for(i1=0;i1<Nk1;i1++){Nk2[i1]=1;I2[i1]=1;}//We know there's at least one triangle

i1=1;
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]!=K2sel[ll+1]){Nk2[i1-1]++;}//Number of different values of K2 explored
   if(K1sel[ll]!=K1sel[ll+1]){i1++;}

}

Nk3 = (int**) calloc(Nk1,sizeof(int*));
for(i1=0;i1<Nk1;i1++)
{
Nk3[i1] = (int*) calloc(Nk2[i1],sizeof(int));
}

for(i1=0;i1<Nk1;i1++)
{
  for(i2=0;i2<Nk2[i1];i2++)
  {
     Nk3[i1][i2]=1;
  }
}

i1=1;
for(ll=0;ll<Nktot_i-1;ll++)
{
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]==K2sel[ll+1] && K3sel[ll]!=K3sel[ll+1]){ Nk3[i1-1][I2[i1-1]-1]++; }//Number of different values of K3 explored
   if(K1sel[ll]==K1sel[ll+1] && K2sel[ll]!=K2sel[ll+1]){I2[i1-1]++;}
   if(K1sel[ll]!=K1sel[ll+1]){i1++;}

}

K1eff = (double*) calloc(Nk1,sizeof(double));
K2eff = (double**) calloc(Nk1,sizeof(double*));
K3eff = (double***) calloc(Nk1,sizeof(double**));


for(i1=0;i1<Nk1;i1++)
{
   K2eff[i1] = (double*) calloc(Nk2[i1],sizeof(double));
   K3eff[i1] = (double**) calloc(Nk2[i1],sizeof(double*));

   for(i2=0;i2<Nk2[i1];i2++)
   {
      K3eff[i1][i2] = (double*) calloc(Nk3[i1][i2],sizeof(double));

   }

}



K1eff[0]=K1sel[0];
i1=1;
for(ll=1;ll<Nktot_i;ll++)
{
   if(K1sel[ll]!=K1sel[ll-1]){K1eff[i1]=K1sel[ll];i1++;}
}

K2eff[0][0]=K2sel[0];
for(i1=0;i1<Nk1;i1++){I2[i1]=1;}
i1=0;
for(ll=1;ll<Nktot_i;ll++)
{

if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]!=K2sel[ll-1]){K2eff[i1][I2[i1]]=K2sel[ll];I2[i1]++;}
if(K1sel[ll]!=K1sel[ll-1]){i1++;K2eff[i1][0]=K2sel[ll];}

}


I3 = (int**) calloc(Nk1,sizeof(int*));
for(i1=0;i1<Nk1;i1++)
{
I3[i1] = (int*) calloc(Nk2[i1],sizeof(int));
}


for(i1=0;i1<Nk1;i1++)
{
  for(i2=0;i2<Nk2[i1];i2++)
  {
     I3[i1][i2]=1;
  }
}

for(i1=0;i1<Nk1;i1++){I2[i1]=0;}
i1=0;
K3eff[0][0][0]=K3sel[0];
for(ll=1;ll<Nktot_i;ll++)
{
   if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]==K2sel[ll-1] && K3sel[ll]!=K3sel[ll-1]){K3eff[i1][I2[i1]][I3[i1][I2[i1]]]=K3sel[ll]; I3[i1][I2[i1]]++;}
   if(K1sel[ll]==K1sel[ll-1] && K2sel[ll]!=K2sel[ll-1]){I2[i1]++;K3eff[i1][I2[i1]][0]=K3sel[ll];}
   if(K1sel[ll]!=K1sel[ll-1]){i1++;K3eff[i1][0][0]=K3sel[ll];}
}


free(K1sel);
free(K2sel);
free(K3sel);

bispectrum_calculator(K1eff, K2eff, K3eff, Nk1, Nk2, Nk3, deltak_re, deltak_im, NULL,NULL,Ninterlacing, kmin, kmax, L1, L2, ngrid_i, Deltakbis, 0, 0, 0, 0, Power_spectrum_shot_noise, n_lines_parallel, binning_type, name_bs_out,name_bs002_out,name_bs020_out,name_bs200_out,triangles_num, write_triangles, triangles_id, input_bin,do_bispectrum2);


free(K1eff);
freeTokens(K2eff,Nk1);
freeTokens2(K3eff,Nk1,Nk2);

free(Nk2);
free(I2);
freeTokensInt(Nk3,Nk1);
freeTokensInt(I3,Nk1);
free(deltak_re);
free(deltak_im);


}


}while(ngrid_i<ngrid);

free(K1);
free(K2);
free(K3);
free(grid);

}

