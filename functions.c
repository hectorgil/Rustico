#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "structures.h"

double Leg1(double x)
{
double f=x;

return f;
}

double Leg3(double x)
{
double f=0.5*(5.*x*x-3.)*x;

return f;
}


double Leg2(double x)
{
double f;
f=0.5*(3.*x*x-1.);
return f;
}

double Leg4(double x)
{
double f;
f=1./8.*(35.*pow(x,4)-30.*pow(x,2)+3.);
return f;
}

double Leg6(double x)
{
double f;
f=1./16.*(231.*pow(x,6)-315.*pow(x,4)+105*pow(x,2)-5.);
return f;
}


double Leg8(double x)
{
double f;
f=1./128.*(6435.*pow(x,8)-12012.*pow(x,6)+6930.*pow(x,4)-1260.*pow(x,2)+35.);
return f;
}

void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
     f_params params_function = *(f_params *) fdata;//cast from void to f_params
     double omega_m= params_function.OMEGA_M;
     double omega_L= params_function.OMEGA_L;
     double omega_k=1-omega_m-omega_L;
     double z=x[0];
//     double c=299792.458;//speed of light
//     fval[0]=c/(100.*sqrt(omega_m*pow(1+z,3)+1.-omega_m));

     fval[0]=1./(sqrt(omega_m*pow(1+z,3)+omega_L+omega_k*pow(1+z,2)));


}

void randstring(char *name,long int seed) {

sprintf(name,"temp_file_%ld.dat",seed);

}

void shuffle_lines(long int npar_used,long int *L,long int seed)
{
srand48(seed);
int fraction_ramdomized_part=90;//fraction (over 100) of randomized particles. The rest are kept in the same order filling the remaining free slots. 
FILE *f;
//char *name,*name2;
char name[2000],name2[2000];
char name_command[2000];
long int i,j,interval,l,c,b;
long int *Lminus;
long int random;
double randomd;

int n_parallel;//number of threads
int tid;//thread id
int *interval_thread;//
int *exit_condition;
Lminus=  (long int*) calloc(npar_used, sizeof(long int)); //for(i=0;i<npar_used;i++)L[i]=-1;

       #pragma omp parallel for private(i,tid) shared(n_parallel)
        for(i=0;i<10000;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){n_parallel=omp_get_num_threads();}
        }
        i=fftw_init_threads();
        fftw_plan_with_nthreads(n_parallel);

interval_thread=  (int*) calloc(n_parallel, sizeof(int));
exit_condition=  (int*) calloc(n_parallel, sizeof(int));

//name=  (char*) calloc(20, sizeof(char));
//name2=  (char*) calloc(20, sizeof(char));

 randstring(name,seed);
i=0;
do
{
randomd=drand48();
random=(long int)(randomd*npar_used);
if(random>=npar_used || random<0){random=0;}
L[i]=random;
interval=1;

if(random==i){interval=0;}
else{

      if(i>0)
      {

          if(i>npar_used*1/100){//only avoid it in the 1st 1% (unlikley to be needed)

                for(tid=0;tid<n_parallel;tid++){interval_thread[tid]=1;exit_condition[tid]=0;}//set all intervals to 1
                #pragma omp parallel for private(j) shared(i,interval_thread,L)
                for(j=0;j<i;j++)//this can be parallelized
                { tid=omp_get_thread_num();
                   if(L[i]==L[j]){interval_thread[tid]=0; for(l=0;l<n_parallel;l++){exit_condition[l]=i;}}//break

                 j=j+exit_condition[tid];
                }

                   for(l=1;l<n_parallel;l++)
                   {
                    interval_thread[0]=interval_thread[0]+interval_thread[l];
                   }
                   if(interval_thread[0]==n_parallel){interval=1;}//all interval_thread to 1 (never exited by L[i]=L[j], we can take this realization
                   else{interval=0;}//at least one L[i]=L[j] was found, we need another realization of drand48

          }else{

                for(j=0;j<i;j++)
                {
                   if(L[i]==L[j]){interval=0;break;}
                }
          }


      }
}

//printf("%d/%d -> %d/%d\n",i,npar_used,L[i],npar_used);
i=i+interval;
}while(i<npar_used*fraction_ramdomized_part/100);//do this to 90% of particles, last 10% assigned in order among the remaining spots

f=fopen(name,"w");//write the L[i] relation
for(l=0;l<i;l++)
{
fprintf(f,"%ld %ld\n",l,L[l]);
}
fclose(f);

randstring(name2,seed*10);

//printf("Temporal files %s and %s (%ld)\n",name,name2,seed);

sprintf(name_command,"sort -k2 -n %s > %s",name,name2);
system(name_command);

sprintf(name_command,"rm %s",name);
system(name_command);

f=fopen(name2,"r");
for(l=0;l<i;l++)
{
fscanf(f,"%*d %ld\n", &Lminus[l]);
}
fclose(f);
//exit(0);
c=0;
b=0;

//for(l=i;l<npar_used;l++){//orignal order
for(l=npar_used-1;l>=i;l--){//reverse order

      for(j=c;j<npar_used;j++){
            
          if(Lminus[j-b]!=j){break;}
      }

 L[ l ]=j;
 c=j+1;
 b++;
}

sprintf(name_command,"rm %s",name2);
system(name_command);

free(Lminus);
//free(name);
//free(name2);
free(interval_thread);
free(exit_condition);
}

void reshuffle(double radata[],double decdata[],double wcoldata[],double wsysdata[], long int npar_used,long int seed)//this re-shuffles the indices on angular-data wrt the unchanged radial-data in case of both-shuffling option.
{
//srand48(seed);

long int i,j;
double *a1;
double *a2;
double *a3;
double *a4;
long int *L;


L=  (long int*) calloc(npar_used, sizeof(long int)); //for(i=0;i<npar_used;i++)L[i]=-1;
a1=  (double*) calloc(npar_used, sizeof(double));
a2=  (double*) calloc(npar_used, sizeof(double));
a3=  (double*) calloc(npar_used, sizeof(double));
a4=  (double*) calloc(npar_used, sizeof(double));

printf("Re-shuffling angles...\t");

shuffle_lines(npar_used,L,seed);

/*
int n_parallel;//number of threads
int tid;//thread id
int *interval_thread;//
int *exit_condition;
Lminus=  (long int*) calloc(npar_used, sizeof(long int)); //for(i=0;i<npar_used;i++)L[i]=-1;

//
        #pragma omp parallel for private(i,tid) shared(n_parallel)
        for(i=0;i<10000;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){n_parallel=omp_get_num_threads();}
        }
        i=fftw_init_threads();
        fftw_plan_with_nthreads(n_parallel);

interval_thread=  (int*) calloc(n_parallel, sizeof(int));
exit_condition=  (int*) calloc(n_parallel, sizeof(int));

name=  (char*) calloc(20, sizeof(char));
name2=  (char*) calloc(20, sizeof(char));

 randstring(name,seed);


i=0;
do
{
randomd=drand48();
random=(long int)(randomd*npar_used);
if(random>=npar_used || random<0){random=0;}
L[i]=random;
interval=1;

if(random==i){interval=0;}
else{

      if(i>0)
      {
 
          if(i>npar_used*1/100){//only avoid it in the 1st 1% (unlikley to be needed)

                for(tid=0;tid<n_parallel;tid++){interval_thread[tid]=1;exit_condition[tid]=0;}//set all intervals to 1
                #pragma omp parallel for private(j) shared(i,interval_thread,L)
                for(j=0;j<i;j++)//this can be parallelized
                { tid=omp_get_thread_num();
                   if(L[i]==L[j]){interval_thread[tid]=0; for(l=0;l<n_parallel;l++){exit_condition[l]=i;}}//break

                 j=j+exit_condition[tid];
                }

                   for(l=1;l<n_parallel;l++)
                   {
                    interval_thread[0]=interval_thread[0]+interval_thread[l];
                   }
                   if(interval_thread[0]==n_parallel){interval=1;}//all interval_thread to 1 (never exited by L[i]=L[j], we can take this realization
                   else{interval=0;}//at least one L[i]=L[j] was found, we need another realization of drand48

          }else{

                for(j=0;j<i;j++)
                {
                   if(L[i]==L[j]){interval=0;break;}
                }
          }
 
          
      }
}


i=i+interval;
}while(i<npar_used*90/100);//do this to 90% of particles, last 10% assigned in order among the remaining spots

f=fopen(name,"w");//write the L[i] relation
for(l=0;l<i;l++)
{
fprintf(f,"%d %d\n",l,L[l]);
}
fclose(f);

randstring(name2,random);

sprintf(name_command,"sort -k2 -n %s > %s",name,name2);
system(name_command);

sprintf(name_command,"rm %s",name);
system(name_command);


f=fopen(name2,"r");
for(l=0;l<i;l++)
{
fscanf(f,"%*d %d\n", &Lminus[l]);
}
fclose(f);

c=0;
b=0;
for(l=i;l<npar_used;l++)
{ 
      for(j=c;j<i;j++)
      {
          if(Lminus[j-b]!=j){break;}
      }   

 L[ l ]=j;
 c=j+1;
 b++;
}

sprintf(name_command,"rm %s",name2);
system(name_command);
*/
printf("Done");


for(i=0;i<npar_used;i++)
{
/*if(L[i]==i){printf("Error when reshuffling...%ld\n",i);exit(0);}*/
a1[L[i]]=radata[i];
a2[L[i]]=decdata[i];
a3[L[i]]=wcoldata[i];
a4[L[i]]=wsysdata[i];
}

for(i=0;i<npar_used;i++)
{
radata[i]=a1[i];
decdata[i]=a2[i];
wcoldata[i]=a3[i];
wsysdata[i]=a4[i];
}
printf(" and complete\n\n");

free(a1);
free(a2);
free(a3);
free(a4);
free(L);
//free(Lminus);
}
double P_interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N-1);
if(i==0)
{
//printf("Interpol %d %lf %lf %lf %lf %lf -> \t",i,k0,k[i],k[i+1],P[i],P[i+1]);
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
else
{
//printf("Interpol %d %lf %lf %lf %lf %lf- >\t",i,k0,k[i-1],k[i],P[i-1],P[i]);
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
//printf("%lf\n",P0);
return P0;
}


double Interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N-1);
if(i==0)
{
//printf("Interpol %d %lf %lf %lf %lf %lf -> \t",i,k0,k[i],k[i+1],P[i],P[i+1]);
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
//P0=pow(10.,P0);
}
else
{
//printf("Interpol %d %lf %lf %lf %lf %lf- >\t",i,k0,k[i-1],k[i],P[i-1],P[i]);
m=( (P[i]) - (P[i-1]) )/( (k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
//P0=pow(10.,P0);
}
//printf("%lf\n",P0);
return P0;
}

int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
int ch, number_of_lines = -1;

do 
{
	    ch = fgetc(myfile);
		    if(ch == '\n')
				        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0) 
    number_of_lines++;

	fclose(myfile);
	return number_of_lines;
}

void check_box_for_yamamoto(double *parameters, int ngrid)
{
double L1,L2;
int warning;
double vector;
int vector_i;
double epsilon,deltaepsilon,L1old,L2old;
int j;

L1old=parameters[0];
L2old=parameters[1];
L1=L1old;
L2=L2old;
deltaepsilon=(L2-L1)/ngrid*1.*0.001;//0.1 percent of gridcell
j=0;
warning=0;
do
{
warning=0;
epsilon=epsilon+j*deltaepsilon;

vector=-L1*ngrid/(L2-L1);
vector_i=(int)(vector);

if(vector==vector_i*1.){warning=-1;}

L1=L1+epsilon;
L2=L2+epsilon;
j++;
}while(warning==-1);

if(L1!=L1old){printf("Box displaced by %lf (iterated %d times) from (L1,L2)=(%lf, %lf) to (%lf, %lf) for better performance\n",epsilon,j,L1old,L2old,L1,L2);}
parameters[0]=L1;
parameters[1]=L2;


}

void check_box_for_yamamoto_old(double *parameters, int ngrid)
{
//printf("%lf %lf %d\n",parameters[0],parameters[1],ngrid);
int i,j;
int warning;
double epsilon,deltaepsilon,L1old,L2old;
double l;
double L1,L2;
L1old=parameters[0];
L2old=parameters[1];
L1=L1old;
L2=L2old;
j=0;
deltaepsilon=(L2-L1)/ngrid*1.*0.001;//0.1 percent of gridcell
//printf("delta=%lf\n",deltaepsilon);
do
{
warning=0;
//epsilon+=j*deltaepsilon;
epsilon=epsilon+j*deltaepsilon;
//printf("j=%d, epsilon=%lf\n",j,epsilon);
for(i=0;i<ngrid;i++)
{
l=L1+i*(L2-L1)/ngrid*1.;
//printf("i=%d, l=%lf\n",i,l);
if(l==0){warning=-1;break;}
}
//L1+=epsilon;
//L2+=epsilon;
L1=L1+epsilon;
L2=L2+epsilon;
j++;
}while(warning==-1);

if(L1!=L1old){printf("Box displaced by %lf (iterated %d times) from (L1,L2)=(%lf, %lf) to (%lf, %lf) for better performance\n",epsilon,j,L1old,L2old,L1,L2);}
parameters[0]=L1;
parameters[1]=L2;

}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);

}




void freeTokens2(double ***tokens, int N1, int *N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2[i];++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

void freeTokens3(double ***tokens, int N1, int N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

void freeTokens4(double ****tokens, int N1, int N2, int N3)
{
int i,j,l;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        for(l=0;l<N3;++l)
        {
            free(tokens[i][j][l]);
        }
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

void freeTokensInt(int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}


void freeTokensLInt(long int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}

void freeTokensInt2(int ***tokens,int N1,int *N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2[i];++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

void freeTokensInt3(int ***tokens,int N1,int N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

void freeTokensLInt3(long int ***tokens,int N1,int N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}

/*
void randstring(char *name) {

int random;
random=drand48();
sprintf(name,"temp_file_%d.dat",random);

}
*/
