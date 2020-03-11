#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void prepare_pointer_to_delta(int ngrid)
{
long int L_in, IJK_in;
int K_in;
FILE *f,*g;
long int l,L;
int i,j,k,I,J,K;
char name[200];
long int ngridtot;
long int count;
sprintf(name,"temp.dat");
f=fopen(name,"w");

ngridtot=pow(ngrid,3);

count=0;
for(l=0;l<ngridtot;l++)
{
i=(int)(l/(ngrid*ngrid*1.));
j=(int)( (l-i*ngrid*ngrid)/(ngrid*1.));
k=l-i*ngrid*ngrid-j*ngrid;

if(i<=ngrid/2){I=i;}
else{I=-(ngrid-i);}

if(j<=ngrid/2){J=j;}
else{J=-(ngrid-j);}

if(k<=ngrid/2){K=k;}
else{K=-(ngrid-k);}

L=I*I+J*J+K*K;

if(L<=ngrid*ngrid/4.)
{
fprintf(f,"%ld %ld %d\n",l,L,K);
count++;
}

}
fclose(f);

system("sort -k2 -n temp.dat > temp2.dat");
system("rm temp.dat");

sprintf(name,"wisdom%d.dat",ngrid);

f=fopen("temp2.dat","r");
g=fopen(name,"w");
for(l=0;l<count;l++) 
{ 
fscanf(f,"%ld %ld %d\n",&L_in,&IJK_in,&K_in);
fprintf(g,"%ld %ld %d\n",L_in,IJK_in,K_in);

} 
fclose(f);
fclose(g);
system("rm temp2.dat");


}


/*
void perpare_pointer_to_delta(long int *L_in, long int *IJK_in, int ngrid)
{
//L and IJK must be of size ngrid*ngrid*ngrid?
long int l,lread,L,c;
int i,j,k,I,J,K;


for(l=0;l<ngrid*ngrid*ngrid;l++)
{
i=(int)(l/(ngrid*ngrid*1.));
j=(int)( (l-i*ngrid*ngrid)/(ngrid*1.));
k=l-i*ngrid*ngrid-j*ngrid;

if(i<=ngrid/2){I=i;}
else{I=-(ngrid-i);}

if(j<=ngrid/2){J=j;}
else{J=-(ngrid-j);}

if(k<=ngrid/2){K=k;}
else{K=-(ngrid-k);}

L=I*I+J*J+K*K;

if(l==0)
{
L_in[l]=l;
IJK_in[l]=L;
}
else//leff is at least 1
{
for(lread=0;lread<l;lread++)
{
if(IJK_in[lread]>L){break;}
}
//end without breaking => l=lread
if(lread==l)
{
L_in[l]=l;
IJK_in[l]=L;
}
else//end breaking at lread<l
{
//displace all elements of L and IJK inversely from l-1 to lread-1 by one unit
for(c=l-1;c>lread-1;c--)
{
//L_in[c+1]=L_in[c];
//IJK_in[c+1]=IJK_in[c];
*(L_in+c+1) = *(L_in+c);
*(IJK_in+c+1) = *(IJK_in+c);

}
//write L and IJK leff
L_in[lread]=l;
IJK_in[lread]=L;
}
}
}

//end of function
}
*/
