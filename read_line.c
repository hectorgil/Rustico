#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structures.h"
void get_line(FILE *f, double params[],int type)//skycut option
{

double RA,dec,redshift,wcp,wnoz,wcol,wsys,wfkp,n_z,dist;
int veto,wcp_int,wnoz_int,wcol_int;
char wfkp_cmass_c[200],wfkp_eboss_c[200];
double wfkp_cmass,wfkp_eboss;
wcol=1;//collision weight
wcp=1;
wnoz=1;
wfkp=1;//FKP weight
n_z=1;//number density of galaxies
veto=1;//whether this galaxy needs to be vetoed (veto=0 for vetoed galaxy, otherwise and by default 1)
wsys=1;//systematic/imaging weight.
redshift=-1;//default redshift off (if both distances and redshifts are provided, redshift is taken as the default input)
dist=-1;//default distances off (if both distances and redshifts are provided, redshift is taken as the default input)
char in_eboss_foot[20],iscmass[200];
double Pi=(4.*atan(1.));
//if some of these inputs doesn't exist in your survey (for both data or random catalogue) do not assigne a value (by default is 1). 

//Total weight which will be applied to each galaxy: wsys*veto*wfkp*wcol.
//Note that for some surveys, such as BOSS DR12, the catalogues provide the close-pair weight (wcp) and the redshift failure weight (wnoz), and the total collision weight is formed out of these weight as: wcol=(wcp+wnoz-1). In these and similar cases you will need to build the collision weight from the individual pieces, just below the fscanf line. 

if(type==0){//data

//fscanf(f,"%lf %lf %lf %*f %lf %*f %d %d\n",&RA,&dec,&redshift, &n_z, &veto, &wcol_int);wsys=1.0;//For BOSS LRG patchy mocks
//wfkp=1./(1+10000.*n_z);
//wcol=wcol_int*1.;

fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wcp,&wnoz,&wsys,&n_z);wcol=wcp+wnoz-1;//For BOSS LRGs data

//fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wcp,&wnoz,&n_z,&wsys);//for eBOSS LRG mocks
//if(wnoz==1 || wnoz==2 || wnoz>2){wcol=wcp+wnoz-1;}
//else{wcol=wcp*wnoz;}


//fscanf(f,"%lf %lf %lf %s %lf %lf %lf %lf %*s %s %*s %s %*s %s %lf\n",&RA,&dec,&redshift,wfkp_eboss_c,&wsys,&wcp,&wnoz,&n_z,iscmass,in_eboss_foot,wfkp_cmass_c,&wfkp);// For eBOSS LRG data
//if(strcmp(iscmass,"T") == 0){wcol=wcp+wnoz-1;}
//else{wcol=wcp*wnoz;}

//fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %*s\n",&RA,&dec,&redshift,&wfkp,&wsys,&wcp,&wnoz,&n_z);wcol=wcp*wnoz;//For eBOSS QSO data

//fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wnoz,&wcp,&n_z,&wsys);wcol=wcp*wnoz;//For eBOSS QSO mocks

}

if(type==1){//random

//  fscanf(f,"%lf %lf %lf %lf %*f %d %d\n",&RA, &dec, &redshift, &n_z,&veto,&wcol_int);//For BOSS LRG patchy mocks
//  wfkp=wcol_int*1./(1.+10000*n_z);
//  wcol=wcol_int*1.;

fscanf(f,"%lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&n_z);//For BOSS LRGs data

//fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wcp,&wnoz,&n_z,&wsys);//For eBOSS LRG mocks
//if(wnoz==1 || wnoz==2 || wnoz>2){wcol=wcp+wnoz-1;}
//else{wcol=wcp*wnoz;}


//fscanf(f,"%lf %lf %lf %s %lf %lf %lf %lf %*s %s %s %s %*s %lf\n",&RA,&dec,&redshift,wfkp_eboss_c,&wsys,&wcp,&wnoz,&n_z,in_eboss_foot,iscmass,wfkp_cmass_c,&wfkp);// For eBOSS LRG data
//if(strcmp(iscmass,"T") == 0){wcol=wcp+wnoz-1;}
//else{wcol=wcp*wnoz;}

//fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wsys,&wcp,&wnoz,&n_z);wcol=wcp*wnoz;//For eBOSS QSO data

//fscanf(f,"%lf %lf %*s %*s %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&wfkp,&wnoz,&wcp,&wsys,&n_z,&redshift);wcol=wcp*wnoz;//For eBOSS QSO mocks

}

params[0]=RA;
params[1]=dec;
params[2]=redshift;
params[3]=wfkp;
params[4]=n_z;
params[5]=wcol;
params[6]=wsys;
params[7]=veto;
params[8]=dist;
}

void get_line_periodic(FILE *f, double params[],int type)//for periodic box
{
double x,y,z,vx,vy,vz,weight;
double scale_factor;
int veto;
weight=1;
veto=1;
//strictly necessary only need to provide x,y,z;

//data
if(type==0){

//fscanf(f,"%lf %lf %lf %*f %*f %*f %*f %*f %*f\n",&x,&y,&z);
fscanf(f,"%lf %lf %lf %*f %*f %lf\n",&x,&y,&z,&vz);

}

//randoms
if(type==1){
fscanf(f,"%lf %lf %lf\n",&x,&y,&z);

}

     scale_factor=0.66667;//z=0.5
     z=z+vz*1.0/(100.*scale_factor*sqrt(0.286*pow(scale_factor,-3)+1.-0.286));//z-space distortion correction (assuming z direction)          //z=z+vz/(100.*sqrt(0.286*pow(scale_factor,-3)+1.-0.286));


        if(z>=2600.){z=z-2600.;}
        if(z<0){z=z+2600.;}
        

params[0]=x;
params[1]=y;
params[2]=z;
params[3]=weight;
params[7]=veto;
}
