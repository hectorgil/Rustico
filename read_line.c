//You will need to modify certain line of this file to match your code structure. This operation must be done BEFORE compiling

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

if(type==0){//data. Skycut option

//modify the line below accordingly to your file structure.
/* for eg. */ fscanf(f,"%lf %lf %lf %*f\n",&RA,&dec,&redshift);// first column: Right Ascencion (in rad); second declination (in rad); third redshift; forth unsued. 

//uncoment for changing deg to rad. 
//RA=RA*180./Pi;
//dec=dec*180/Pi;
//dec=dec-90.;//theta to dec


}

if(type==1){//random Skycut option

//modify the line below accordingly to your file structure.
fscanf(f,"%lf %lf %lf %*f\n",&RA,&dec,&redshift);

//uncoment for changing deg to rad.
//RA=RA*180./Pi;
//dec=dec*180/Pi;
//dec=dec-90.;//theta to dec


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
long int i;
//strictly necessary only need to provide x,y,z;

//data
if(type==0){
//modify the line below accordingly to your file structure.
/* for eg. */ fscanf(f,"%lf %lf %lf %lf %lf %lf\n",&x,&y,&z,&vx,&vy,&vz);//x-position, y-position; z-position; x-velocity, y-velocity; z-velocity
}

//randoms
if(type==1){
//modify the line below accordingly to your file structure.
fscanf(f,"%lf %lf %lf\n",&x,&y,&z);

}


//This is just an example of how to apply RSD along the line of sight. These formula are code-dependent! Make sure you input the right values and formula if you want to apply RSD. 
/*
     scale_factor=0.66667;//z=0.5
     z=z+vz*1.0/(100.*scale_factor*sqrt(0.286*pow(scale_factor,-3)+1.-0.286));//z-space distortion correction (assuming z direction)     

	    //This ensures the particles re-enter the box. In this case of size 2600 (hardcoded value). 
        if(z>=2600.){z=z-2600.;}
        if(z<0){z=z+2600.;}
*/      

params[0]=x;
params[1]=y;
params[2]=z;
params[3]=weight;
params[7]=veto;
}
