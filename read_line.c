#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void get_line(FILE *f, double params[],int type)//skycut option
{

double RA,dec,redshift,wcp,wnoz,wcol,wsys,wfkp,n_z;
int veto,wcp_int,wnoz_int,wcol_int;
wcol=1;//collision weight
wcp=1;
wnoz=1;
wfkp=1;//FKP weight
n_z=1;//number density of galaxies
veto=1;//whether this galaxy needs to be vetoed (veto=0 for vetoed galaxy, otherwise and by default 1)
wsys=1;//systematic/imaging weight.
//if some of these inputs doesn't exist in your survey (for both data or random catalogue) do not assigne a value (by default is 1). 

//Total weight which will be applied to each galaxy: wsys*veto*wfkp*wcol.
//Note that for some surveys, such as BOSS DR12, the catalogues provide the close-pair weight (wcp) and the redshift failure weight (wnoz), and the total collision weight is formed out of these weight as: wcol=(wcp+wnoz-1). In these and similar cases you will need to build the collision weight from the individual pieces, just below the fscanf line. 

if(type==0){//data
fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&wcp,&wnoz,&wsys,&n_z);wcol=wcp+wnoz-1;
//wfkp=1;
//wsys=1;
//wcp=1;
//wnoz=1;
//wcol=1;
}

if(type==1){//random
fscanf(f,"%lf %lf %lf %lf %lf\n",&RA,&dec,&redshift,&wfkp,&n_z);
//wfkp=1;
//wsys=1;
//wcp=1;
//wnoz=1;
//wcol=1;
}

params[0]=RA;
params[1]=dec;
params[2]=redshift;
params[3]=wfkp;
params[4]=n_z;
params[5]=wcol;
params[6]=wsys;
params[7]=veto;

}


void get_line_periodic(FILE *f, double params[])//for periodic box
{
double x,y,z,vx,vy,vz,weight;
double scale_factor;
weight=1;

//strictly necessary only need to provide x,y,z;

fscanf(f,"%lf %lf %lf %lf %lf %lf %*s %*d %*d %*d %*d %*f %*f\n",&x,&y,&z,&vx,&vy,&vz);//tinker

        //In case you need to introduce RSD-shifts like, do it here. You need to perform the shifts 'by hand'
        //scale_factor=1.0;//z=0
        //scale_factor=2./3.;//z=0.5
        //scale_factor=0.5;//z=1.0
        //scale_factor=0.4;//z=1.5
        //scale_factor=1./3.;//z=4
        scale_factor=1./1.695; //z=0.695
        //scale_factor=1./1.865;//z=0.865
//        z=z+vz*1.0/(100.*scale_factor*sqrt(0.26479*pow(scale_factor,-3)+1.-0.26479));//z-space distortion correction (assuming z direction)
        z=z+vz/(100.*sqrt(0.26479*pow(scale_factor,-3)+1.-0.26479));
        if(z>=1000.){z=z-1000.;}
        if(z<0){z=z+1000.;}
        

params[0]=x;
params[1]=y;
params[2]=z;
params[3]=weight;
}
