#include <math.h>

void ngc_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	
	int  xindex=(int)(ngrid*s_x/(L2-L1));
	int  yindex=(int)(ngrid*s_y/(L2-L1));
	int  zindex=(int)(ngrid*s_z/(L2-L1));
        long int index;
	double sx;
	double sy;
	double sz;

	sx=(xindex+0.5-s_x*ngrid*1./(L2-L1));
	sy=(yindex+0.5-s_y*ngrid*1./(L2-L1));
	sz=(zindex+0.5-s_z*ngrid*1./(L2-L1));

	if( fabs(sx)<=1./2. && fabs(sy)<=1./2. && fabs(sz)<=1./2.)//only contributing to the adjacent grids
	{
                        index=(pow(ngrid,2)*xindex+ngrid*yindex+zindex);
		        delta[index]=delta[index]+weight;
	}


}

void cic_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	

	int i,j,l;
	 int  xindex=(int)(ngrid*s_x/(L2-L1));
	  int  yindex=(int)(ngrid*s_y/(L2-L1));
	   int  zindex=(int)(ngrid*s_z/(L2-L1));

	    int  xindex2;
		 int  yindex2;
		  int  zindex2;

		  double sx;
		  double sy;
		  double sz;
                  long int index;
		  for(i=-1;i<=1;i++)
		  {
			  xindex2=xindex+i;
			  if(xindex2<0){xindex2=xindex2+ngrid;}
			  if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
			  sx=(xindex+i+0.5-s_x*ngrid*1./(L2-L1));

			  for(j=-1;j<=1;j++)
			  {
				  yindex2=yindex+j;
				  if(yindex2<0){yindex2=yindex2+ngrid;}
				  if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
				  sy=(yindex+j+0.5-s_y*ngrid*1./(L2-L1));


				  for(l=-1;l<=1;l++)
				  {
					  zindex2=zindex+l;
					  if(zindex2<0){zindex2=zindex2+ngrid;}
					  if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
					  sz=(zindex+l+0.5-s_z*ngrid*1./(L2-L1));


					  if( fabs(sx)<=1 && fabs(sy)<=1 && fabs(sz)<=1)//only contributing to the adjacent grids
					  {
                                                          index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
                                                          delta[index]=delta[index]+(1.-fabs(sx))*(1.-fabs(sy))*(1.-fabs(sz))*weight;

					  }

				  }
			  }
		  }

}
void tsc_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	

	int i,j,l;
	 int  xindex=(int)(ngrid*s_x/(L2-L1));
	  int  yindex=(int)(ngrid*s_y/(L2-L1));
	   int  zindex=(int)(ngrid*s_z/(L2-L1));

                  int  xindex2;
	          int  yindex2;
		  int  zindex2;
                  long int index;
		  double sx;
		  double sy;
		  double sz;

		  double factor_x,factor_y,factor_z;
		  factor_x=0;
		  factor_y=0;
		  factor_z=0;
		  int number_of_cells_explored=1;

		  for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
		  {
			  xindex2=xindex+i;
			  if(xindex2<0){xindex2=xindex2+ngrid;}
			  if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
			  sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));

			  factor_x=0;
			  if( fabs(sx)<1./2.)
			  {
				  factor_x=3./4.-sx*sx;
			  }
			  if(fabs(sx)>=1./2. && fabs(sx)<3./2.)
			  {
				  factor_x=1./2.*pow(3./2.-fabs(sx),2);
			  }


			  for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
			  {
				  yindex2=yindex+j;
				  if(yindex2<0){yindex2=yindex2+ngrid;}
				  if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
				  sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
				  factor_y=0;
				  if( fabs(sy)<1./2.)
				  {
					  factor_y=3./4.-sy*sy;
				  }
				  if(fabs(sy)>=1./2. && fabs(sy)<3./2.)
				  {
					  factor_y=1./2.*pow(3./2.-fabs(sy),2);
				  }

				  for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
				  {
					  zindex2=zindex+l;
					  if(zindex2<0){zindex2=zindex2+ngrid;}
					  if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
					  sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
					  factor_z=0;
					  if( fabs(sz)<1./2.)
					  {
						  factor_z=3./4.-sz*sz;
					  }
					  if(fabs(sz)>=1./2. && fabs(sz)<3./2.)
					  {
						  factor_z=1./2.*pow(3./2.-fabs(sz),2);
					  }

                                          index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);

					  delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

				  }
			  }
		  }

}

void pcs_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	

	int i,j,l;
	 int  xindex=(int)(ngrid*s_x/(L2-L1));
	  int  yindex=(int)(ngrid*s_y/(L2-L1));
	   int  zindex=(int)(ngrid*s_z/(L2-L1));
           long int index;
	    int  xindex2;
		 int  yindex2;
		  int  zindex2;

		  double sx;
		  double sy;
		  double sz;

		  double factor_x,factor_y,factor_z;
		  factor_x=0;
		  factor_y=0;
		  factor_z=0;
		  int number_of_cells_explored=2;

		  for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
		  {
			  xindex2=xindex+i;
			  if(xindex2<0){xindex2=xindex2+ngrid;}
			  if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
			  sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));

			  factor_x=0;
			  if( fabs(sx)<1)//only contributing to the adjacent grids
			  {
				          factor_x=2./3.-sx*sx+1./2.*pow(fabs(sx),3);
			  }
			  if( fabs(sx)>=1 && fabs(sx)<2.)//only contributing to the adjacent grids
			  {
				          factor_x=4./3.-2.*fabs(sx)+sx*sx-1./6.*pow(fabs(sx),3);
			  }

			  for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
			  {
				  yindex2=yindex+j;
				  if(yindex2<0){yindex2=yindex2+ngrid;}
				  if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
				  sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
				  factor_y=0;
				  if( fabs(sy)<1)//only contributing to the adjacent grids
				  {
					          factor_y=2./3.-sy*sy+1./2.*pow(fabs(sy),3);
				  }
				  if( fabs(sy)>=1 && fabs(sy)<2.)//only contributing to the adjacent grids
				  {
					          factor_y=4./3.-2.*fabs(sy)+sy*sy-1./6.*pow(fabs(sy),3);
				  }


				  for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
       			{
				zindex2=zindex+l;
				if(zindex2<0){zindex2=zindex2+ngrid;}
				if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
				sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
				factor_z=0;
				if( fabs(sz)<1)//only contributing to the adjacent grids
				{
					        factor_z=2./3.-sz*sz+1./2.*pow(fabs(sz),3);
				}
				if( fabs(sz)>=1 && fabs(sz)<2.)//only contributing to the adjacent grids
				{
					        factor_z=4./3.-2.*fabs(sz)+sz*sz-1./6.*pow(fabs(sz),3);
				}
                                index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
				delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

			}
			  }
		  }

}

void pq4s_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	
	int i,j,l;
	 int  xindex=(int)(ngrid*s_x/(L2-L1));
	  int  yindex=(int)(ngrid*s_y/(L2-L1));
	   int  zindex=(int)(ngrid*s_z/(L2-L1));
long int index;
	    int  xindex2;
		 int  yindex2;
		  int  zindex2;

		  double sx;
		  double sy;
		  double sz;

		  double factor_x,factor_y,factor_z;
		  factor_x=0;
		  factor_y=0;
		  factor_z=0;
		  int number_of_cells_explored=2;

		  for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
		  {
			  xindex2=xindex+i;
			  if(xindex2<0){xindex2=xindex2+ngrid;}
			  if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
			  sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));

			  factor_x=0;
			  if( fabs(sx)<1./2.)//only contributing to the adjacent grids
			  {
				          factor_x=115./192.-sx*sx*5./8.+sx*sx*sx*sx/4.;
			  }
			  if( fabs(sx)>=1./2. && fabs(sx)<3./2.)//only contributing to the adjacent grids
			  {
				          factor_x=55./96.+5./24.*fabs(sx)-sx*sx*5./4.+5./6.*pow(fabs(sx),3)-pow(sx,4)/6.;
			  }
			  if( fabs(sx)>=3./2. && fabs(sx)<5./2.)//only contributing to the adjacent grids
			  {
				          factor_x=625./384.-125./48.*fabs(sx)+25./16.*sx*sx-5./12.*pow(fabs(sx),3)+pow(sx,4)/24.;
			  }


			  for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
			  {
				  yindex2=yindex+j;
				  if(yindex2<0){yindex2=yindex2+ngrid;}
				  if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
				  sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
				  factor_y=0;
				  if( fabs(sy)<1./2.)//only contributing to the adjacent grids
				  {
					          factor_y=115./192.-sy*sy*5./8.+sy*sy*sy*sy/4.;
				  }
				  if( fabs(sy)>=1./2. && fabs(sy)<3./2.)//only contributing to the adjacent grids
				  {
					          factor_y=55./96.+5./24.*fabs(sy)-sy*sy*5./4.+5./6.*pow(fabs(sy),3)-pow(sy,4)/6.;
				  }
				  if( fabs(sy)>=3./2. && fabs(sy)<5./2.)//only contributing to the adjacent grids
				  {
					          factor_y=625./384.-125./48.*fabs(sy)+25./16.*sy*sy-5./12.*pow(fabs(sy),3)+pow(sy,4)/24.;
				  }
				  for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
				  {
					  zindex2=zindex+l;
					  if(zindex2<0){zindex2=zindex2+ngrid;}
					  if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
					  sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
					  factor_z=0;
					  if( fabs(sz)<1./2.)//only contributing to the adjacent grids
					  {
						          factor_z=115./192.-sz*sz*5./8.+sz*sz*sz*sz/4.;
					  }
					  if( fabs(sz)>=1./2. && fabs(sz)<3./2.)//only contributing to the adjacent grids
					  {
						          factor_z=55./96.+5./24.*fabs(sz)-sz*sz*5./4.+5./6.*pow(fabs(sz),3)-pow(sz,4)/6.;
					  }
					  if( fabs(sz)>=3./2. && fabs(sz)<5./2.)//only contributing to the adjacent grids
					  {
						          factor_z=625./384.-125./48.*fabs(sz)+25./16.*sz*sz-5./12.*pow(fabs(sz),3)+pow(sz,4)/24.;
					  }
                                          index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
					  delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

				  }
			  }
		  }

}

void pq5s_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);
	
	int i,j,l;
	 int  xindex=(int)(ngrid*s_x/(L2-L1));
	  int  yindex=(int)(ngrid*s_y/(L2-L1));
	   int  zindex=(int)(ngrid*s_z/(L2-L1));

	    int  xindex2;
		 int  yindex2;
		  int  zindex2;
                  long int index;
		  double sx;
		  double sy;
		  double sz;

		  double factor_x,factor_y,factor_z;
		  factor_x=0;
		  factor_y=0;
		  factor_z=0;
		  int number_of_cells_explored=3;

		  for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
		  {
			  xindex2=xindex+i;
			  if(xindex2<0){xindex2=xindex2+ngrid;}
			  if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
			  sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));

			  factor_x=0;
			  if( fabs(sx)<1.)//only contributing to the adjacent grids
			  {
				          factor_x=11./20.-0.5*sx*sx+0.25*sx*sx*sx*sx-1./12.*pow(fabs(sx),5);
			  }
			  if( fabs(sx)>=1. && fabs(sx)<2.)//only contributing to the adjacent grids
			  {
				          factor_x=17./40.+5./8.*fabs(sx)-7./4.*sx*sx+5./4.*pow(fabs(sx),3)-3./8.*pow(sx,4)+1./24.*pow(fabs(sx),5);
			  }
			  if( fabs(sx)>=2. && fabs(sx)<3.)//only contributing to the adjacent grids
			  {
				          factor_x=243./120.-81./24.*fabs(sx)+9./4.*pow(sx,2)-3./4.*pow(fabs(sx),3)+1./8.*pow(sx,4)-1./120.*pow(fabs(sx),5);
			  }


			  for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
			  {
				  yindex2=yindex+j;
				  if(yindex2<0){yindex2=yindex2+ngrid;}
				  if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
				  sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
				  factor_y=0;
				  if( fabs(sy)<1.)//only contributing to the adjacent grids
				  {
					          factor_y=11./20.-0.5*sy*sy+0.25*sy*sy*sy*sy-1./12.*pow(fabs(sy),5);
				  }
				  if( fabs(sy)>=1. && fabs(sy)<2.)//only contributing to the adjacent grids
				  {
					          factor_y=17./40.+5./8.*fabs(sy)-7./4.*sy*sy+5./4.*pow(fabs(sy),3)-3./8.*pow(sy,4)+1./24.*pow(fabs(sy),5);
				  }
				  if( fabs(sy)>=2. && fabs(sy)<3.)//only contributing to the adjacent grids
				  {
					          factor_y=243./120.-81./24.*fabs(sy)+9./4.*pow(sy,2)-3./4.*pow(fabs(sy),3)+1./8.*pow(sy,4)-1./120.*pow(fabs(sy),5);
				  }


				  for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
				  {
					  zindex2=zindex+l;
					  if(zindex2<0){zindex2=zindex2+ngrid;}
					  if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
					  sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
					  factor_z=0;
					  if( fabs(sz)<1.)//only contributing to the adjacent grids
					  {
						          factor_z=11./20.-0.5*sz*sz+0.25*sz*sz*sz*sz-1./12.*pow(fabs(sz),5);
					  }
					  if( fabs(sz)>=1. && fabs(sz)<2.)//only contributing to the adjacent grids
					  {
						          factor_z=17./40.+5./8.*fabs(sz)-7./4.*sz*sz+5./4.*pow(fabs(sz),3)-3./8.*pow(sz,4)+1./24.*pow(fabs(sz),5);
					  }
					  if( fabs(sz)>=2. && fabs(sz)<3.)//only contributing to the adjacent grids
					  {
						          factor_z=243./120.-81./24.*fabs(sz)+9./4.*pow(sz,2)-3./4.*pow(fabs(sz),3)+1./8.*pow(sz,4)-1./120.*pow(fabs(sz),5);
					  }
                                          index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
					  delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

				  }
			  }
		  }

}

