double Leg1(double x);

double Leg2(double x);

double Leg3(double x);

double Leg4(double x);

void write_power_spectrum_skyscuts(double kmin,double kmax,double kx[], double ky[], double kz[], double P0r[], double P0i[],double P1r[], double P1i[], double P2r[], double P2i[], double P3r[], double P3i[], double P4r[], double P4i[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type, char *do_anisotropy, char *do_odd_multipoles, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type,char *type,int N_interlacing);

void write_power_spectrum_periodic(double kmin, double kmax, double deltak_re[], double deltak_im[], double bin_ps, long int  ngrid, double L1, double L2, int Ninter, char *name_ps_out, double P_shot_noise, char *binning_type, char *do_odd_multipoles, char *do_anisotropy);

void write_power_spectrum_periodic2D(double kmin, double kmax, double deltak_re[], double deltak_im[], double Deltak, int mubin, int  ngrid, double L1, double L2, int Ninterlacing, char *name_ps_out, double P_shot_noise, char *binning_type, char *file_for_mu);
