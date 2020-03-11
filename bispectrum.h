
double Minimum(double a, double b);

int Sign(double x);

void write_triangles_bin(int ngrid,double *Triangles_k, long int *Weight_k, double k_input, double Deltakbis, double  kf, long int *L,long int *IJK);

double Number_triangles_exact(double K_eff1,double K_eff2,double K_eff3,double Deltakbis);

double Number_triangles_exact_equi(double k,double Dk);

void FFT_k(int ngrid, double *out_k1, double k_input, double Deltakbis, double kf, long int *L,long int *IJK,long int *Kz,double *deltak_re, double *deltak_im, double keff_out[], long int neff_out[]);

void FFT_kt(int ngrid, double *out_k1,double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK, long int *Kz,double *deltak_re, double *deltak_im, double *keff_out, long int *neff_out);

void FFT_t(int ngrid, double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK);

void select_triangles(double kmin,double kmax, double Deltak, double L1, double L2, double *k1, double *k2, double *k3, int *which_grid, long int N, int ngrid, char *do_multigrid, char *triangle_shapes);

long int count_triangles(double kmin,double kmax, double Deltak, char *triangle_shapes);

void loop_bispectrum_skycut_caller(double kmin,double kmax, int Ninterlacing,  double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, int ngrid, double P_shot_noise, double Deltakbis, double I33,double I22, double IN, double Bsn,  double alpha, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *do_multigrid, char *triangle_shapes,char *do_bispectrum2);

void loop_bispectrum_periodic_for_gadget_caller(double kmin,double kmax, int Ninterlacing, double L1, double L2, int ngrid, double Deltakbis, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *name_data_in,int gadget_files, char *do_multigrid, char *triangle_shapes,char *RSD, char *do_bispectrum2);

void bispectrum_calculator(double *K1,double **K2, double ***K3, int Nk1, int *Nk2, int **Nk3, double *deltak_re, double *deltak_im, double *deltak2_re, double *deltak2_im,int Ninterlacing, double kmin,double kmax,double L1,double L2, int ngrid, double Deltakbis, double I33,double I22, double IN, double Bsn, double Psn, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *triangles_num, char *write_triangles, char *triangles_id, long int *input_bin,char *do_bispectrum2);
