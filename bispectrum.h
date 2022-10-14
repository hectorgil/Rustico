
double Minimum(double a, double b);

int Sign(double x);

void write_triangles_bin(int ngrid,double *Triangles_k, long int *Weight_k, double k_input, double Deltakbis, double  kf, long int *L,long int *IJK);

double Number_triangles_exact(double K_eff1,double K_eff2,double K_eff3,double Deltakbis);

double Number_triangles_exact_equi(double k,double Dk);

void FFT_k(int ngrid, double *out_k1, double k_input, double Deltakbis, double kf, long int *L,long int *IJK,long int *Kz,double *deltak_re, double *deltak_im, double keff_out[], long int neff_out[],int mode);

void FFT_kt(int ngrid, double *out_k1,double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK, long int *Kz,double *deltak_re, double *deltak_im, double *keff_out, long int *neff_out);

void FFT_t(int ngrid, double *out_k1_NT, double k_input, double Deltakbis, double kf, long int *L,long int *IJK);

void select_triangles(double kmin,double kmax, double Deltak, double L1, double L2, double *k1, double *k2, double *k3, int *which_grid, long int N, int ngrid, char *do_multigrid, char *triangle_shapes);

long int count_triangles(double kmin,double kmax, double Deltak, char *triangle_shapes);

void loop_bispectrum_skycut_caller(double kmin,double kmax, int Ninterlacing,  double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double Ndataw, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double *s_xB, double *s_yB, double *s_zB, double *weightB, long int NdataB,double NdataBw, double *s_x_ranB, double *s_y_ranB, double *s_z_ranB, double *weight_ranB, long int NrandB, double L1, double L2, int ngrid, double P_shot_noise, double P_shot_noiseB, double Deltakbis, double I33,double I33B,double I22,double I22B, double IN,double INB, double Bsn,double BsnB,  double alpha,double alphaB, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out, char *name_bsAAB_out,char *name_bsABA_out,char *name_bsBAA_out,char *name_bsABB_out,char *name_bBABs_out,char *name_bsBBA_out,char *name_bsBBB_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *do_multigrid, char *bis_opt, char *triangle_shapes,char *do_bispectrum2, char *type_of_code, char *name_density,char *type_of_input,char *type_of_survey);


void loop_bispectrum_periodic_for_gadget_caller(double kmin,double kmax, int Ninterlacing, double L1, double L2, int ngrid, double Deltakbis, int mode_correction, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bsAAB_out,char *name_bsABA_out,char *name_bsBAA_out,char *name_bsABB_out,char *name_bsBAB_out,char *name_bsBBA_out,char *name_bsBBB_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *name_data_in,int gadget_files,char *name_dataB_in,int gadget_filesB, char *do_multigrid, char *bis_opt, char *triangle_shapes,char *RSD,char *RSDB, char *do_bispectrum2,char *type_of_code,char *type_of_survey);

void loop_bispectrum_periodic_for_gadget_x_ascii_caller(double kf,double kny,int Ninterlacing, double L1, double  L2, int ngrid, double  Deltakbis, int  mode_correction, int n_lines_parallel, char *binning_type, char *name_bsAAA_out, char *name_bsAAB_out, char *name_bsABA_out, char *name_bsBAA_out,char *name_bsABB_out, char *name_bsBAB_out,char *name_bsBBA_out, char *name_bsBBB_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *type_of_mass_assigment, char *triangles_num, char *write_triangles, char *triangles_id, char *name_data_in, int gadget_files,  double *pos_xB, double *pos_yB, double *pos_zB, double *weightB, long int Ndata2B,double Ndata2Bw, double P_shot_noiseB,char *do_multigrid, char *bis_opt, char *triangle_shapes,char *RSDB,char *do_bispectrum2,char *type_of_code,int reverse,char *type_of_survey);


void bispectrum_calculator(double *K1,double **K2, double ***K3, int Nk1, int *Nk2, int **Nk3, double *deltak_re, double *deltak_im,double *deltak_reB, double *deltak_imB,double *deltak2_re, double *deltak2_im,int Ninterlacing, double kmin,double kmax,double L1,double L2, int ngrid, double Deltakbis, double I33,double I22, double IN, double Bsn, double Psn,double I33B,double I22B, double INB, double BsnB, double PsnB, int n_lines_parallel, char *binning_type, char *name_bs_out,char *name_bs002_out,char *name_bs020_out,char *name_bs200_out, char *triangles_num, char *write_triangles, char *triangles_id, long int *input_bin, char *do_bispectrum2, char *type_cross,char *type_of_survey, char *bis_opt, double *Kuni, int Nk, int *vector_pointerk1, int **vector_pointerk2, int ***vector_pointerk3);

