#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>


long int count_particles_gadget(char *name_data_in);

void load_snapshot(char *fname, int files, int *params);

void loop_interlacing_periodic_gadget(double kmin,double kmax, int Ninterlacing, char *name_data_in ,int gadget_files, double L1, double L2, int ngrid, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out, char *type_of_mass_assigment,double Shot_noise_factor,char *grid_correction_string, char *RSD,char *do_odd_multipoles,char *do_anisotropy, int Nmu, char *file_for_mu);

void loop_interlacing_skycut_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand, double L1, double L2, int ngrid, double alpha, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re0, double *deltak_im0, double *deltak_re2, double *deltak_im2,char *do_bispectrum2);

void modify_input_for_yamamoto(int mode_yamamoto,int mode_yamamoto2, double in[], int ngrid, double L1, double L2);

void fftw_yamamoto_skycut(int mode_yamamoto, double in[], double deltak_re[], double deltak_im[], int ngrid, double L1, double L2, int mode_mass_ass);

void loop_interlacing_skycut(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight,long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand, long int Nrand, double L1, double L2, int ngrid,  double P_shot_noise, double bin_ps, double I22, double alpha, int mode_correction, int n_lines_parallel, char *binning_type, char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out, char *type_of_mass_assigment, char *do_bispectrum, int window_function);

void loop_interlacing_skycut2(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight,long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand, long int Nrand, double L1, double L2, int ngrid,  double P_shot_noise, double bin_ps, double I22, double alpha, int mode_correction, int n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out, char *type_of_mass_assigment,char *do_bispectrum, int window_function);

void loop_directsum_yamamoto_skycut(double kmin,double kmax, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata2, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand, long int Nrand2, double L1, double L2, int ngrid, double P_shot_noise, double bin_ps, double I22, double alpha, int n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out);

void loop_interlacing_periodic(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double L1, double L2, int ngrid, double P_shot_noise, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out, char *type_of_mass_assigment, char *do_odd_multipoles,char *do_anisotropy, char *do_bispectrum,int Nmu, char *file_for_mu);

void loop_directsum_exact_skycut(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, int ngrid, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type,  char *Quadrupole_type, char *Octopole_type, char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out);

void loop_directsum_yamamoto_skycut_caller(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type, char *Quadrupole_type, char *Octopole_type,char *Hexadecapole_type, char *do_odd_multipoles,char *do_anisotropy, char *name_ps_out, char *type_of_computation);

void loop_interlacing_periodic_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im);

void loop_interlacing_periodic_gadget_for_bispectrum(int Ninterlacing, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im, char *name_data_in,int gadget_files, double *params, char *RSD);

