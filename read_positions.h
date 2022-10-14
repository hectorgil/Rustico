//void get_skycuts_write_density_randoms(char *filename, double parameter_value[], double alpha, char *name_den_out,double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, long int Ndata, char *shuffle);

long int get_number_used_lines_periodic(char *filename, double parameter_value[],int type);

double get_number_used_lines_weighted_periodic(char *filename, double parameter_value[],int type);

void get_skycuts_write_density_data(char *filename, double parameter_value[],char *name_den_out);

void get_periodic_data(char *filename_data, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],int type);

void get_skycuts_data(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[], char *type_normalization_mode, double radata[], double decdata[], double zdata[], double wcoldata[], double wsysdata[],double wfkpdata[], double nzdata[],char *shuffle, double Rsmooth);

//void get_skycuts_randoms(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],char *type_normalization_mode, char *type_normalization_mode2, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata,long int Ndata, char *shuffle);

//void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

long int get_number_used_lines_data(char *filename, double parameter_value[]);

long int get_number_used_lines_randoms(char *filename, double parameter_value[]);

//void get_skycuts_randoms(char *path,char *id, char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],char *type_normalization_mode, char *type_normalization_mode2, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata,long int Ndata, char *shuffle,char *write_shuffled_randoms);


void get_skycuts_write_density_randoms(char *filename, double parameter_value[], double alpha, char *name_den_out, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata, long int Ndata,double alpha_true, char *shuffle);

void get_skycuts_randoms(char *path, char *id, char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],char *type_normalization_mode, char *type_normalization_mode2, double *radata, double *decdata, double *zdata, double *wcoldata, double *wsysdata,double *wfkpdata, double *nzdata,long int Ndata, double alpha_true, char *shuffle, char *write_shuffled_randoms, char *name_out_randoms,double Rsmooth);

