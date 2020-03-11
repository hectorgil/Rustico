
//void reshuffle(double radata[],double decdata[],double wcoldata[],double wsysdata[],long int npar_used,long int seed);//this re-shuffles the indices on angular-data wrt the unchanged radial-data in case of both-shuffling option.
void reshuffle(double radata[],double decdata[],double wcoldata[],double wsysdata[], long int npar_used,long int seed);

double P_interpol(double k0, double *k, double *P, int N);
int countlines(char *filename);

void freeTokens(double** tokens, int N);

void freeTokensInt(int** tokens, int N);

void freeTokensLInt(long int** tokens, int N);

void freeTokensLInt3(long int*** tokens, int N1, int N2);

void freeTokens2(double ***tokens, int N1, int *N2);

void freeTokensInt2(int ***tokens,int N1,int *N2);

void freeTokens3(double ***tokens, int N1, int N2);

void freeTokensInt3(int ***tokens,int N1,int N2);

void check_box_for_yamamoto_old(double *parametersL1, int ngrid);

void check_box_for_yamamoto(double *parametersL1, int ngrid);

//void free3D(double ***arr, int p, int c);

//double *** Create3D(int p, int c, int r);
