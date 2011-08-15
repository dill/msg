void wood_path(int*, int*, double*, double*, int*, double*, double*, double*, double*, int*, double*, double*, int*, int*, double*, double*, double*,int*, int*);
int make_bnd_path(double[2], double[2], int nbnd, double**, node**, int);
void delete_step(node**, int nbnd, double**);
void alter_step(node**, int nbnd, double**);
int has_converged(node*, node*);
int iter_path(node**, int, double**);
void get_euc_path(double*, double*, int, double**, int, double*, int);
int append_path(node**, node**, double[2], double[2],int, int, double**);


