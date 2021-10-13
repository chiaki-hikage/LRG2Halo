namespace field {
  void ztoangdis(Cpara cp, double *z, double *r, const int &np);
  void zhist(double *z, double *w, const int &np, const double& zmin, const double& zmax, const int& nz, Cpara cp, double *zdist);
  int datanum_LRG(const char *inf);
  void ifreadasc_LRG(const char *inf, double **d, const int& np, const char *type);
  void datarange(double **d, const int& n, double* len, double *dmin);
  void fofgroup_Reid(double **d, long *indx, const int &np, int *g);
  void groupcatalog(double **d, int *g, int& np, const int& jmax, char *type, Cpara cp, int *mul);
  double lumbias(const double& m);
  void biashist(double *z, double *bg, const int &np, const double& zmin, const double& zmax, const int& nz, double *bz, double *bz2);
  double allocbias_z(const double& zmin, const double& zmax, const int& nz, double *bz, const double& z);
  double optwei_LRG(const double& bg, const double& mbg2, const double& nden, const double& comp, const double& fcol);
  void region(double **d, int &np, const int &jmax, char*string);
  void pk3d_1d(double ***f, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, const double& shotn, char* outf);
  void densassign_3d_wb(double **d, double *w, double *b, double *len, const int& n, double ***f, int *npix);
  void denfrac_3d_runsub(double ***f, double ***g, const double& alpha, const double& norm, int *npix);
}
