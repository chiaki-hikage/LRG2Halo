float ran2(long *idum);

namespace utils {
  void polar2cart_3d(double *d, double *e);
  void cart2polar_3d(double *d, double *e);
  double sqsum_div(double *f, double *g, const int& np);
  double sqsum_wei(double *f, double *w, const int& np);
  double invsum(double *f, const int& np);
  double effvol(double *nden, const int& np);
  double weisum(double **d, const int& np);
  void zweight(const double& zmin, const double& zmax, double *zwei, const int& nzwei, double *z, double *nden, double *w, double *bg, const int& np);
  void ofwriteasc_3d(char*string,double **d, const int& np);
  int lineqbin_num(const double& xmin, const double& xmax, const int& nbin, const double& x);
  double boxsize(double *len);
  void indexx(long n, double*arr, long*indx);
  void spline(double*x, double*y, const int n, double*y2);
  void splint(double*xa, double*ya, double*y2a, const int n, const double& x, double *y);
  double**initmat2(const int& n1,const int& n2);
  double***initmat3(const int& n1,const int& n2,const int& n3);
  double gauss(const double& k);
  double norm(double *x, const int& n);
  double dist(double*x,double*y,const int& n);
  void SWAP_long(long& a, long& b);
  void ofwritefunc(char*string, double*x, double*y, const int& n);
  void ofwritefunc_ncol(char*string, double*x, double**y, const int& n, const int& m);
}
