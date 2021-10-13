#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define pi 3.141592653589793238462643
#define degrad 0.0174532925199

//using namespace std;

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

namespace utils {

  float ran2(long *idum);
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


  void polar2cart_3d(double *d, double *e) {
    // ra,dec,r --> x,y,z
    double ra=d[0]*degrad;
    double dec=d[1]*degrad;
    double r=d[2];
    e[0]=r*cos(ra)*cos(dec);
    e[1]=r*sin(ra)*cos(dec);
    e[2]=r*sin(dec);
  }

  void cart2polar_3d(double *d, double *e) {
    // x,y,z -> ra,dec,r
    e[0]=atan(d[1]/d[0])/degrad; if (d[0]<0) e[0]+=180; if (e[0]<0) e[0]+=360;
    e[1]=atan(d[2]/sqrt(d[0]*d[0]+d[1]*d[1]))/degrad;
    e[2]=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
  }

  double sqsum_div(double *f, double *g, const int& np) {
    double y=0.;
    for (int i=0;i<np;i++) y+=(f[i]/g[i])*(f[i]/g[i]);
    return y;
  }

  double sqsum_wei(double *f, double *w, const int& np) {
    double y=0.;
    for (int i=0;i<np;i++) y+=w[i]*f[i]*f[i];
    return y;
  }

  double invsum(double *f, const int& np) {
    double y=0.;
    for (int i=0;i<np;i++) y+=1/f[i];
    return y;
  }

  double effvol(double *nden, const int& np) {
    const double pk=2e4;
    double sum=0.;
    for (int i=0;i<np;i++) sum+=pk/(1+nden[i]*pk);
    return sum;
  }

  double weisum(double **d, const int& np) {
    double y=0.;
    //for (int i=0;i<np;i++) y+=d[i][7]/d[i][5];
    for (int i=0;i<np;i++) y+=d[i][7];
    return y;
  }

  void zweight(const double& zmin, const double& zmax, double *zwei, const int& nzwei, double *z, double *nden, double *w, double *bg, const int& np) {
    double zmean[nzwei];
    for (int i=0;i<nzwei;i++) {zwei[i]=0; zmean[i]=0;}
    for (int i=0;i<np;i++) {
      int iz=int((z[i]-zmin)/(zmax-zmin)*nzwei);
      zwei[iz]+=nden[i]*w[i]*w[i]/bg[i]/bg[i];
      zmean[iz]+=z[i]*nden[i]*w[i]*w[i]/bg[i]/bg[i];
    }
    double tot=0;
    for (int i=0;i<nzwei;i++) zmean[i]/=zwei[i];
    for (int i=0;i<nzwei;i++) tot+=zwei[i];
    for (int i=0;i<nzwei;i++) zwei[i]/=tot;
    for (int i=0;i<nzwei;i++) cout << zmean[i] << " " << zwei[i] << endl;
  }


  void ofwriteasc_3d(char*string,double **d, const int& np) {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int i=0;i<np;i++) ofs << d[i][0] << " " << d[i][1] << " " << d[i][2] << endl;
    ofs.close();
  }

  int lineqbin_num(const double& xmin, const double& xmax, const int& nbin, const double& x) {
    return int(floor((x-xmin)/(xmax-xmin)*nbin));
  }

  double boxsize(double *len) {
    return len[0]*len[1]*len[2];
  }

  void indexx(long n, double*arr, long*indx)  {
    const int NSTACK=50,M=7;
    long i,indxt,ir=n-1,itemp,j,k,l=0;
    int jstack=0,*istack=new int[NSTACK+1];
    double a;
    
    for (j=0;j<n;j++) indx[j]=j;
    for (;;) {
      if (ir-l < M) {
	for (j=l+1;j<=ir;j++) {
	  indxt=indx[j];
	  a=arr[indxt];
	  for (i=j-1;i>=0;i--) {
	    if (arr[indx[i]] <= a) break;
	    indx[i+1]=indx[i];
	  }
	  indx[i+1]=indxt;
	}
	if (jstack == 0) break;
	ir=istack[jstack--];
	l=istack[jstack--];
      } else {
	k=(l+ir) >> 1;
	SWAP_long(indx[k],indx[l+1]);
	if (arr[indx[l+1]] > arr[indx[ir]]) SWAP_long(indx[l+1],indx[ir]);
	if (arr[indx[l]] > arr[indx[ir]]) SWAP_long(indx[l],indx[ir]);
	if (arr[indx[l+1]] > arr[indx[l]]) SWAP_long(indx[l+1],indx[l]);
	i=l+1;
	j=ir;
	indxt=indx[l];
	a=arr[indxt];
	for (;;) {
	  do i++; while (arr[indx[i]] < a);
	  do j--; while (arr[indx[j]] > a);
	  if (j < i) break;
	  SWAP_long(indx[i],indx[j]);
	}
	indx[l]=indx[j];
	indx[j]=indxt;
	jstack += 2;
	if (jstack > NSTACK) throw ios::failure("NSTACK too small in indexx.");
	if (ir-i+1 >= j-l) {
	  istack[jstack]=ir;
	  istack[jstack-1]=i;
	  ir=j-1;
	} else {
	  istack[jstack]=j-1;
	  istack[jstack-1]=l;
	  l=i;
	}
      }
    }
    delete [] istack;
  }

  void spline(double*x, double*y, const int n, double*y2) {
    double *u=new double[n];
    double yp1=1.e30,ypn=1.e30;

    //*y2 = -0.5;
    //*u=(3.0/(*(x+1)-*x))*((*(y+1)-*y)/(*(x+1)-*x)-yp1);
    *y2=0., *u=0.;
    for (int i=1;i<n-1;i++) {
      double sig=(*(x+i)-*(x+i-1))/(*(x+i+1)-*(x+i-1));
      double p=sig*(*(y2+i-1))+2.0;
      *(y2+i)=(sig-1.0)/p;
      *(u+i)=(*(y+i+1)-*(y+i))/(*(x+i+1)-*(x+i)) 
	- (*(y+i)-*(y+i-1))/(*(x+i)-*(x+i-1));
      *(u+i)=(6.0*(*(u+i))/(*(x+i+1)-*(x+i-1))-sig*(*(u+i-1)))/p;
    }
    double qn=0.,un=0.; 
    //qn=0.5;
    //un=(3.0/(*(x+n-1)-*(x+n-2)))*(ypn-(*(y+n-1)-*(y+n-2))/(*(x+n-1)-*(x+n-2)));
    *(y2+n-1)=(un-qn*(*(u+n-2)))/(qn*(*(y2+n-2))+1.0);
    for (int k=n-2;k>=0;k--) *(y2+k)=*(y2+k)*(*(y2+k+1))+*(u+k);
    delete [] u;
  }

  void splint(double*xa, double*ya, double*y2a, const int n, const double& x, double *y) {
    int klo=0,khi=n-1,k;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (*(xa+k) > x) khi=k;
      else klo=k;
    }
    double h=*(xa+khi)-*(xa+klo);
    if (h == 0.0) throw ios::failure("Bad xa input to routine splint");
    double a=(*(xa+khi)-x)/h, b=(x-*(xa+klo))/h;
    *y=a*(*(ya+klo))+b*(*(ya+khi))+((a*a*a-a)*(*(y2a+klo))+(b*b*b-b)*(*(y2a+khi)))*(h*h)/6.0;
  }

  double**initmat2(const int& n1,const int& n2) {
    double**f=new double*[n1];
    f[0]=new double[n1*n2];
    for (int i=1;i<n1;i++) f[i]=f[0]+i*n2;
    for (int i=0;i<n1;i++) for (int j=0;j<n2;j++) f[i][j]=0;
    return f;
  }

  double***initmat3(const int& n1,const int& n2,const int& n3) {
    // pointers to first level
    double***f=new double**[n1];

    // pointers to second level
    f[0]=new double*[n1*n2];
    for (int i=1;i<n1;i++) f[i]=f[i-1]+n2;

    // pointers to third level
    f[0][0]=new double[n1*n2*n3];
    for (int i=1;i<n1;i++) f[i][0]=f[i-1][0]+n2*n3;
    for (int i=0;i<n1;i++) for(int j=1;j<n2;j++) f[i][j]=f[i][j-1]+n3;
    for (int i=0;i<n1;i++) for(int j=0;j<n2;j++) for(int k=0;k<n3;k++) f[i][j][k]=(double)0;
    return f;
  }

  double gauss(const double& k) {
    return exp(-pow(k,2)/2.);
  }

  double norm(double*x,const int& n) {
    double sm=0.;
    for (int i=0;i<n;i++) sm+=x[i]*x[i];
    return sqrt(sm);
  }

  double dist(double*x,double*y,const int& n) {
    double sm=0.;
    for (int i=0;i<n;i++) sm+=(x[i]-y[i])*(x[i]-y[i]);
    return sqrt(sm);
  }

  void SWAP_long(long& a, long& b) {
    long c;
    c=a;
    a=b;
    b=c;
  }

  void ofwritefunc(char*string, double*x, double*y, const int& n) {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int i=0;i<n;i++) {
      ofs << x[i] << " " << y[i] << " " << endl;
    }
    ofs.close();
  }

  void ofwritefunc_ncol(char*string, double*x, double**y, const int& n, const int& m) {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int i=0;i<n;i++) {
      ofs << x[i] << " ";
      for (int j=0;j<m;j++) ofs << y[i][j] << " ";
      ofs << endl;
    }
    ofs.close();
  }

}
