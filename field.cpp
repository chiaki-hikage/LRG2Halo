#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define pi 3.141592653589793238462643
#define degrad 0.0174532925199

#include"evolve.h"
#include"cpara.h"
#include"utils.h"
#include"fourier_ver3.h"

namespace field {

  void ztoangdis(Cpara cp, double *z, double *r, const int &np) {
    int nzmax=500;
    double zlim=10., zp[nzmax],rp[nzmax],rp2[nzmax];
    for (int i=0;i<nzmax;i++) {
      zp[i]=zlim/(nzmax-1)*i;
      rp[i]=evolve::angdis(cp,zp[i]);
    }
    utils::spline(zp,rp,nzmax,rp2);
    for (int n=0;n<np;n++) {
      if (z[n]>zp[nzmax-1]||z[n]<zp[0]) {
	cout << n << " " << z[n] << endl;
	throw ios::failure("set zlim to larger in function ztoangdis");
      }
      utils::splint(zp,rp,rp2,nzmax,z[n],r+n);
    }
  }

  void zhist(double *z, double *w, const int &np, const double& zmin, const double& zmax, const int& nz, Cpara cp, double *zdist) {
    for (int i=0;i<nz;i++) zdist[i]=0;
    for (int i=0;i<np;i++) {
      int iz=utils::lineqbin_num(zmin,zmax,nz,z[i]);
      if (iz>=0&&iz<nz) zdist[iz]+=w[i];
    }
    double zmean[nz];
    for (int i=0;i<nz;i++) {
      double z1=zmin+(zmax-zmin)/nz*i, rmin=evolve::angdis(cp,z1);
      double z2=zmin+(zmax-zmin)/nz*(i+1), rmax=evolve::angdis(cp,z2);
      zmean[i]=zmin+(zmax-zmin)/nz*(i+0.5);
      //!!! check effective sky fraction
      double fsky=pi;
      //double fsky=0.721434*pi;
      zdist[i]/=fsky/3*(pow(rmax,3)-pow(rmin,3));
    }
  }

int datanum_LRG(const char *inf) {
  ifstream ifs;
  ifs.open(inf,ios::in); 
  if (!ifs) throw ios::failure("Failed to open input file");
  int np;
  ifs >> np;
  cout << "Number of objs:" << np << endl;
  ifs.close();
  return np;
}

void ifreadasc_LRG(const char *inf, double **d, const int& np, const char* type) {
  ifstream ifs;
  ifs.open(inf,ios::in); 
  if (!ifs) throw ios::failure("Failed to open input file");
  int ilss,icomb,sector;
  double rwei,fogt,ignore;
  int np2; ifs >> np2;
  for (int i=0;i<np;i++) {
    if (type[0]=='r') {
      ifs >> d[i][0] >> d[i][1] >> d[i][3] >> d[i][5] >> d[i][6] >> rwei >> icomb >> sector;
      d[i][6]*=1e-4;
      d[i][4]=0;  // Mg
      d[i][7]=1;  // fcol
    } else {
      ifs >> d[i][0] >> d[i][1] >> d[i][3] >> d[i][4] >> d[i][5] >> d[i][6] >> rwei >> d[i][7] >> fogt >> ilss >> icomb >> sector;
      d[i][6]*=1e-4;
    }
  }
  cout << "wi sum:" << utils::weisum(d,np) << endl;
  ifs.close();
}

  void datarange(double **d, const int& n, double*len, double*dmin) {
    double dmax[3];
    for (int j=0;j<3;j++) {
      for (int i=0;i<n;i++) {
	if (i==0) {
	  dmin[j]=d[i][j]; 
	  dmax[j]=d[i][j];
	} else {
	  if (d[i][j]<dmin[j]) dmin[j]=d[i][j];
	  if (d[i][j]>dmax[j]) dmax[j]=d[i][j];
	}
      }
    }
    for (int j=0;j<3;j++) cout << j << " " << dmax[j]-dmin[j] << " " << len[j] << endl;
    for (int j=0;j<3;j++) for (int i=0;i<n;i++) {
      if (d[i][j]>len[j]+dmin[j]||d[i][j]<dmin[j]) {
	cout << j << " " << d[i][j] << " " << dmin[j] << endl;
	pause();
      }
    }
  }

  void fofgroup_Reid(double **d, long *indx, const int &np, int *g) {
    int ig=0,jmax=3;
    cout << np << endl;
    for (int i=0;i<np;i++) {
      if (g[indx[i]]==-1) {g[indx[i]]=ig; ig++;
	//cout << indx[i] << " " << ig << endl;
	//pause();
      }
      if (i%10000==0) cout << ig << endl;
      const double drlim=0.8;
      double z0=d[indx[i]][3],r0=d[indx[i]][2];
      double ra0=d[indx[i]][0], dec0=d[indx[i]][1];
      int j=i+1; 
      double z1;if (j<np) z1=d[indx[j]][3]; else z1=z0;
      double dzlim=0.006*(1+z1);
      while (z1-z0<dzlim&&j<np) {
	double ra1=d[indx[j]][0], dec1=d[indx[j]][1];
	double dra=fabs(ra1-ra0)*degrad, ddec=fabs(dec1-dec0)*degrad;
	if (r0*cos(dec0*degrad)*dra<2*drlim&&r0*ddec<2*drlim) {
	  double d0[3]={ra0,dec0,r0},e0[jmax]; utils::polar2cart_3d(d0,e0);
	  double d1[3]={ra1,dec1,r0},e1[jmax]; utils::polar2cart_3d(d1,e1);
	  if (utils::dist(e0,e1,jmax)<drlim) {
	    //g[indx[j]]=g[indx[i]]; 
	    if (g[indx[j]]==-1) g[indx[j]]=g[indx[i]]; 
	    else {
	      for (int i1=0;i1<i;i1++) if (g[indx[i1]]==g[indx[i]]) g[indx[i1]]=g[indx[j]];
	      for (int i1=i+1;i1<np;i1++) if (g[indx[i1]]==g[indx[i]]) g[indx[i1]]=g[indx[j]];
	      g[indx[i]]=g[indx[j]];
	    }
	  }
	}
	j++;
	if (j<np) z1=d[indx[j]][3]; else z1=z0;
	dzlim=0.006*(1+z1);
      }
    }

    int mul[ig],ig0=0;
    for (int i=0;i<ig;i++) mul[i]=0;
    for (int i=0;i<np;i++) {
      if (g[i]<0) cout << i << " " << g[i] << endl;
      mul[g[i]]++;
    }
    for (int i=ig;i>=0;i--) {
      if (mul[i]==0) {
	ig0++;
	for (int j=0;j<np;j++) if (g[j]>i) g[j]--;
      }
    }
    ig-=ig0;

    cout << "nhalo " << ig << endl;
    cout << "ndata " << np << endl;

    for (int i=0;i<ig;i++) mul[i]=0;
    for (int i=0;i<np;i++) mul[g[i]]++;

    ofstream ofs;
    ofs.open("output/groupnum.dat",ios::out);
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int i=0;i<np;i++) {
      ofs << i << " " << g[i] << endl;
    }
    ofs.close();
  }

  void groupcatalog(double **d, int *g, int& np, const int& jmax, char *type, Cpara cp, int *mul) {
    double **h=utils::initmat2(np,jmax);
    for (int i=0;i<np;i++) mul[i]=0;
    int gnum=0;
    for (int i=0;i<np;i++) {
      if (mul[g[i]]==0) {
	for (int j=0;j<jmax;j++) h[g[i]][j]=d[i][j];
	//!!! fcol = 1
	h[g[i]][7]=1;  // fcol
	gnum++; 
      } else {
	//!!! Keep Brightest data in each group
	if (type[0]=='b') {
	  if (d[i][4]<h[g[i]][4]) for (int j=0;j<jmax;j++) h[g[i]][j]=d[i][j];
	} else if (type[0]=='f') {
	  if (d[i][4]<h[g[i]][4]) {
	    h[g[i]][4]=d[i][4];  // Mg
	  } else {
	    for (int j=0;j<4;j++) h[g[i]][j]=d[i][j];
	    for (int j=5;j<jmax;j++) h[g[i]][j]=d[i][j];
	  }
	} else if (type[0]=='m') {
	  if (d[i][4]<h[g[i]][4]) h[g[i]][4]=d[i][4];  // Mg
	  for (int j=3;j<4;j++) h[g[i]][j]=(h[g[i]][j]*mul[g[i]]+d[i][j])/(mul[g[i]]+1);
	  for (int j=5;j<jmax;j++) h[g[i]][j]=(h[g[i]][j]*mul[g[i]]+d[i][j])/(mul[g[i]]+1);
	} else {
	  cout << "no type" << endl;
	  throw ios::failure("Failed to open output file");
	}
	//!!! fcol = 1
	h[g[i]][7]=1;  // fcol
      }
      mul[g[i]]++;
    }

    if (type[0]=='m') for (int i=0;i<np;i++) utils::cart2polar_3d(h[g[i]]+jmax-3,h[g[i]]);

    /*
    int flag[np]; for (int i=0;i<np;i++) flag[i]=0;
    for (int i=0;i<np;i++) {
      if (mul[g[i]]>1&&flag[g[i]]==0) {
	ofs << h[g[i]][0] << " " << h[g[i]][1] << " " << h[g[i]][3] << " "; //ra,dec,z
	ofs << h[g[i]][4] << " ";  // Mg
	ofs << mul[g[i]] << endl; // number of LRGs in the halo
	flag[g[i]]=1;
      }
    }
    ofs.close();
    */

    int hist[10];
    for (int j=0;j<10;j++) hist[j]=0;
    for (int j=0;j<gnum;j++) hist[mul[j]]++;
    for (int j=0;j<10;j++) {
      cout << j << " " << hist[j] << endl;
    }

    np=gnum;
    for (int i=0;i<np;i++) for (int j=0;j<jmax;j++) d[i][j]=h[i][j];
    cout << "number of groups: " << np << endl;
    delete [] *h; delete [] h;
  }

  double lumbias(const double& m) {
    //Tegmark et al. 2004;
    const double mstar=-20.83;
    //!!! check
    const double bstar=1.7;
    return bstar*(0.895+0.15*pow(10,-0.4*(m-mstar))-0.04*(m-mstar));
  }

  void biashist(double *z, double *bg, const int &np, const double& zmin, const double& zmax, const int& nz, double *bz, double *bz2) {
    int nbz[nz];
    double zmean[nz],mbias=0.;
    for (int i=0;i<nz;i++) {nbz[i]=0; bz[i]=0.; bz2[i]=0.;}
    for (int i=0;i<np;i++) {
      int iz=utils::lineqbin_num(zmin,zmax,nz,z[i]);
      if (iz>=0&&iz<nz) {
	bz[iz]+=bg[i];
	bz2[iz]+=bg[i]*bg[i];
	nbz[iz]++;
      }
      mbias+=bg[i];
    }
    for (int i=0;i<nz;i++) {
      zmean[i]=zmin+(zmax-zmin)/nz*(i+0.5);
      bz[i]/=nbz[i];
      bz2[i]/=nbz[i];
    }
    mbias/=np;
    cout << "Mean bias:" << mbias << endl;
    utils::ofwritefunc("output/bz.dat",zmean,bz,nz);
    //for (int i=0;i<nz;i++) cout << zmin+(zmax-zmin)/nz*(i+0.5) << " " << bz[i] << endl;
  }

  double allocbias_z(const double& zmin, const double& zmax, const int& nz, double *bz, const double& z) {
    return bz[min(max(utils::lineqbin_num(zmin,zmax,nz,z),0),nz-1)];
  }

  double optwei_LRG(const double& bg, const double& mbg2, const double& nden, const double& comp, const double& fcol) {
    const double plrg=4e4;
    return bg*bg*fcol/comp/(1+nden*mbg2*plrg);
  }

  void region(double **d, int &np, const int &jmax, char*string) {
    int n=0;
    double dmin=0,dmax=360;
    if (string[0]=='e') {dmin=90; dmax=180;}
    if (string[0]=='w') {dmin=180; dmax=270;}
    if (string[0]=='n') {dmin=90; dmax=270;}
    for (int i=0;i<np;i++) {
      if (d[i][0]>dmin&&d[i][0]<dmax) {
	for (int j=0;j<jmax;j++) d[n][j]=d[i][j];
	n++; 
      }
    }
    cout << "number of data:" << n << endl;
    np=n;
  }

  void densassign_3d_wb(double **d, double *w, double *b, double *len, const int& n, double ***f, int *npix) {
    double dmin[3]; datarange(d,n,len,dmin);
    int ip[3],ip1[3];
    double fp[3];
    for (int i=0;i<n;i++) {
      for (int j=0;j<3;j++) {
	ip[j]=int(floor((d[i][j]-dmin[j])/len[j]*npix[j]-0.5));
	fp[j]=(d[i][j]-dmin[j])/len[j]*npix[j]-0.5-ip[j];
	ip[j]=(ip[j]+npix[j])%npix[j];
	ip1[j]=(ip[j]+1)%npix[j];
      }
      double fac=w[i]/b[i];
      f[ip[0]][ip[1]][ip[2]]+=(1-fp[0])*(1-fp[1])*(1-fp[2])*fac;
      f[ip1[0]][ip[1]][ip[2]]+=fp[0]*(1-fp[1])*(1-fp[2])*fac;
      f[ip[0]][ip1[1]][ip[2]]+=(1-fp[0])*fp[1]*(1-fp[2])*fac;
      f[ip[0]][ip[1]][ip1[2]]+=(1-fp[0])*(1-fp[1])*fp[2]*fac;
      f[ip1[0]][ip1[1]][ip[2]]+=fp[0]*fp[1]*(1-fp[2])*fac;
      f[ip[0]][ip1[1]][ip1[2]]+=(1-fp[0])*fp[1]*fp[2]*fac;
      f[ip1[0]][ip[1]][ip1[2]]+=fp[0]*(1-fp[1])*fp[2]*fac;
      f[ip1[0]][ip1[1]][ip1[2]]+=fp[0]*fp[1]*fp[2]*fac;
    }
  }

  void denfrac_3d_runsub(double ***f, double ***g, const double& alpha, const double& norm, int *npix) {
    for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) {
      f[i][j][k]/=norm;
      g[i][j][k]*=alpha/norm;
      f[i][j][k]-=g[i][j][k];
    }
  }

  void pk3d_1d(double ***f, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, const double& shotn, char* outf) {
    // logarithmic equal binning between kmin and kmax with number of nkbin

    fourier::westward_ho_3d(f,npix[0],npix[1],npix[2]);
    double nbin[nkbin], dnbin[nkbin], **dpk=utils::initmat2(nkbin,3), **pk=utils::initmat2(nkbin,3);
    for (int i=0;i<nkbin;i++) {nbin[i]=0;}
    double knyq=2*pi/len[0]*npix[0]/2;
      
    for (int k=0;k<=npix[2]/2;k++) {
      for (int i=0;i<nkbin;i++) {dnbin[i]=0; for (int j=0;j<3;j++) dpk[i][j]=0.;}
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) {
	int i1=i; if (i>npix[0]/2) i1=i-npix[0];
	int j1=j; if (j>npix[1]/2) j1=j-npix[1];
	double kv[3]={(float)i1/(npix[0]/2),(float)j1/(npix[1]/2),(float)k/(npix[2]/2)}, kabs=utils::norm(kv,3);
	if (kabs>=kmin&&kabs<kmax) {
	  double mu=((float)k/(npix[2]/2))/kabs;
	  int ibin=int((kabs-kmin)/(kmax-kmin)*nkbin);
	  double wkx,wky,wkz,wk2;
	  if (i1*j1*k==0) wk2=1.; else {
	    wkx=sin(pi*i1/npix[0])/(pi*i1/npix[0]);
	    wky=sin(pi*j1/npix[1])/(pi*j1/npix[1]);
	    wkz=sin(pi*k/npix[2])/(pi*k/npix[2]);
	    wk2=pow(wkx*wky*wkz,3.5);
	    //wk2=1.;
	  }
	  double pl[3]={1,5./2.*(3*mu*mu-1),9./8.*(35*pow(mu,4)-30*mu*mu+3)};
	  double s=0.;
	  for (int l=0;l<3;l++) dpk[ibin][l]+=(pow(f[i][j][2*k],2)+pow(f[i][j][2*k+1],2))/wk2*pl[l]*pow(utils::gauss(kabs*2*pi*s/len[0]),2);
	  dnbin[ibin]++;
	} 
      }
      int fac=1; if (k>0) fac=2;
      for (int i=0;i<nkbin;i++) {nbin[i]+=dnbin[i]*fac; for (int j=0;j<3;j++) pk[i][j]+=dpk[i][j]*fac;}
    }
    delete [] *dpk; delete [] dpk;
    cout << "P(k=0)=" << pow(f[0][0][0],2)+pow(f[0][0][1],2) << endl;
    double x[nkbin];
    for (int i=0;i<nkbin;i++) {
      x[i]=knyq*(kmin+(kmax-kmin)/nkbin*(i+0.5));
      if (nbin[i]>0) {
	for (int j=0;j<3;j++) pk[i][j]/=nbin[i];
	pk[i][0]-=shotn;
      }
      if (pk[i][0]!=0) {
	pk[i][1]=pk[i][1]/pk[i][0];
	pk[i][2]=pk[i][2]/pk[i][0];
      }
      //cout << x[i] << " " << pk[i][0] << " "  << pk[i][1] << " " << pk[i][2] << endl;
    }
    utils::ofwritefunc_ncol(outf,x,pk,nkbin,3);
    delete [] *pk; delete [] pk;

  }

}

