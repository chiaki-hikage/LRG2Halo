#include<iostream>
#include<time.h>
#include"main.h"
#include"cpara.h"
#include"utils.h"
#include"field.h"

#define pi 3.141592653589793238462643
#define degrad pi/180

using std::cout;
using std::endl;

int main(int argc, char** argv){

  Cpara cp;
  cp.defaultval();
  long idum2=-1;

  cout << cp.ombh2 << endl;
  cout << ran2(&idum2) << endl;
  return 0;

  Parlist par;
  par.defaultval();
  par.argument(argc,argv);

  double len[3]={par.len, par.len, par.len};
  int npix[3]={par.npix,int(par.npix*len[1]/len[0]),int(par.npix*len[2]/len[0])};
  double ***f=utils::initmat3(npix[0],npix[1],npix[2]+2);
  double ***g=utils::initmat3(npix[0],npix[1],npix[2]+2);
  double alpha;
  double fran; // Ngal/Nran
  double beta=1; // Nhalo/Ngal
  double pshot; // shot noise power
  double vol;  // survey volume
  double norm; // normalization of weighted galaxy fluctuation

  const int nz=30,nzwei=3;
  double bz[nz],bz2[nz],zdist_lrg[nz],zdist_halo[nz],zwei[nzwei];

  for (int k=0;k<2;k++) {
    double wisum_d,wisum_r;
    char inf[100], datatype[1];
    if (k==0) strcpy(datatype,"d"); else strcpy(datatype,"r");  
    if (datatype[0]=='d') strcpy(inf,par.infd); 
    if (datatype[0]=='r') strcpy(inf,par.infr); 
    cout << inf << endl;

    int np=field::datanum_LRG(inf);
    int jmax=11;
    double **d=utils::initmat2(np,jmax); // 0: ra, 1: dec, 2: dis, 3: z, 4: Mg, 5: comp, 6: nden, 7: fcol, 8: xi, 9: yi, 10: zi
    field::ifreadasc_LRG(inf,d,np,datatype);
    field::region(d,np,jmax,par.dir);
    double z[np]; for (int i=0;i<np;i++) z[i]=d[i][3];
    double r[np]; field::ztoangdis(cp,z,r,np);
    for (int i=0;i<np;i++) d[i][2]=r[i];
    for (int i=0;i<np;i++) utils::polar2cart_3d(d[i],d[i]+jmax-3);

    double nden[np]; for (int i=0;i<np;i++) nden[i]=d[i][6];
    cout << "survey volume:" << utils::invsum(nden,np) << endl;
    cout << "effective volume:" << utils::effvol(nden,np) << endl;

    if (datatype[0]=='d') wisum_d=utils::weisum(d,np);
    if (datatype[0]=='r') {
      wisum_r=utils::weisum(d,np);
      fran=wisum_d/wisum_r;
      cout << "fran (=Ngal/Nran):" << fran << endl;
      for (int i=0;i<np;i++) d[i][6]*=fran;
    }

    int mul[np]; 
    if (par.type[0]=='h'&&datatype[0]=='d') {
      double w[np]; 
      for (int i=0;i<np;i++) w[i]=d[i][7]/d[i][5];
      field::zhist(z,w,np,par.zmin,par.zmax,nz,cp,zdist_lrg);
      int g[np]; for (int i=0;i<np;i++) g[i]=-1;
      long indx[np]; utils::indexx((long)np,z,indx);
      field::fofgroup_Reid(d,indx,np,g);
      beta/=np;
      field::groupcatalog(d,g,np,jmax,par.type,cp,mul);
      beta*=np;
      cout << "Nhalo/Ngal:" << beta << endl;
      for (int i=0;i<np;i++) {z[i]=d[i][3]; w[i]=d[i][7]/d[i][5];}
      field::zhist(z,w,np,par.zmin,par.zmax,nz,cp,zdist_halo);
    }

    // change nden for halo catalogs
    if (par.type[0]=='h') {
      double x[nz],y[nz],y2[nz];
      for (int i=0;i<nz;i++) {
	x[i]=par.zmin+(par.zmax-par.zmin)/nz*(i+0.5);
	y[i]=zdist_halo[i]/zdist_lrg[i];
      }
      utils::spline(x,y,nz,y2);
      int np2=np;
      long idum=-1;
      np=0;
      for (int i=0;i<np2;i++) {
	double y3;
	if (z[i]<x[0]) y3=y[0]+(y[1]-y[0])/(x[1]-x[0])*(z[i]-x[0]);
	else if (z[i]>x[nz-1]) y3=y[nz-1]+(y[nz-2]-y[nz-1])/(x[nz-2]-x[nz-1])*(z[i]-x[nz-1]);
	utils::splint(x,y,y2,nz,z[i],&y3);
	if (datatype[0]=='r') {
	  if (ran2(&idum)<y3) {
	    for (int j=0;j<jmax;j++) d[np][j]=d[i][j];
	    d[np][6]=d[np][6]*y3;
	    np++;
	  }
	} else {
	  d[np][6]=d[np][6]*y3;
	  np++;
	}
      }
    }
    cout << "np:" << np << endl;
    for (int i=0;i<np;i++) nden[i]=d[i][6];
    
    double bg[np],mbg[np],mbg2[np];
    for (int i=0;i<np;i++) z[i]=d[i][3];
    if (datatype[0]=='d') {
      for (int i=0;i<np;i++) bg[i]=field::lumbias(d[i][4]);
      field::biashist(z,bg,np,par.zmin,par.zmax,nz,bz,bz2);
    }

    for (int i=0;i<np;i++) {
      mbg[i]=field::allocbias_z(par.zmin,par.zmax,nz,bz,z[i]);
      mbg2[i]=field::allocbias_z(par.zmin,par.zmax,nz,bz2,z[i]);
    }
    if (datatype[0]=='r') for (int i=0;i<np;i++) bg[i]=mbg[i];
    double w[np]; for (int i=0;i<np;i++) w[i]=field::optwei_LRG(bg[i],mbg2[i],d[i][6],d[i][5],d[i][7]);
   

     
    // compute alpha
    double sa=0,sa0;
    //for (int i=0;i<np;i++) sa+=d[i][7]/d[i][5]/bg[i];
    for (int i=0;i<np;i++) sa+=1./d[i][5]/bg[i];
    cout << "sa:" << sa << endl;

    if (datatype[0]=='d') sa0=sa; 
    else {
      cout << sa0/sa << endl;
      alpha=sa0/sa;
    }

    if (datatype[0]=='r') {
      vol=fran*utils::invsum(nden,np);
      cout << "survey volume:" << vol << endl;
      norm=sqrt(alpha*utils::sqsum_wei(w,nden,np));
      cout << "normalization:" << norm << endl;
      pshot=alpha*(1+alpha)*utils::sqsum_div(w,bg,np)/pow(norm,2);
      cout << "P_shot:" << pshot << endl;
      utils::zweight(par.zmin,par.zmax,zwei,nzwei,z,nden,w,bg,np);
    }
    // ra,dec,dis -> x,y,z
    for (int i=0;i<np;i++) for (int j=0;j<3;j++) d[i][j]=d[i][j+jmax-3];
    utils::ofwriteasc_3d("map.dat",d,np);
    //utils::unitboxin(d,np,len);
    
    if (datatype[0]=='d') field::densassign_3d_wb(d,w,bg,len,np,f,npix);  
    else field::densassign_3d_wb(d,w,bg,len,np,g,npix);

    delete [] *d; delete [] d;
  }

  double fvol=vol/utils::boxsize(len);
  cout << "volume fraction:" << fvol << endl;

  field::denfrac_3d_runsub(f,g,alpha,norm,npix);
  //for (int k=0;k<npix[2];k+=8) utils::slicemap(f,npix,k);

  field::pk3d_1d(g,len,npix,par.kmin,par.kmax,par.nkbin,pshot*alpha/(1+alpha),par.outfwin);
  delete [] **g; delete [] *g; delete[] g;

  field::pk3d_1d(f,len,npix,par.kmin,par.kmax,par.nkbin,pshot,par.outf);
  delete [] **f; delete [] *f; delete[] f;

  return 0;
}

#undef pi
#undef degrad
