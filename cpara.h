//cpara.h
#ifndef CPARA_H
#define CPARA_H

#include<math.h>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string.h>

using std::cout;
using std::endl;
using std::string;

class Cpara {

 public:
  
  Cpara(){}
  ~Cpara(){}

  int ifid;
  double ombh2,omch2,omegak,ns,tau,omegad,ampk,omnuh2,nu_l,nu_m,nrun,w0,wa;
  double k0,Tcmb,sigma8,gcont,rhoc,gconst;
  double omegam,h0,omegac,omegab,omeganu,rho0,gamma,ch0;
  double da_norm,dz_norm;
  char *infpkl, *infpknl;
  string para;

  void defaultval(){
    // fiducial values (WMAP7)
    ifid=1;
    //    ombh2=0.0226;
    ombh2=0.0228;
    //omch2=0.1108;
    omch2=0.1125;
    omegak=0.;
    ns=0.963;
    tau=0.089;
    //omegad=0.7383;
    omegad=0.727;
    ampk=2.43e-9; // primordial amplitude at k=0.002/Mpc
    omnuh2=0.; nu_l=3.04; nu_m=0;
    nrun=0.;  
    w0=-1;
    wa=0.;
    infpkl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_lin.dat";
    infpknl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_nl.dat";
    k0=0.002;
    Tcmb=2.726;
    //    sigma8=0.84;
    sigma8=0.817;
    gconst=4.3e-9;
    rhoc=3e4/(8*3.141592653589793238462643*gconst); // unit: (M_solar/h)/(Mpc/h)^3
    da_norm=1;
    dz_norm=1;
    ch0=2997.;
    para="fid";

    relpara();
  }

  void defaultval_lasdamas(){
    // fiducial values (WMAP7)
    ifid=1;
    ombh2=0.0196;
    omch2=0.1029;
    omegak=0.;
    ns=1.0;
    tau=0.089;
    omegad=0.75;
    ampk=2.43e-9; // primordial amplitude at k=0.002/Mpc
    omnuh2=0.; nu_l=3.04; nu_m=0;
    nrun=0.;  
    w0=-1;
    wa=0.;
    infpkl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_lin.dat";
    infpknl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_nl.dat";
    k0=0.002;
    Tcmb=2.726;
    sigma8=0.8;
    gconst=4.3e-9;
    rhoc=3e4/(8*3.141592653589793238462643*gconst); // unit: (M_solar/h)/(Mpc/h)^3
    da_norm=1;
    dz_norm=1;
    ch0=2997.;
    para="fid";

    relpara();
  }

  void defaultval_wrong(){
    // change w0 from -1 to -0.7
    ifid=1;
    ombh2=0.0196;
    omch2=0.1029;
    omegak=0.;
    ns=1.0;
    tau=0.089;
    omegad=0.75;
    ampk=2.43e-9; // primordial amplitude at k=0.002/Mpc
    omnuh2=0.; nu_l=3.04; nu_m=0;
    nrun=0.;  
    w0=-0.7;
    wa=0.;
    infpkl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_lin.dat";
    infpknl="/home/hikage/work/FoG/output/camb/wmap_fid_matterpower_z0.0_nl.dat";
    k0=0.002;
    Tcmb=2.726;
    sigma8=0.8;
    gconst=4.3e-9;
    rhoc=3e4/(8*3.141592653589793238462643*gconst); // unit: (M_solar/h)/(Mpc/h)^3
    da_norm=1;
    dz_norm=1;
    ch0=2997.;
    para="fid";

    relpara();
  }

  void argument(int argc,char*argv[]) {
    char c;
    double dy;
    optind=1;
    while((char)EOF!=(c=getopt(argc,argv,"b:c:k:s:t:d:a:f:r:w:z:o:i:W:S:")))
      switch(c) {
      case 'b': 
	ombh2+=atof(optarg); 
	para="ombh2";
	ifid=0;
	break;
      case 'c': 
	omch2+=atof(optarg); 
	para="omch2";
	ifid=0;
	break;
      case 'k': 
	omegak+=atof(optarg);
	para="omk";
	ifid=0;
	break;
      case 's': 
	dy=atof(optarg); 
	para="ns";
	ns+=dy;
	ampk*=pow(k0/0.05,-dy);
	ifid=0;
	break;
      case 't': 
	tau+=atof(optarg); 
	para="tau";
	ifid=0;
	break;
      case 'd': 
	omegad+=atof(optarg);
	para="omde";
	ifid=0;
	break;
      case 'a': 
	ampk*=(1+atof(optarg)); 
	para="amp";
	ifid=0;
	break;
      case 'f': 
	omnuh2+=(ombh2+omch2+omnuh2)*atof(optarg); 
	para="nu1";
	nu_l-=1.;
	nu_m+=1.;
	ifid=0;
	break;
      case 'r': 
	nrun+=atof(optarg); 
	para="nrun";
	ns-=nrun*log(k0/0.05);
	ifid=0;
	break;
      case 'w': 
	w0+=atof(optarg); 
	para="w0";
	ifid=0;
	break;
      case 'z': 
	break;
      case 'o': 
	break;
      case 'i': 
	infpkl=optarg; 
	break;
      case 'W': 
	para="off";
	ifid=0;
	break;
      case 'S':
	break;
      }
    relpara();
  }

  void relpara() {
    omegac=omch2/h0/h0;
    omegab=ombh2/h0/h0;
    omeganu=omnuh2/h0/h0;
    omegam=1-omegak-omegad;
    h0=sqrt((ombh2+omch2+omnuh2)/omegam);
    rho0=rhoc*omegam;
    //!!! Peacock & Dodds 94, Sugiyama 95
    gamma=omegam*h0*exp(-omegab*(1+sqrt(2*h0)/omegam));
  }

  void outpara() {
    cout << "ombh2 " << ombh2 << endl;
    cout << "omch2 " << omch2 << endl;
    cout << "omegam " << omegam << endl;
    cout << "ns " << ns << endl;
    cout << "tau " << tau << endl;
    cout << "ampk " << ampk << endl;
    cout << "omnuh2 " << omnuh2 << endl;
    cout << "nu_l " << nu_l << endl;
    cout << "nu_m " << nu_m << endl;
    cout << "nrun " << nrun << endl;
    cout << "w0 " << w0 << endl;
    cout << "wa " << wa << endl;
    cout << "omegad " << omegad << endl;
    cout << "omegak " << omegak << endl;
    cout << "h0 " << h0 << endl;
    cout << "omegac " << omegac << endl;
    cout << "omeganu " << omeganu << endl;
    cout << "infpkl " << infpkl << endl;
  }

};

#endif // CPARA_H
