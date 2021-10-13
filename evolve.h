//evolve.h
#ifndef EVOLVE_H
#define EVOLVE_H
#define pi 3.141592653589793238462643

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include"cpara.h"

namespace evolve {

  double mnfw(const double c);
  double mnfw2(const double c);
  double rs(const Cpara& cp, const double& m, const double& z);
  double rvir(const Cpara& cp, const double& m, const double& z);
  double delv(const Cpara& cp, const double& z);
  double rhos(const Cpara& cp, const double& m, const double& z);
  double sigv(const Cpara& cp, const double& x, const double& m, const double& z);
  double conpara(const double& m, const double& z);
  double sigv_vir(const Cpara& cp, const double& m,const double& z);
  double growth_deriv(const Cpara& cp, const double& z);
  double growth_deriv_full(const Cpara& cp, const double& z);
  double growth(const Cpara& cp, const double& z);
  double growth_full(const Cpara& cp, const double& z);
  double omegaf(const Cpara& cp, const double& z);
  double omegamz(const Cpara& cp, const double& z);
  double omegakz(const Cpara& cp, const double& z);
  double omegadz(const Cpara& cp, const double& z);
  double hrate(const Cpara& cp, const double& z);
  double hratederiv(const Cpara& cp, const double& z);
  double zvol(const Cpara& cp, const double& zmin, const double& zmax, const double& omega);
  double angdis(const Cpara& cp, const double& z);
  double* angdis_out(const Cpara& cp, const double& zmax, const int& nzmax);
}

namespace evolve {

  double mnfw(const double c) {
    return log(1+c)-c/(1+c);
  }

  double mnfw2(const double c) {
    return 1-log(1+c)/c;
  }

  double rs(const Cpara& cp, const double& m, const double& z) {
    //unit: Mpc/h
    return rvir(cp,m,z)/conpara(m,z);
  }

  double rvir(const Cpara& cp, const double& m, const double& z) {
    //virial radius at redshift z in unit of Mpc/h
    return pow(3.*m/4./pi/delv(cp,z)/cp.rho0,1./3.)/(1+z);
  }

  double delv(const Cpara& cp, const double& z) {
    // over density of a virialized object at redshift z
    double wf=1/omegaf(cp,z)-1;
    return 18.*pi*pi*(1+0.4093*pow(wf,0.9052));
  }

  double rhos(const Cpara& cp, const double& m, const double& z) {
    //unit: M_solar * h * h
    return m/(4*pi*pow(rs(cp,m,z),3)*mnfw(conpara(m,z)));
  }

  double sigv(const Cpara& cp, const double& x, const double& m, const double& z) {
    // x in unit of rs;
    double c=conpara(m,z);
    return sigv_vir(cp,m,z)*sqrt((mnfw(x)/x)/(mnfw(c)/c));
  }

  double conpara(const double& m, const double& z) {
    // unit of m : h^{-1}*m_solar
    // ref. Duffy, Schaye, Kay, dalla Vecchia, MNRAS 390, L64, 2008
    double afac=7.85,bfac=-0.081,cfac=-0.71;
    return afac*pow(m/2.e12,bfac)*pow(1+z,cfac);
  }

  double sigv_vir(const Cpara& cp, const double& m,const double& z) {
    // unit: km/s
    return sqrt(cp.gconst*m/2./rvir(cp,m,z));
  }

  double growth_deriv(const Cpara& cp, const double& z) {
    // derivative of linear growth rate  dln D/dln a 
    // a fitting function by Lahav et al. 1991
    double ommz = omegamz(cp,z)/pow(hrate(cp,z),2);
    double omdz = omegadz(cp,z)/pow(hrate(cp,z),2);
    return pow(ommz,4./7.)+(1+ommz/2.)*omdz/70.;
  }

  double growth_deriv_full(const Cpara& cp, const double& z) {
    // derivative of linear growth rate  dln D/dln a 
    // full calculation
    double ommz = omegamz(cp,z)/pow(hrate(cp,z),2);
    double omdz = omegadz(cp,z)/pow(hrate(cp,z),2);
    return -1-0.5*ommz-0.5*omdz*(1+3*(cp.w0+cp.wa*z/(1+z)))
      +2.5*ommz/(growth_full(cp,z)*(1+z));
  }

  double growth(const Cpara& cp, const double& z) {
    // linear growth rate  D(z)
    // fitting function by Carrol Press & Turner 1992, Peacock and Dodds. 1996
    double ommz = omegamz(cp,z)/pow(hrate(cp,z),2);
    double omdz = omegadz(cp,z)/pow(hrate(cp,z),2);
    double gz = 2.5*ommz/(pow(ommz,4./7.)-omdz+(1+ommz/2.)*(1+omdz/70.)); 
    return gz/(1+z);
  }

  double growth_full(const Cpara& cp, const double& z) {
    // linear growth rate D(z), full calculation
    const int imax=1000;
    double a=1/(1+z),da=a/imax,gz=0;
    for (int i=0;i<imax;i++) {
      double at=da*(i+0.5);
      double zt=1/at-1;
      gz+=da/pow(at*hrate(cp,zt),3);
    }
    gz*=2.5*cp.omegam*hrate(cp,z);
    return gz;
  }

  double omegaf(const Cpara& cp, const double& z) {
    return omegamz(cp,z)/(omegamz(cp,z)+omegakz(cp,z)+omegadz(cp,z));
  }

  double omegamz(const Cpara& cp, const double& z) {
    // Omegam at z
    return cp.omegam*pow(1+z,3);
  }

  double omegakz(const Cpara& cp, const double& z) {
    // Omegak at z
    return (1-cp.omegam-cp.omegad)*pow(1+z,2);
  }

  double omegadz(const Cpara& cp, const double& z) {
    // omegadz=omegad*(1+z)**(3*(1+w(z)))
    // w(z)=w0+wa*(1-a)
    return cp.omegad*pow(1+z,3*(1+cp.w0+cp.wa))*exp(-3*cp.wa*z/(1+z));
  }
  
  double hrate(const Cpara& cp, const double& z) {
    // H(z)/H0: Hubble parameter normalized to unity at z=0;
    return sqrt(omegamz(cp,z)+omegakz(cp,z)+omegadz(cp,z));
  }

  double hratederiv(const Cpara& cp, const double& z) {
    // dln H(a)/dln a: derivative of Hubble parameter by a 
    return sqrt(omegamz(cp,z)+omegakz(cp,z)+omegadz(cp,z));
  }
  
  double zvol(const Cpara& cp, const double& zmin, const double& zmax, const double& omega) {
    // unit: Gpc/h;
    const int imax=1000;
    double vol=0,dz=(zmax-zmin)/imax,r=angdis(cp,zmin)*1e-3;
    for (int i=0;i<imax;i++) {
      double z=zmin+dz*(i+0.5);
      double dr=3*dz/hrate(cp,z);
      vol+=omega*r*r*dr;
      r+=dr;
    }
    return vol;
  }

  double angdis(const Cpara& cp, const double& z) {
    //  angular diameter distance in unit of Gpc/h
    const int imax=1000;
    double r=0.,dz=z/imax;
    for (int i=0;i<imax;i++) {
      double zt=dz*(i+0.5);
      r+=2.997e3*dz/hrate(cp,zt);
    }
    return r;
  }

  double* angdis_out(const Cpara& cp, const double& zmax, const int& nzmax) {
    double *dis=new double[nzmax];
    for (int i=0;i<nzmax;i++) {
      double z=zmax/nzmax*(i-0.5);
      dis[i]=angdis(cp,z);
    }
    return dis;
  }

};
#endif // EVOLVE_H

