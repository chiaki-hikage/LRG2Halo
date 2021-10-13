//main.h
#ifndef MAIN_H
#define MAIN_H
#define pi 3.141592653589793238462643
#define degrad pi/180

#include<iostream>
#include<math.h>
#include<fstream>
#include<stdlib.h>
#include<string.h>
using namespace std;

class Parlist {

 public:
  
  Parlist(){}
  ~Parlist(){}

  double len;          // side lengh
  double zmin,zmax;    // redshift range
  int npix,njack;        // grid number of each side
  char *outf,*outfwin;         // output filename
  char *infd, *infr;       // filename of input image (ascii)
  double kmin, kmax; // kmin, kmax in unit of knyq
  int nkbin;       // binning number
  char *type, *space, *dir;     // data or random

  /* set default values for global variables */				   
  void defaultval(){
    len=3072.;               
    npix=512;
    outf="pow_1d.dat";
    outfwin="pow_win.dat";
    infd="data/DR7-Full.ascii";
    infr="data/random-DR7-Full.ascii";
    type="m";
    space="r";
    kmin=0.002;
    kmax=0.6;
    nkbin=300;
    dir="n";
    zmin=0.16;
    zmax=0.44;
    njack=9999;
  }
  
  /* random command line into the global parameter list */
  void argument(int argc,char*argv[]) {
    char c;
   
    while((char)EOF!=(c=getopt(argc,argv,"l:p:o:w:i:r:t:s:d:b:z:Z:j:"))) 
      switch(c) {
      case 'l': len      = atof(optarg); break;
      case 'p': npix     = atoi(optarg); break;
      case 'b': nkbin    = atoi(optarg); break;
      case 'o': outf     =      optarg ; break;
      case 'w': outfwin  =      optarg ; break;
      case 'i': infd     =      optarg ; break;
      case 'r': infr     =      optarg ; break;
      case 't': type     =      optarg ; break;
      case 's': space    =      optarg ; break;
      case 'd': dir      =      optarg ; break;
      case 'z': zmin     = atof(optarg); break;
      case 'Z': zmax     = atof(optarg); break;
      case 'j': njack    = atoi(optarg); break;
      }
  }
  
};

#undef degrad
#undef  pi
#endif // MAIN_H
