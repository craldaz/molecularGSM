#ifndef TC_H
#define TC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <math.h>

#include "stringtools.h"
#include "pTable.h"
#include "constants.h"
#include "utils.h"

class TC
{
  private:
  
  int nscffail;
  int firstrun;

   int runNum;
   int runend;
   string outfile;
   string scrdir;
   string scrBaseDir;
   string runName;
   string runName0;
   string fileloc;
   double* xyz;
    double* grad;
    double* dvec;
    double* dgrad;

   int natoms;
   int* anumbers;
   string* anames;
   string basis;

   int nclosed;
   int nstates;
   int nactive;
   double* E;

  public:

   int CHARGE;
   int MULT;
   int wstate; 
	 int wstate2; 
	 int wstate3; 
   string method; //for TeraChem
   bool docoupling;
   bool readOrb;

   void grads(int runend,int runnum);
   void calc_dvec(int runend,int runnum);
   void get_dvec(int runend,int runnum,double* dvec);
   double calc_energy();
   void alloc(int natoms);
   void init(string infilename, int natoms, int* anumbers, string* anames, int run, int rune);
   void init(int natoms0,int nstates0, int nclosed0, int nactive0, string* anames0, double* xyz0, string basis0,string method0,int run, int rune);
   void freemem();
   void write_xyz_grad(double* coords, double* grad, string filename);

   int ncpu;
   int gradcalls;
   void reset(double* xyz1);
   double energy0;
   double energy;
   double get_grad(int choosestate,int runend,int runnum, double* grad);
   double getE(int choosestate);
   int readE();

};

#endif
