#ifndef GRAD_H
#define GRAD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <math.h>

#include "stringtools.h"
#include "qchem.h"
#include "gaussian.h"
#include "ase.h"
#include "knnr.h"
#include "orca.h"
#include "molpro.h"
#include "terachem.h"
#include "qchemsf.h"

class Gradient 
{
  
  private:
  
   int runNum; //the job number
   int runend; //the node number
   string runends;
   string runName;
   string runName0;

   int natoms;
   int N3;
   string* anames;
   int* anumbers;
   int CHARGE; //total system charge
   int MULT; //multiplicity

   QChem qchem1;
   QChemSF qchemsf1;
   GAUSSIAN gaus1;
   ASE ase1;
   ORCA orca1;

   int knn_k;
   KNNR knnr1;
   int knnr_inited;

   int nforce;
   int* fa;
   double* fv;
   double* fk;

   int read_nstates();
   void read_molpro_settings(int& nstates, int& nclosed, int& nocc, int& nelec, string& basis);
   int read_molpro_init(string* &hf_lines);
   void read_tc_settings(int& nstates0, int& nclosed, int& nactive, string& basis,string& method);
   int force_init(string ffile);
   
   void read_qchem_settings(int& sstates, int& tstates);

  public:

   int hessian(double* H);
   double grads(double* coords, double* grad, double* Ut, int type);
   void add_force(double* coords, double* grad);
   void init(int natoms, int* anumbers, string* anames, double* coords, int run, int rune, int ncpu, int use_knnr, int q1);
   void update_knnr();
   void freemem();
   void write_xyz_grad(double* coords, double* grad, string filename);
   int external_grad(double* coords, double* grads);

   int knnr_active;
   int always_do_exact;
   int write_on;
   int wrote_grad;
   int xyz_grad;
   int gradcalls;
   int nscffail;
   double V0;
   double fdE; //force * distance energy

   double energy0;
   double energy;

  //for multistate
   int nstates;
   int wstate;
   int wstate2;
   int wstate3;
   int sstates;
   int tstates;
   int swstate;
   int swstate2;
   int swstate3;
   int twstate;
   int twstate2;
   int twstate3;


	 double* dE;
   double** grada; //multistate gradients
   double* E; //for multiple states
   Molpro mp1;
   TC tc1;

		double dDmw_dR;

   int seedType;
	 int dvec_calc(double* coords, double* dvec, int run, int rune);
   double energy_initial(double* coords,int run, int rune,int penalty,double sigma);
	double levine_penalty(double* coords, double* grad, double* Ut, int type,double sigma,double dmw);
   
   int res_t; //restart found files

};

#endif
