#ifndef CONICAL_H
#define CONICAL_H

#include "stringtools.h"
#include "pTable.h"
#include "grad.h"
#include "icoord.h"


class Conical {
	private:
	 ICoord conical;

   double* dgrad; 
   //double* dgrad_q; 
   //double* dgrad_U; 
   double* dvec;
   //double* dvec_q; 
   //double* dvec_U; 
	 double* grad; 

	double energy; 

	int nnodes;
	int n_cpu;
	int runEnd;
	int run;
	int natoms;

	double dE;
	int wstate;
	int wstate2;
	int nstates;
	int opt_steps;
	//void form_MECI_space(int node);
	//void constrain_bp(double* dgrad_U,double* dvec_U,int node);
	//void project(double* gradq,double* gradq_U);
	//void vec_to_q(double* grad1,double* gradq1);
	//double dgrot_mag(double* dgradq,double* dvecq);

	public:
	 Conical(int nnodes, ICoord icoord,int ncpu, int runNum, int runend,int nsteps,int isMECI); //constructs the object
	 void print_bp();
	 //void print_xyz();
	 void opt_meci(); // function optimizes geometry to meci
	 void calc_dgrad();
	 void calc_dvec();
	 double calc_BP();
	 double form_MECI_space();
	 double* coords;	
};

#endif
