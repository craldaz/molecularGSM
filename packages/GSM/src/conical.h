#ifndef CONICAL_H
#define CONICAL_H

#include "stringtools.h"
#include "pTable.h"
#include "grad.h"
#include "icoord.h"


class Conical {
	private:
	 ICoord* ics;
	 double** coords;	
   double** dgrada; 
   double** dgrada_q; 
   double** dgrada_U; 
   double** dveca; 
   double** dveca_q; 
   double** dveca_U; 

   double** grada; 

	 string* anames;
  int* anumbers;		//array of atomic indices (for looking up period table stuff)
	 int nnodes;
	int n_cpu;
	int runEnd;
	int run;
	int natoms;
	int len_d;
	int size_ic;
	double dE;
	int wstate;
	int wstate2;

	public:
	 Conical(int nnodes, ICoord* icoords,int ncpu, int runNum, int runend); //constructor
	 void print_bp();
	 void print_xyz();
	 void opt_meci(int node,int runNum);
};
#endif
