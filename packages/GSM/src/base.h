#ifndef BASE_H
#define BASE_H


#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <sys/stat.h>

#include "GitSHA1.h"
#include "icoord.h"
#include "utils.h"
#include "stringtools.h"
#include "pTable.h"
#include "constants.h"

using namespace std;


class Base 
{
  protected:
     //MECP mecp;
    //MECI meci;
    //GString gstr;
    

   // ==> String parameters
      int nnnax;
      int runNum;
      int runend;
      string runends;
      int nn;
      int nnR;
      int nnP;
      int nnmax;
      int ncpu;
  int tstype;
  double prodelim;
  int lastOpt;
  int initialOpt;
  double DQMAG_SSM_MAX;
  double DQMAG_SSM_MIN;
  double QDISTMAX;
  double PEAK4_EDIFF;
  double SCALING;
  int GROWD;


   // ==> Method parameters
      int isMECP;
  int isRestart;
  int hessSSM; //starting SSM hessian given
  int isFSM; //freezing string flag
  int isMECI; //MECI opt
  int isPRODUCT; //MECI opt
        int isMAP_DE;
        int isMAP_SE;
        int isSE_ESSM;
  int use_exact_climb; //whether to climb or do exact TS search
        int restart_wfn;                                                //don't seed MOLPRO wfn, use existing files
  int STEP_OPT_ITERS;
  int MAX_OPT_ITERS;
  double CONV_TOL;
  double ADD_NODE_TOL;
  int bondfrags;

 // ==> structure_init
  int natoms;
  int CHARGE;                   //charge of the molecular complex
  int* anumbers;                //array of atomic indices (for looking up period table stuff)
  double* amasses;              //array of atomic masses (used for mass-weighting coordinates)
  string* anames;               //array of atomic symbols (for creating input QC file)
  int* frozen;
  double** coords;
  double** tangents;
  double** grads;
  double** perp_grads;
  double* V_profile;

 // isomer_init
  int nfound;
  int nbond;
  int nadd;
  int nbrk;
  int nangle;
  int ntors;
  int* bond;
  int* add;
  int* brk;
  int* angles;
  double* anglet;
  int* tors;
  double* tort;




  public:
    int nnmax0; //input value of nnmax
    int isomer_init(string isofilename);
    void parameter_init(string infilename);
 //**** Initialize functions *****
  void structure_init(string xyzfile);
  void initialize_icoords();
  ICoord* icoords;
  int isSSM; //shooting string flag
};

#endif