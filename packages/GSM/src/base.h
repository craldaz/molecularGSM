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

#include "GitSHA1.h"
#include "constants.h"
//#include "mecp.h"
#include "gstring.h"

using namespace std;

// This file is the new "driver" class

class Base 
{
  private:
     //MECP mecp;
    //MECI meci;
    //GString gstr;
    

   // ==> String parameters
      int nnnax;
      int runNum;
      int runend;
      string runends;
      string infile0;
      int nn;
      int nnR;
      int nnP;
      int nnmax;
      int nnmax0; //input value of nnmax
      int ncpu;

   // ==> Method parameters
      int isMECP;


  public:
    void driver();
    void init(string infilename, int run, int nprocs);
    int isomer_init(string isofilename);
    void parameter_init(string infilename);
};

#endif
