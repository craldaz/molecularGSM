#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <sys/stat.h>
    
#include "utils.h"
#include "mopac.h"
#include "stringtools.h"
#include "align.h"   
     
      
using namespace std;
      
void xyz_read(int natoms, string* anames, double* coords, string xyzfile);
void xyz_read_last(int natoms1, double* coords, string xyzfile);
int get_natoms(string filename);
int get_charge(string filename);
void get_all_xyz(int natoms, string* anames, vector<double*> &xyzs, string xyzfile);
void align_and_opt(int natoms1, int natoms2, string* anames, string* anamesm, string* anamest, int* anumbers, int* anumbersm, int charget, int nstruct, int* unique, vector<double*> xyzall, double* xyzm);
void write_all_xyz(int natoms, string* anames, double* E, vector<double*> xyzs, string xyzfile_string);
void write_all_xyz(int natoms, string* anames, int nstruct, double* E, double** xyzs, string xyzfile_string);
void write_gsm(int natoms, string* anames, int charge, int nstruct, double* E, double** xyzs, int nadd, int* adds);
void do_gsm(int nstruct);


