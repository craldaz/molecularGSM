#include "mecp.h"


void Mecp::calc_V0()
{
  //printf("Hello\n");
    V0 = icoords[0].grad1.grads(coords[0], grads[0], icoords[0].Ut, 3);
  printf("  setting V0 to: %8.1f (%12.8f au) \n",V0,V0/627.5);
}

void Mecp::driver()
{
    //preopt if turned on
   
    for (int n=0;n<nnmax0;n++)
        icoords[n].grad_init(ncpu,runNum,n,0,0); 
 
    calc_V0(); 

   double* dq = new double[len_d+100];
   //TODO addNode
   //addNode(0,...);
   //TODO print_string(...); //put in base class
}

Mecp::Mecp(int run, int nprocs,int NNODES)/*, float DQMAG_SSM_MAX, float DQMAG_SSM_MIN,
          int MAX_OPT_ITERS, int STEP_OPT_ITERS, int INITIAL_OPT,
          float SSM_DQMAX, float CONV_TOL, float ADD_NODE_TOL,
          int BOND_FRAGMENTS, int NNODES)*/ //TODO put parameters here
{ 
  
  //read xyz file
  isSSM=1;
  string nstr=StringTools::int2str(run,4,"0");
  ncpu = nprocs;
  runNum = run;
  runends = nstr;
  string xyzfile = "scratch/initial"+nstr+".xyz";
  printf(" %i\n",isSSM);
  structure_init(xyzfile);
 
  //initialize ISOMER if necessary
  string isomerfile = "ISOMERS"+nstr;
  struct stat sts;
  if (stat(isomerfile.c_str(), &sts) != -1)
    printf(" using ISOMERS file: %s \n",isomerfile.c_str());
  else
    isomerfile = "scratch/ISOMERS"+nstr;
  int nfound = isomer_init(isomerfile);

  isomer_init(isomerfile);
  initialize_icoords();
}

