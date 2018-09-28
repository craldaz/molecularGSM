#include "mecp.h"


void Mecp::calc_V0()
{
  printf("Hello\n");
    V0 = icoords[0].grad1.grads(coords[0], grads[0], icoords[0].Ut, 3);
  printf("  setting V0 to: %8.1f (%12.8f au) \n",V0,V0/627.5);
}

void Mecp::driver()
{
    //preopt if turned on
    init()
    for (int n=0;n<nnmax0;n++)
    {
 	    icoords[n].grad_init(infile0,ncpu,runNum,n,0,0); //TODO take infile0 out of everything called by grad_init
    }
    calc_V0(); //TODO eventually we need to both singlet and triplet (this will be calculator dependent, e.g. qchem)
    //if MECP 
      //set variables, parameters, etc
      //do opt,etc

}

Mecp::Mecp(int run, int nprocs) //TODO put parameters here
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

