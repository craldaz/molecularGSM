#include "conical.h"
#include "utils.h"
#define DG_ROT 1
#define UPDATE_BP 0

using namespace std;


Conical::Conical(int nnodes_value, ICoord icoord, int ncpu, int runNum, int runend,int nsteps,int isMECI) : nnodes(nnodes_value), n_cpu(ncpu), run(runNum), runEnd(runend), opt_steps(nsteps)
{
	printf(" Default constructor\n");

	natoms = icoord.natoms;
	conical.alloc(icoord.natoms);
	conical.reset(natoms,icoord.anames,icoord.anumbers,icoord.coords);
	conical.copy_ic(icoord);
	
	//allocate the bmatrix and gradient
	string infilename="inpfileq";
	conical.bmat_alloc();
	conical.bmatp_create();
	conical.bmatp_to_U();
	conical.bmat_create();
	if (isMECI)
		conical.grad1.seedType = 3;
	else 
		conical.grad1.seedType = icoord.grad1.seedType;
 	conical.grad_init(infilename,n_cpu,run,runEnd,0,0);
	
	grad = new double[3*natoms];
	dgrad = new double[3*natoms];
	dvec = new double[3*natoms];
	coords = new double[3*natoms];

 for (int i=0;i<3*natoms;i++)
   coords[i]=icoord.coords[i];


}

void Conical::calc_dgrad()
{
	printf(" Calculating dgrad\n");
	energy=conical.grad1.grads(conical.coords, grad,conical.Ut, 3); 
	for (int i=0;i<3*natoms;i++)
  	dgrad[i] = conical.grad1.grada[1][i] - conical.grad1.grada[0][i]; 
	return;
}
void Conical::calc_dvec()
{
	printf(" Calculating dvec\n");
	conical.grad1.dvec_calc(conical.coords, dvec,run,runEnd); 
	return;
}

double Conical::calc_BP()
{
	calc_dgrad();
	calc_dvec();
	return energy;
}

void Conical::print_bp()
{
  printf(" printing dgrad\n");
  for (int i=0;i<3*natoms;i++)
    printf("%1.3f\t",dgrad[i]); 
  printf("\n");
  printf(" printing dvec\n");
  for (int i=0;i<3*natoms;i++)
    printf("%1.3f\t",dvec[i]);  
  printf("\n");
	
}

double Conical::form_MECI_space()
{

	energy = calc_BP();
	//print_bp();
	printf(" Forming the 3N-6 dimensional space defined by a MECI\n");
	conical.dgrad_to_dgradq(dgrad);
	conical.dvec_to_dvecq(dvec);

#if DG_ROT
	double norm_dg=conical.dgrot_mag();
	conical.project_dgradq();
	conical.project_dvecq();
#else
	//need to code dvec rot
	//project(dvecq,dvecq_U);
	//norm_dg=project(dgradq,dgradq_U);
	printf(" norm_dg = %1.2f",norm_dg); 
#endif
	conical.constrain_bp();
	conical.bmat_create();
	//conical.print_q();

	return energy;
}


void Conical::opt_meci()
{
	printf(" Optimizing node to MECI using Combined-Step Optimizer\n");
	energy = form_MECI_space();	
	printf(" Initial energy is %1.4f\n",energy);
	conical.V0 = energy;
   
	string runName0 = StringTools::int2str(run,4,"0")+"."+StringTools::int2str(runEnd,4,"0");	
	//use combined step optimizer
	conical.gradrms = 100.;
	conical.make_Hint();
	conical.opt_MECI("scratch/cfile_"+runName0+".xyz",opt_steps,runEnd,run,grad,dvec,dgrad);
	
  printf(" %s",conical.printout.c_str());
 	for (int i=0;i<3*natoms;i++)
 	  coords[i]=conical.coords[i];
	
	return;
}
