#include "conical.h"
#include "utils.h"

using namespace std;


Conical::Conical(int nnodes_value, ICoord* icoords) : nnodes(nnodes_value)
{
	printf(" Default constructor\n");
	//dynamic memopry allocation
	
	natoms = icoords[0].natoms;
	len_d = icoords[0].nicd0;
	size_ic = icoords[0].nbonds + icoords[0].nangles+ icoords[0].ntor;

	//Branching plane (BP) in cartesian coordinates
	dgrada = new double*[nnodes]; 
  for (int i=0;i<nnodes;i++)
      dgrada[i] = new double[3*natoms];                                                   
  for (int n=0;n<nnodes;n++)
    for (int i=0;i<3*natoms;i++)                                                          
      dgrada[n][i] = 0.; 
  dveca = new double*[nnodes];                                                    
  for (int i=0;i<nnodes;i++)                                                               
      dveca[i] = new double[3*natoms];                                                    
  for (int n=0;n<nnodes;n++)                                                               
    for (int i=0;i<3*natoms;i++)
      dveca[n][i] = 0.;

	//Branching plane (BP) in representation of delocalized ICs (U)
  dgrada_q = new double*[nnodes];                                                  
  for (int i=0;i<nnodes;i++)
      dgrada_q[i] = new double[len_d+100];                                                 
  for (int n=0;n<nnodes;n++) 
    for (int i=0;i<len_d+100;i++)
      dgrada_q[n][i]=0.;
  dveca_q = new double*[nnodes];                                                   
  for (int i=0;i<nnodes;i++)
      dveca_q[i] = new double[3*natoms];
  for (int n=0;n<nnodes;n++)
    for (int i=0;i<3*natoms;i++)
      dveca_q[n][i] = 0.;

	//Branching plane (BP) delocalized vectors in representation of ICs
  dveca_U = new double*[nnodes];
  for (int i=0;i<nnodes;i++)
      dveca_U[i] = new double[size_ic+100];
  for (int n=0;n<nnodes;n++)
    for (int i=0;i<size_ic;i++)
      dveca_U[n][i] = 0.;
	dgrada_U = new double*[nnodes];
  for (int i=0;i<nnodes;i++)
      dgrada_U[i] = new double[size_ic+100];
  for (int n=0;n<nnodes;n++)
    for (int i=0;i<size_ic;i++)
      dgrada_U[n][i]=0.;

	coords = new double*[nnodes];
	for (int i=0;i<nnodes;i++)
		coords[i] = new double[3*natoms];
	for (int n=0;n<nnodes;n++)
		for (int i=0;i<3*natoms;i++)
			coords[n][i] = icoords[n].coords[i];	

	anames = new string[natoms];
	for (int i=0;i<natoms;i++)
		anames[i]=icoords[0].anames[i];
			
}


void Conical::print_bp()
{
  printf(" printing dgrad\n");
	for (int n=0;n<nnodes;n++)
	{
  	for (int i=0;i<3*natoms;i++)
  	  printf("%1.3f\t",dgrada[n][i]); 
  	printf("\n");
	}
  printf(" printing dvec\n");
	for (int n=0;n<nnodes;n++)
	{
  	for (int i=0;i<3*natoms;i++)
  	  printf("%1.3f\t",dveca[n][i]);  
  	printf("\n");
	}

}


void Conical::print_xyz()
{
	printf(" %i \n",natoms);
	printf("\n");
	for (int n=0;n<nnodes;n++)
	{
		for (int i=0;i<natoms;i++) 
		{
				cout << "  " << anames[i];
				printf(" %f %f %f \n",coords[n][3*i+0],coords[n][3*i+1],coords[n][3*i+2]);
		}
		printf("\n");
	}
}
