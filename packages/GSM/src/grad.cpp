#include "grad.h"
#include "mopac.h"

using namespace std;

#define NOGRAD 0


void Gradient::add_force(double* coords, double* grad)
{
  fdE = 0.;
 //add spring forces to bonds
  for (int i=0;i<nforce;i++)
  {
    int a1 = fa[2*i+0];
    int a2 = fa[2*i+1];
   //keep a.u. for now
    double dx = (coords[3*a1+0] - coords[3*a2+0])*ANGtoBOHR;
    double dy = (coords[3*a1+1] - coords[3*a2+1])*ANGtoBOHR;
    double dz = (coords[3*a1+2] - coords[3*a2+2])*ANGtoBOHR;
    double d = sqrt(dx*dx+dy*dy+dz*dz);
    double d0 = fv[i] * ANGtoBOHR;

   //spring force
    double t = (1-d0/d) * fk[i];
   //for "zero" equilibrium distance, use constant force
    if (fv[i]<=0.) t = fk[i] / d / 2.;  
    fdE += fk[i] * d;
    printf(" d: %8.5f t: %8.5f fk: %8.5f fk*d: %8.5f \n",d,t,fk[i],fdE);

    t *= ANGtoBOHR; //back to HT/A
    if (grad!=NULL)
    {
      grad[3*a1+0] += t*dx;
      grad[3*a1+1] += t*dy;
      grad[3*a1+2] += t*dz;
      grad[3*a2+0] -= t*dx;
      grad[3*a2+1] -= t*dy;
      grad[3*a2+2] -= t*dz;
      printf("f");
    }
  }

  return;
}

int Gradient::hessian(double* H)
{
  printf("\n ERROR: entering untested function \n");
  exit(1);
#if QCHEM
  return qchem1.read_hess(H);
#else
  return 0;
#endif
}

void Gradient::update_knnr()
{
#if !USE_KNNR
  return;
#endif
  if (knnr1.npts>0 && knnr_inited)
    int newpts = knnr1.add_extra_points();
  else
  {
    printf("  loading KNNR files \n");
    int kfound = knnr1.begin(runNum,natoms);
    if (kfound>0)
      res_t = 1000;
    knnr_inited = 1;
  }
  return;
}


//type 0: QM gradient only
//type 1: knnr if accurate
//type 2: knnr only
//type 3: knnr if exact

double Gradient::grads(double* coords, double* grad, double* Ut, int type)
{
  //printf(" grads(%i) k: %i \n",type,knn_k);
#if USE_KNNR
  if (type==2 && knnr_active==0)
  {
    printf("\n ERROR: grad type 2 (knnr only) while knnr_active==0 \n");
    exit(1);
  }
#endif

  if (knnr_active && knnr_inited==0 && type>0) update_knnr();

  if (knnr_active==3) type = 3; //force type 3

  wrote_grad = 0;

  double errknn = 10.;
#if USE_KNNR
  if (knnr_active && type>0)
  {
    //printf(" about to call grad_knnr() \n");
    errknn = knnr1.grad_knnr(coords,energy,grad,Ut,knn_k);
    energy *= 627.5;
    //printf(" E(kNN): %4.3f",energy);
    xyz_grad = 0;
    if (energy - V0 < -500. || energy - V0 > 500.)
    {
      printf(" kf");
      errknn = 99.;
    }
  }
#else
  energy = -99.;
  type = 0;
#endif

  //if (type==3) printf(" t3errknn: %5.4f",errknn);

  int do_exact = 0;
  if (type==0) do_exact = 1;
  else if (type==1 && errknn > KNNR_MAX_DIST) do_exact = 1;
  else if (type==2) do_exact = 0;
  else if (type==3 && errknn > 0.01) do_exact = 1;

  if (always_do_exact) do_exact = 1;
  //do_exact = 1; //debug
  if (do_exact)
  {
    printf(" eg"); fflush(stdout);
#if !NOGRAD
    int success = external_grad(coords,grad);
#endif
    xyz_grad = 1;
  }
  else
    printf(" kg"); 

  add_force(coords,grad);
  //printf(" E: %12.8f",energy);
  energy += fdE*627.5;
  //printf(" E+fdE: %12.8f \n",energy);

#if 0
  if (xyz_grad)
  {
    printf(" Grad: \n");
    for (int i=0;i<natoms;i++)
      printf(" %s %12.10f %12.10f %12.10f \n",anames[i].c_str(),grad[3*i+0],grad[3*i+1],grad[3*i+2]);
  }
#endif

#if 0
  printf(" XYZ: \n");
  for (int i=0;i<natoms;i++)
    printf(" %s %4.3f %4.3f %4.3f \n",anames[i].c_str(),coords[3*i+0],coords[3*i+1],coords[3*i+2]);
#endif

  return energy;
}

int Gradient::external_grad(double* coords, double* grad)
{

#if QCHEM
  energy = qchem1.grads(coords,grad);
#elif QCHEMSF
	if (wstate<0)
	{
		printf("wstate must be greater than 1!\n"); 
		exit(-1); 
	} 
	else
	{
  	qchemsf1.calc_grads(coords);
  	energy = qchemsf1.getE(wstate-1);
  	if (V0==0.) V0 = energy;
  	qchemsf1.getGrad(wstate,grada[0]);
  	for (int i=0;i<N3;i++)
  	  grad[i] = grada[0][i];
	}
	if (wstate2>0) 
   {
    energy += qchemsf1.getE(wstate2-1);
    qchemsf1.getGrad(wstate2,grada[1]);
    if (wstate3==0)
    	{
				energy /= 2;
    		for (int i=0;i<N3;i++) 
      		grad[i] = (grada[0][i] + grada[1][i])/2.;
			} 
  }
  if (wstate3>0)
  {
    energy += qchemsf1.getE(wstate3);
    qchemsf1.getGrad(2,grada[2]);
		energy /=3.;
    for (int i=0;i<N3;i++) 
      grad[i] = (grada[0][i] + grada[1][i] + grada[2][i])/3.;
  }
  for (int i=0;i<nstates;i++)
    E[i] = qchemsf1.getE(i);
	
	//for (int i=0;i<nstates;i++)
	//	printf(" E[%i]=%1.4f\n",i,E[i]);
	
#elif USE_MOLPRO
  mp1.reset(coords);
	//cout << " seedType is "<< seedType << endl;
  if (gradcalls==0 && seedType<1)
  {
    //copy wfn because hasn't been set yet
    int runendCopy = runend;
    if (seedType==0) 
    {
      printf(" ERROR: seedType not set! \n");
      exit(1);
    }
    if (seedType==-1) runendCopy--;
    if (seedType==-2) runendCopy++;

		if (seedType==-1 || seedType==-2)
		{
    	string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
    	string cmd = "cp scratch/mp_"+runNameCopy+" scratch/mp_"+runName0;
    	printf(" copying file from mp_%s to mp_%s \n",runNameCopy.c_str(),runName0.c_str());
    	system(cmd.c_str());
		}
  }
  int error = mp1.run(); //grad and energy
	if (error == 1)
	{
		error = mp1.run();
		if (error==1)
		{	
			cout << " Calculation failed twice! Exiting." << endl;
			exit(-1);
		}
	}
  energy = mp1.getE(wstate) * 627.5;
  error = mp1.getGrad(grada[0],wstate);
  for (int i=0;i<N3;i++)
    grada[0][i] *= ANGtoBOHR;
  if (wstate2==0)
  for (int i=0;i<N3;i++)
    grad[i] = grada[0][i];

 //average energy and gradient
  if (wstate2>0)
  {
    //error = mp1.run(wstate2,0); //grad and energy
    energy += mp1.getE(wstate2) * 627.5;
    energy /= 2.; //average energy
    error = mp1.getGrad(grada[1],wstate2);
    for (int i=0;i<N3;i++)
      grada[1][i] *= ANGtoBOHR;
    if (wstate3==0)
    for (int i=0;i<N3;i++) 
      grad[i] = (grada[0][i] + grada[1][i])/2.;
  }
  if (wstate3>0)
  {
    //error = mp1.run(wstate3,0); //grad and energy
    energy = 2*energy + mp1.getE(wstate3) * 627.5;
    energy /= 3.; //average energy
    error = mp1.getGrad(grada[2],wstate3);
    for (int i=0;i<N3;i++)
      grada[2][i] *= ANGtoBOHR;
    for (int i=0;i<N3;i++) 
      grad[i] = (grada[0][i] + grada[1][i] + grada[2][i])/3.;
  }
  for (int i=0;i<nstates;i++)
    E[i] = mp1.getE(i+1) * 627.5;
	

#elif USE_TC
  tc1.reset(coords);
  if (gradcalls==0 && seedType<1)
  {
    //copy wfn 
    int runendCopy = runend;
    if (seedType==0) 
    {
      printf(" ERROR: seedType not set! \n");
      exit(1);
    }
    if (seedType==-1) runendCopy--;
    if (seedType==-2) runendCopy++;

		if (seedType==-1 || seedType==-2)
		{
     string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
   	 printf(" copying ORBKEY from ORBFILE%s to ORBFILE%s \n",runNameCopy.c_str(),runName0.c_str());
     tc1.readOrb=true;
   	 string cmd = "cp scratch/ORBFILE"+runNameCopy+" scratch/ORBFILE"+runName0;
   	 system(cmd.c_str());
		}
  }
  tc1.grads(runend,runNum); //grad and energy
  int error = tc1.get_grad(wstate,runend,runNum,grada[0]);
  if (error)
    printf(" Error = %i\n",error);

  energy = tc1.getE(wstate)*627.5; 
  for (int i=0;i<N3;i++)
    grada[0][i] *=ANGtoBOHR;
  if (wstate2==0)
  for (int i=0;i<N3;i++)
    grad[i] = grada[0][i];

  if (wstate2>0)
  {
	  //tc1.grads(wstate2,runend,runNum);
    error = tc1.get_grad(wstate2,runend,runNum,grada[1]);
    energy += tc1.getE(wstate2)*627.5;
    energy /= 2.; //average energy
    for (int i=0;i<N3;i++)
      grada[1][i] *= ANGtoBOHR;
    if (wstate3==0)
    for (int i=0;i<N3;i++) 
      grad[i] = (grada[0][i] + grada[1][i])/2.;
  }
  if (wstate3>0)
  {
	  //tc1.grads(wstate3,runend,runNum);
    error=tc1.get_grad(wstate3,runend,runNum,grada[2]);
    energy = 2*energy + tc1.getE(wstate3)*627.5;
    energy /= 3.; //average energy
    for (int i=0;i<N3;i++)
      grada[2][i] *= ANGtoBOHR;
    for (int i=0;i<N3;i++) 
      grad[i] = (grada[0][i] + grada[1][i] + grada[2][i])/3.;
  }
  for (int i=0;i<nstates;i++)
    E[i] = tc1.getE(i+1)*627.5;
#elif USE_ORCA
  energy = orca1.grads(coords,grad);
#elif USE_GAUSSIAN
  energy = gaus1.grads(coords,grad);
#elif USE_ASE
  energy = ase1.grads(coords,grad);
#else
  char* pbsPath;
  pbsPath = getenv ("PBSTMPDIR");
  string pdir = "";
  if (pbsPath!=NULL)
  {
    string pstr(pbsPath);
    pdir = pstr + "/";
  }
  Mopac mop1; 
  mop1.alloc(natoms);
  mop1.reset(natoms,anumbers,anames,coords);
  mop1.set_charge(CHARGE);
  string cmd = "mkdir "+pdir+"scratch";
  system(cmd.c_str());
  energy = mop1.grads(pdir+"scratch/mxyzfile"+runends);
  if (nforce>0)
    energy -= 100000.;
  for (int i=0;i<N3;i++)
    grad[i] = mop1.grad[i];
  mop1.freemem();
#endif
  gradcalls++;

  int success = 1;
  if (V0==0.0)  V0 = energy;
  if (energy-V0>1000. || energy-V0 < -1000.)  success = 0;

//  printf(" write_on: %i success: %i energy: %6.5f \n",write_on,success,energy);

#if WRITE_FILES
  if (write_on && success)
  {
    string nstr = StringTools::int2str(gradcalls+res_t,4,"0");
    string filename = "scratch/qcsave"+runName0+"_"+nstr;
    write_xyz_grad(coords,grad,filename);
  }
#endif

  return success;
}


void Gradient::write_xyz_grad(double* coords, double* grad, string filename)
{
  wrote_grad = 1;

  ofstream xyzfile;
  string xyzfile_string = filename+".xyz";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(15);

  xyzfile << natoms << endl << energy/627.5 << endl;
  for (int i=0;i<natoms;i++)
    xyzfile << anames[i] << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2] << endl;

  xyzfile.close();

  ofstream gradfile;
  string gradfile_string = filename+".grad";
  gradfile.open(gradfile_string.c_str());
  gradfile.setf(ios::fixed);
  gradfile.setf(ios::left);
  gradfile << setprecision(15);

  gradfile << natoms << endl << energy/627.5 << endl;
  for (int i=0;i<natoms;i++)
    gradfile << anames[i] << " " << grad[3*i+0] << " " << grad[3*i+1] << " " << grad[3*i+2] << endl;

  gradfile.close();

  return;
}

int Gradient::read_nstates()
{
  string filename = "NSTATES";
  wstate = 1; //default 
  wstate2 = 0;
  wstate3 = 0;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf(" couldn't find NSTATES file \n");
    exit(1);
  }

  string line;
  int nstates0 = 0;
  bool success =true; 
  int nf = 0;
	
  while (!infile.eof())
  {
		success=getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ");
    if (nf==0)
      nstates0 = atoi(tok_line[1].c_str());
    else if (nf==1)
    {
      wstate = atoi(tok_line[1].c_str());
      if (tok_line.size()>2) wstate2 = atoi(tok_line[2].c_str());
      if (tok_line.size()>3) wstate3 = atoi(tok_line[3].c_str());
      printf("  wstate: %i %i %i \n",wstate,wstate2,wstate3);
    }
    nf++;
	}

  infile.close();

  if (nstates0>4)
  {
    printf(" nstates too high (>4). exiting \n");
    exit(1);
  }
  if (nstates0<0)
  {
    printf(" nstates must be greater than 0. exiting \n");
    exit(1);
  }
  return nstates0;
}

int Gradient::read_molpro_init(string* &hf_lines)
{
  printf("   in read_molpro_init \n");
  if (seedType<1) return 0;
  if (seedType!=1 && seedType!=2) 
  {
    printf(" ERROR: seedType in read_molpro_init must be 0, 1 or 2 \n");
    exit(1);
  }

  string filename = "MOLPRO_INIT1";
  if (seedType==2)
    filename = "MOLPRO_INIT2";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    hf_lines = new string[1];
    printf("    couldn't find MOLPRO file %s. Proceeding without it. \n",filename.c_str());
    return 0;
  }
  printf("    reading file %s \n",filename.c_str());

  int nhf = 0;
  string* hf_lines0 = new string[100];
  string line;
  bool success=true;
  while (!infile.eof())
  {
    success=getline(infile, line);
    if (success)
      hf_lines0[nhf++] = line;
  }
  infile.close();

 //copy it over to save
  hf_lines = new string[nhf];
  for (int i=0;i<nhf;i++)
    hf_lines[i] = hf_lines0[i];
  delete [] hf_lines0;

  
  return nhf;
}

void Gradient::read_molpro_settings(int& nstates0, int& nclosed, int& nocc, int& nelec, string& basis)
{
  wstate = 1; //default 
  wstate2 = 0;
  wstate3 = 0;

  string filename = "MOLPRO";
  //printf(" reading file: %s \n",filename.c_str());


  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf(" Error: couldn't open settings file: %s \n",filename.c_str());
    exit(-1);
  }
  string line;
  bool success=true;
  int nf = 0;
  while (!infile.eof())
  {
    success=getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ");
    //cout << "RR0: " << line << endl; fflush(stdout);

    if (nf==0)
      nstates0 = atoi(tok_line[1].c_str());
    else if (nf==1)
    {
      wstate = atoi(tok_line[1].c_str());
      if (tok_line.size()>2) wstate2 = atoi(tok_line[2].c_str());
      if (tok_line.size()>3) wstate3 = atoi(tok_line[3].c_str());
      printf("  wstate: %i %i %i \n",wstate,wstate2,wstate3);
    }
    else if (nf==2)
      nclosed = atoi(tok_line[1].c_str());
    else if (nf==3)
      nocc = atoi(tok_line[1].c_str());
    else if (nf==4)
      nelec = atoi(tok_line[1].c_str());
    else if (nf==5)
      basis = tok_line[1];
    nf++;
  }
  infile.close();

  printf("   settings. nstates/wstate/nclosed/nocc/nelec/basis: %2i %2i %2i %2i %2i %s \n",nstates0,wstate,nclosed,nocc,nelec,basis.c_str());
  
  if (wstate2>0 && nstates0<2)
  {
    printf("  wstate2: %i nstates: %i. exiting \n",wstate2,nstates0);
    exit(1);
  }
  if (wstate3>0 && nstates0<3)
  {
    printf("  wstate3: %i nstates: %i. exiting \n",wstate3,nstates0);
    exit(1);
  }

	mp1.wstate=wstate; 
	mp1.wstate2=wstate2; 
	mp1.wstate3=wstate3;

  return;
}

void Gradient::read_tc_settings(int& nstates0, int& nclosed, int& nactive,  string& basis,string& method)
{
  wstate = 1; //default 
  wstate2 = 0;
  wstate3 = 0;

  string filename = "TeraChem";
  printf(" reading file: %s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf(" Error: couldn't open settings file: %s \n",filename.c_str());
    exit(-1);
  }
  string line;
  bool success=true;
  int nf = 0;
  while (!infile.eof())
  {
    success=getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ");
    if (nf==0)
      method = tok_line[1].c_str();
    else if (nf==1)
      nstates0 = atoi(tok_line[1].c_str());
    else if (nf==2)
    {
      wstate = atoi(tok_line[1].c_str());
      if (tok_line.size()>2) wstate2 = atoi(tok_line[2].c_str());
      if (tok_line.size()>3) wstate3 = atoi(tok_line[3].c_str());
      printf("  wstate: %i %i %i \n",wstate,wstate2,wstate3);
    }
    else if (nf==3)
      nclosed = atoi(tok_line[1].c_str());
    else if (nf==4)
      nactive = atoi(tok_line[1].c_str());
    else if (nf==5)
      basis = tok_line[1];
    nf++;
  }
  infile.close();
  if (method == "FOMO")
    printf("   settings. nstates/wstate/nclosed/nactive/basis: %2i %2i %2i %2i %s \n",nstates0,wstate,nclosed,nactive,basis.c_str());
  else if (method == "DFT")
    printf("   settings. basis: %s \n",basis.c_str());
  
  if (wstate2>0 && nstates0<2)
  {
    printf("  wstate2: %i nstates: %i. exiting \n",wstate2,nstates0);
    exit(1);
  }
  if (wstate3>0 && nstates0<3)
  {
    printf("  wstate3: %i nstates: %i. exiting \n",wstate3,nstates0);
    exit(1);
  }

	tc1.wstate=wstate; 
	tc1.wstate2=wstate2; 
	tc1.wstate3=wstate3;

  return;
}

void Gradient::init(string infilename, int natoms0, int* anumbers0, string* anames0, double* coords0, int run, int rune, int ncpu, int knnr_level, int q1)
{

#if QCHEM && QCHEMSF
  printf(" Cannot use both QCHEM and QCHEMSF \n");
  exit(-1);
#endif
#if QCHEM && USE_ASE
  printf(" Cannot use both QCHEM and ASE \n");
  exit(-1);
#endif
#if QCHEM && USE_GAUSSIAN
  printf(" Cannot use both QCHEM and GAUSSIAN \n");
  exit(-1);
#endif
#if QCHEM && USE_ORCA
  printf(" Cannot use both QCHEM and ORCA \n");
  exit(-1);
#endif
#if USE_GAUSSIAN && USE_ASE
  printf(" Cannot use both ASE and GAUSSIAN \n");
  exit(-1);
#endif

  V0 = 0.; 
  fdE = 0.;

  xyz_grad = 0;
  always_do_exact = 0;
  write_on = 1;
  wrote_grad = 0;
  res_t = 0;
  gradcalls = 0;
  nscffail = 0;
  nstates = 1;
  CHARGE = q1;

  natoms = natoms0;
  N3 = natoms*3;
  anumbers = new int[natoms+1];
  anames = new string[natoms+1];

  for (int i=0;i<natoms0;i++)
    anumbers[i] = anumbers0[i];
  for (int i=0;i<natoms0;i++)
    anames[i] = anames0[i];

  runNum = run;
  runend = rune;
  string nstr = StringTools::int2str(run,4,"0");
  runName0 = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runend,4,"0");
  runends = nstr;

#if QCHEM
  qchem1.init(infilename,natoms,anumbers,anames,run,rune);
  qchem1.ncpu = ncpu;
#endif
#if QCHEMSF
  nstates = read_nstates();
  qchemsf1.init(infilename,natoms,anumbers,anames,run,rune);
  qchemsf1.ncpu = ncpu;
#endif
#if USE_MOLPRO
  int nstates0;
  int nclosed,nocc,nelec;
  string basis;
  read_molpro_settings(nstates0,nclosed,nocc,nelec,basis);
  string* hf_lines;
  int nhf_lines = 0;
  if (seedType==3)
    printf("   assuming wfn already seeded in scratch \n");
  else if (seedType>0)
    read_molpro_init(hf_lines);
  nstates = nstates0;
  mp1.init(nstates,nclosed,nocc,nelec,natoms,anames,coords0,ncpu,basis,run);
  mp1.runname("mp_"+runName0);
  printf("   runname set to mp_%s \n",runName0.c_str());
  if (seedType>0 && seedType!=3)
  {
    mp1.init_hf(nhf_lines,hf_lines);
    mp1.seed();
    delete [] hf_lines;
  }
#endif
#if USE_TC
  int nstates0;
  int nclosed,nactive;
  string basis;
  string method;
  read_tc_settings(nstates0,nclosed,nactive,basis,method);
  nstates = nstates0;
  tc1.init(natoms,nstates,nclosed,nactive,anames,coords0,basis,method,run,rune);
  if (seedType==3)
  {
    printf("   assuming wfn already seeded in scratch \n");
    tc1.readOrb = true;
  }
  if (seedType>0 && seedType!=3)
    tc1.readOrb = false;
#endif
#if USE_ORCA
  orca1.init(infilename,natoms,anumbers,anames,run,rune);
  orca1.ncpu = ncpu;
#endif
#if USE_GAUSSIAN
  gaus1.init(infilename,natoms,anumbers,anames,run,rune);
  gaus1.CHARGE = CHARGE;
  gaus1.ncpu = ncpu;
#endif
#if USE_ASE
  ase1.init(infilename,natoms,anumbers,anames,run,rune);
  ase1.CHARGE = CHARGE;
  ase1.ncpu = ncpu;
#endif

  E = new double[nstates];
  for (int i=0;i<nstates;i++) E[i] = 0.;
#if QCHEMSF || USE_MOLPRO || USE_TC
  grada = new double*[nstates];
  for (int i=0;i<nstates;i++)
    grada[i] = new double[N3];
#endif

  dE = new double[nstates-1];
  for (int i=0;i<nstates-1;i++) dE[i]=1000.;

  knnr_active = 0;
  knn_k = KNN_K;
#if USE_KNNR
  if (knnr_level)
  {
    knnr1.printl = 0;
    int kfound = 0;
#if 1
    knnr_inited = 0;
#else
    printf("  loading KNNR files \n");
    kfound = knnr1.begin(runNum,natoms);
#endif
    knnr_active = knnr_level;
   // knnr1.test_points();
    if (kfound>0)
      res_t = 1000;
  }
#endif

  string ffile = "FORCE"+nstr;
  nforce = force_init(ffile);

#if QCHEM
  printf("  grad initiated: Q-Chem mode \n");
#elif QCHEMSF
  printf("  grad initiated: Q-Chem Spin-Flip mode \n");
#elif USE_ORCA
  printf("  grad initiated: ORCA mode \n");
#elif USE_GAUSSIAN
  printf("  grad initiated: Gaussian mode \n");
#elif USE_ASE
  printf("  grad initiated: ASE mode \n");
#elif USE_MOLPRO
  printf("  grad initiated: MOLPRO mode \n");
#elif USE_TC
  printf("  grad initiated: TeraChem mode \n");
#else
  printf("  grad initiated: Mopac mode \n");
#endif

  //printf(" grad init knnr_active: %i \n",knnr_active);

  return;
}

int Gradient::force_init(string ffile)
{
  int use_force = 0;

  int* fa1 = new int[20];
  double* fv1 = new double[10];
  double* fk1 = new double[10];

  ifstream infile;
  infile.open(ffile.c_str());
  if (!infile)
  {
    string ffile1 = "scratch/"+ffile;
    infile.open(ffile1.c_str());
    if (!infile)
      return 0;
  }


  string line;
  bool success=true;
  int nf = 0;
  while (!infile.eof())
  {
    success=getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ");
    if (tok_line.size()>3)
    {
      //cout << "RR: " << line << endl;
      fa1[2*use_force+0] = atoi(tok_line[0].c_str()) -1;
      fa1[2*use_force+1] = atoi(tok_line[1].c_str()) -1;
      fv1[use_force] = atof(tok_line[2].c_str());
      fk1[use_force] = atof(tok_line[3].c_str());
      use_force++;
    }
  }
  //printf("  found %i forces \n",use_force); fflush(stdout);

 //save forces
  fa = new int[2*use_force];
  fv = new double[use_force];
  fk = new double[use_force];

  for (int i=0;i<2*use_force;i++)
    fa[i] = fa1[i];
  for (int i=0;i<use_force;i++)
    fv[i] = fv1[i];
  for (int i=0;i<use_force;i++)
    fk[i] = fk1[i];

  delete [] fa1;
  delete [] fv1;
  delete [] fk1;

#if 1
  printf("    forces:");
  for (int i=0;i<use_force;i++)
    printf("   %i-%i  d: %8.5f k: %8.5f ",fa[2*i+0]+1,fa[2*i+1]+1,fv[i],fk[i]);
  printf("\n");
#endif

  return use_force;
}

int Gradient::dvec_calc(double* coords, double* dvec,int run,int rune)
{ 
  runNum = run;
  runend = rune;
  string nstr = StringTools::int2str(run,4,"0");
  runName0 = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runend,4,"0");
	int error=1;
#if USE_MOLPRO
  mp1.reset(coords);
	mp1.runname("mp_"+runName0);
  mp1.calc_dvec();
  error = mp1.getDVec(dvec);
	if (error==1)
	{
    int runendCopy = runend-1;
    string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
    printf(" copying wfn from mp_%s to mp_%s \n",runNameCopy.c_str(),runName0.c_str());
    string cmd = "cp scratch/mp_"+runNameCopy+" scratch/mp_"+runName0;
    system(cmd.c_str());
    //if using text files
    //cmd = "cp mp_"+runNameCopy+".orb mp_"+runName0+".orb";
    //printf(" copying orbs from mp_%s.orb to mp_%s.orb \n",runNameCopy.c_str(),runName0.c_str());
    //system(cmd.c_str());
		mp1.calc_dvec();
  	error = mp1.getDVec(dvec);
		if (error==1)
			{	
				printf(" Dvec Calculation failed twice!");
				runendCopy--;
    		string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
				if (runendCopy<0)
					exit(-1);
    		//printf(" copying wfn from mp_%s to mp_%s \n",runNameCopy.c_str(),runName0.c_str());
    		string cmd = "cp scratch/mp_"+runNameCopy+" scratch/mp_"+runName0;
    		system(cmd.c_str());
				mp1.calc_dvec();
  			error = mp1.getDVec(dvec);
				if (error==1)
				{
					printf(" Exiting.");
					exit(-1);
				}
			}
	}
	for (int i=0;i<N3;i++)
    dvec[i] *= ANGtoBOHR;
#elif USE_TC
  tc1.reset(coords);
  //assuming dvec already calcd!
  //tc1.calc_dvec(rune,run);
  tc1.get_dvec(rune,run,dvec);
	for (int i=0;i<N3;i++)
    dvec[i] *= ANGtoBOHR;
#elif QCHEMSF
	qchemsf1.calc_dvec(coords);
	qchemsf1.get_dvec(dvec);
#else
	printf(" not yet implemented \n");
  exit(-1);
#endif
	
#if 0
  printf(" dvec (1/Ang): \n");
  for (int i=0;i<natoms;i++)
    printf(" %s %12.10f %12.10f %12.10f \n",anames[i].c_str(),dvec[3*i+0],dvec[3*i+1],dvec[3*i+2]);
#endif

	return error;
}

double Gradient::energy_initial(double* coords,int run, int rune,int penalty, double sigma)
{
  printf(" Calculating initial energy\n");
  runNum = run;
  runend = rune;
  string nstr = StringTools::int2str(run,4,"0");
  runName0 = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runend,4,"0");
  int error=1;
  //cout << " seedType is "<< seedType << endl;
#if USE_MOLPRO
  mp1.reset(coords);
  mp1.runname("mp_"+runName0);
  error = mp1.calc_energy();
  if (error == 1)
  {
    error = mp1.calc_energy();
    if (error==1)
      {
        cout << " Calculation failed twice! Exiting." << endl;
        exit(-1);
      }
  }
#else
  printf(" not yet implemented \n");
  return error;
#endif

  if (penalty)
  {
    if (wstate2<2)
    {
      printf(" Need to have more than 2 states\n");
      exit(-1);
    }
    energy = mp1.getE(wstate) * 627.5;
    energy += mp1.getE(wstate2) * 627.5;
    energy /= 2.0;
    printf(" sigma is %1.2f\n",sigma);
    //double sigma = 1.0; ///BUUUUUG
    double alpha = 0.02*627.5; //kcal/mol
    double deltaE = mp1.getE(wstate2)*627.5 - mp1.getE(wstate)*627.5;
    double G = pow(deltaE,2)/(deltaE+alpha);
    energy += sigma*G;
  }
  else
  {
    if (wstate2)
      energy = mp1.getE(wstate2) * 627.5;
    else
      energy = mp1.getE(wstate) * 627.5;
  }

  for (int i=0;i<nstates;i++)
    E[i] = mp1.getE(i+1) * 627.5;

  //for (int i=0;i<nstates;i++)
  //  printf(" E[%i]=%1.4f\n",i,E[i]);

  return energy;
}

double Gradient::levine_penalty(double* coords, double* grad, double* Ut, int type,double sigma,double dmw)
{
  wrote_grad = 0;
  energy = -99.;
	double alpha = 0.02*627.5; //kcal/mol
  double deltaE=0.0;
	//printf(" In Levine penalty fxn\n");
  int do_exact = 1;
  if (do_exact)
  {
    printf(" eg"); fflush(stdout);
		if (wstate2<2)
		{
			printf(" Need to have more than 2 states\n");
			exit(-1);
		}		

#if USE_MOLPRO
  	mp1.reset(coords);
  	if (gradcalls==0 && seedType<1)
  	{
  	  //copy wfn 
  	  int runendCopy = runend;
  	  if (seedType==0) 
  	  {
  	    printf(" ERROR: seedType not set! \n");
  	    exit(1);
  	  }
  	  if (seedType==-1) runendCopy--;
  	  if (seedType==-2) runendCopy++;

			if (seedType==-1 || seedType==-2)
			{
  	 	 string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
  	 	 printf(" copying wfn from mp_%s to mp_%s \n",runNameCopy.c_str(),runName0.c_str());
  	 	 string cmd = "cp scratch/mp_"+runNameCopy+" scratch/mp_"+runName0;
  	 	 system(cmd.c_str());
			}
  	}
  	int error = mp1.run(); //grad and energy
		if (error == 1)
		{
			error = mp1.run();
			if (error==1)
			{	
				cout << " Calculation failed twice! Exiting." << endl;
				exit(-1);
			}
		}
  	error = mp1.getGrad(grada[0],wstate);
		if (error==1)
		{
			printf(" error in get grad\n");
			exit(-1);
		}
  	error = mp1.getGrad(grada[1],wstate2);
		if (error==1)
		{
			printf(" error in get grad\n");
			exit(-1);
		}
  	for (int i=0;i<N3;i++)
  	  grada[0][i] *= ANGtoBOHR;
  	for (int i=0;i<N3;i++)
  	  grada[1][i] *= ANGtoBOHR;
  	for (int i=0;i<nstates;i++)
  	  E[i] = mp1.getE(i+1) * 627.5;
#elif QCHEMSF
  	qchemsf1.calc_grads(coords);
  	qchemsf1.getGrad(wstate,grada[0]);
    qchemsf1.getGrad(wstate2,grada[1]);
  	for (int i=0;i<nstates;i++)
  	  E[i] = qchemsf1.getE(i);
#elif USE_TC
  tc1.reset(coords);
  tc1.docoupling=false;
  if (gradcalls==0 && seedType<1)
  {
    //copy wfn 
    int runendCopy = runend;
    if (seedType==0) 
    {
      printf(" ERROR: seedType not set! \n");
      exit(1);
    }
    if (seedType==-1) runendCopy--;
    if (seedType==-2) runendCopy++;

		if (seedType==-1 || seedType==-2)
		{
   	 string runNameCopy = StringTools::int2str(runNum,4,"0")+"_"+StringTools::int2str(runendCopy,4,"0");
   	 printf(" copying ORBKEY from ORBFILE%s to ORBFILE%s \n",runNameCopy.c_str(),runName0.c_str());
     tc1.readOrb=true;
   	 string cmd = "cp scratch/ORBFILE"+runNameCopy+" scratch/ORBFILE"+runName0;
   	 system(cmd.c_str());
		}
  }
  tc1.grads(runend,runNum); //grad and energy
	//tc1.grads(wstate2,runend,runNum);
  int error = tc1.get_grad(wstate,runend,runNum,grada[0]);
  error = tc1.get_grad(wstate2,runend,runNum,grada[1]);
  for (int i=0;i<N3;i++)
    grada[0][i] *= ANGtoBOHR;
  for (int i=0;i<N3;i++)
    grada[1][i] *= ANGtoBOHR;
  for (int i=0;i<nstates;i++)
    E[i] = tc1.getE(i+1)*627.5;
#else
   printf(" NOT IMPLEMENTED\n");
   exit(-1);
#endif

		if (type==1)
		{
#if USE_MOLPRO
  		energy = mp1.getE(wstate) * 627.5;
  		energy+= mp1.getE(wstate2) * 627.5; 
			energy /= 2.0;
			deltaE = mp1.getE(wstate2)*627.5 - mp1.getE(wstate)*627.5;
#elif QCHEMSF
			energy = qchemsf1.getE(wstate-1);
			energy+= qchemsf1.getE(wstate2-1);
			energy /= 2.0;
			deltaE = qchemsf1.getE(wstate2-1)- qchemsf1.getE(wstate-1);
#elif USE_TC
      energy=(tc1.getE(wstate)*627.5+tc1.getE(wstate2)*627.5)/2.0;
      deltaE = tc1.getE(wstate2)*627.5 - tc1.getE(wstate)*627.5;
#endif
			//penalty energy and gradient
			double G = pow(deltaE,2)/(deltaE+alpha);	
			energy += sigma*G;                        	
  		for (int i=0;i<N3;i++) 
  			grad[i] = (grada[0][i] + grada[1][i])/2.;
			double factor = sigma*(pow(deltaE,2) + 2*alpha*deltaE)/pow(deltaE+alpha,2);
  		for (int i=0;i<N3;i++) 
  		  grad[i] += factor*(grada[1][i] - grada[0][i]);
		}
  	  xyz_grad = 1;
  }

#if 0
  if (xyz_grad)
  {
    printf( "Energy: %1.4f\n",energy);
    printf(" dE: %1.4f\n",deltaE);
    printf(" Grad: \n");
    for (int i=0;i<natoms;i++)
      printf(" %s %12.10f %12.10f %12.10f \n",anames[i].c_str(),grad[3*i+0],grad[3*i+1],grad[3*i+2]);
  }
#endif
  gradcalls++;
  int success = 1;
  if (V0==0.0)  V0 = energy;
  if (energy-V0>2000. || energy-V0 < -2000.)  success = 0;
	if (success == 0)
	{
		printf(" Energy is too large.\n");
		//printf(" Energy is %1.4f\n",energy);
		//printf(" V0= %1.4f\n",V0);
		//exit(-1);
	}

  return energy;
}
