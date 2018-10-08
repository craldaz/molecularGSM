#include "base.h"



void Base::initialize_icoords()
{
  icoords = new ICoord[nnmax+1];
  for (int i=0;i<nnmax;i++)
    icoords[i].alloc(natoms);
  for (int i=0;i<nnmax;i++)
    icoords[i].gradrms = 1.;
  for (int i=1;i<nnmax-1;i++)
    icoords[i].reset(natoms,anames,anumbers,coords[0]);
  icoords[0].reset(natoms,anames,anumbers,coords[0]);
  icoords[nnmax-1].reset(natoms,anames,anumbers,coords[nnmax-1]);
}

void Base::parameter_init(string infilename)
{

//  printf(" Hello 2\n");

  nnmax = 9;
  hessSSM = 0; //default no Hessian given
  tstype = 0;
  prodelim = 1000.;
  lastOpt = 0; //default do not optimize last SSM node 
  initialOpt = 0; //default do not optimize first node 
  DQMAG_SSM_MAX = 0.8;  
  DQMAG_SSM_MIN = 0.2;  
  QDISTMAX = 5.0;  
  PEAK4_EDIFF = 2.0;
  isRestart = 0; 
  isSSM = 0;
  isFSM = 0;
	isMECI =0;
	isPRODUCT =0;
	isMAP_DE =0;
	isMAP_SE=0;
	isSE_ESSM =0; //these are for excited state methods
  use_exact_climb = 2;
	restart_wfn=0;  //restart molpro calc with wfu

  printf("Initializing Tolerances and Parameters... \n");
  printf("  -Opening %s \n",infilename.c_str());

  ifstream infile;
  infile.open(infilename.c_str());
  if (!infile){
    printf("\n Error opening %s \n",infilename.c_str());
    exit(-1);
  }

  printf("  -reading file... \n");

  // pass infile to stringtools and get the line containing tag
  string tag="String Info";
  bool found=StringTools::findstr(infile, tag);
  if (!found) { cout << "Could not find tag for Default Info" << endl; exit(-1);}
  string line, templine, tagname;

  // parse the input section
  bool stillreading=true;
  while (stillreading)
  {
    stillreading=false;
    // set to false and set back to true if we read something
    // get filename
    getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ,\t");
    templine=StringTools::newCleanString(tok_line[0]);
    tagname=StringTools::trimRight(templine);
    // these variables are denoted by strings with same name
    if (tagname=="MAX_OPT_ITERS") { 
      MAX_OPT_ITERS=atoi(tok_line[1].c_str());
      stillreading=true;
      cout <<"  -MAX_OPT_ITERS: " << MAX_OPT_ITERS << endl;
    }
    if (tagname=="STEP_OPT_ITERS") {
      STEP_OPT_ITERS = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -STEP_OPT_ITERS: " << STEP_OPT_ITERS << endl;
    }
    if (tagname=="RESTART") {
      isRestart = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -RESTART: " << isRestart << endl;
    }
    if (tagname=="TS_FINAL_TYPE") {
      tstype = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -TS_FINAL_TYPE: " << tstype << endl;
      if (tstype!=0 && tstype!=1 && tstype!=2)
      {
        printf("  TS_FINAL_TYPE must be 0, 1 or 2 \n");
        exit(1);
      }
    }
    if (tagname=="PRODUCT_LIMIT") {
      prodelim = atof(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -PRODUCT_LIMIT: " << prodelim << endl;
    }
    if (tagname=="FINAL_OPT") {
      lastOpt = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -FINAL_OPT: " << lastOpt << endl;
    }
    if (tagname=="INITIAL_OPT") {
      initialOpt = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -INITIAL_OPT: " << initialOpt << endl;
    }
    if (tagname=="SSM_DQMAX") {
      DQMAG_SSM_MAX = atof(tok_line[1].c_str());
      DQMAG_SSM_MIN = DQMAG_SSM_MAX/4.;
      stillreading = true;
      cout <<"  -SSM_DQMAX: " << DQMAG_SSM_MAX << endl;
      cout <<"  -SSM_DQMIN: " << DQMAG_SSM_MIN << endl;
    }
    if (tagname=="MIN_SPACING") {
      QDISTMAX = atof(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -SSM_MIN_SPACING: " << QDISTMAX << endl;
    }
    if (tagname=="INT_THRESH") {
      PEAK4_EDIFF = atof(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -INT_THRESH: " << PEAK4_EDIFF << endl;
    }
	  if (tagname=="RESTART_WFN"){
			restart_wfn = atoi(tok_line[1].c_str());
			stillreading=true;
			cout <<"  -restart_wfn = " << restart_wfn << endl; 
		}
    if (tagname=="SM_TYPE") {
      if (tok_line[1]=="SSM")
      {
        printf("  -using SSM \n");
        isFSM = 0;
        isSSM = 1;
      }
      else if (tok_line[1]=="SSMFR")
      {
        printf("  -using SSM with initial Hessian \n");
        isFSM = 0;
        isSSM = 1;
        hessSSM = 1;
      }
      else if (tok_line[1]=="FSM")
      {
        printf("  -using FSM \n");
        isFSM = 1;
        isSSM = 0;
      }
      else if (tok_line[1]=="OPT")
      {
        printf("  -using OPT \n");
        isFSM = -1;
        isSSM = -1;
      }
			else if (tok_line[1]=="MECI")
			{
				printf("  -using MECI \n");
				isMECI = 1;
			}
			else if (tok_line[1]=="PRODUCT")
			{
				printf("  -using PRODUCT\n");
				isPRODUCT = 1;
			}
			else if (tok_line[1]=="SE-ESSM")
			{
				printf("  -using SE-ESSM \n");
				isSSM = 1; //both are turned on 
				isSE_ESSM=1;
			}
			else if (tok_line[1]=="DE-MAP")
			{
				printf("  -using DE-MAP \n");
				isMAP_DE=1;
			}
			else if (tok_line[1]=="SE-MAP")
			{
				printf("  -using SE-MAP \n");
				isMAP_SE=1;
			}
      else
      {
        printf("  -using GSM \n");
        isFSM = 0;
        isSSM = 0;
      }
      stillreading=true;
    }
    if (tagname=="CONV_TOL") {
      CONV_TOL=atof(tok_line[1].c_str());
      stillreading=true;
      cout <<"  -CONV_TOL = " << CONV_TOL << endl;
    }
    if (tagname=="ADD_NODE_TOL"){
      ADD_NODE_TOL=atof(tok_line[1].c_str());
      stillreading=true;
      cout <<"  -ADD_NODE_TOL = " << ADD_NODE_TOL << endl;
    }
    if (tagname=="SCALING"){
      SCALING = atof(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -SCALING = " << SCALING << endl;
    }
    if (tagname=="BOND_FRAGMENTS"){
      bondfrags = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -BOND_FRAGMENTS = " << bondfrags << endl;
    }
    if (tagname=="GROWTH_DIRECTION"){
      GROWD = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -GROWTH_DIRECTION = " << GROWD << endl;
    }
    if (tagname=="CLIMB_TS"){
      use_exact_climb = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -CLIMB_TS = " << use_exact_climb << endl;
    }
    if (tagname=="nnodes" || tagname=="NNODES"){
      nnmax = atoi(tok_line[1].c_str());
      stillreading = true;
      cout <<"  -NNODES = " << nnmax << endl;
      if (nnmax < 3)
      {
        printf("\n\n ERROR: NNODES cannot be less than 3 \n");
        exit(1);
      }
    }

  } //while stillreading
  infile.close();
 
  if (tstype==2)
  {
    printf("  TS_FINAL_TYPE == 2, turning off climbing image and TS search \n");
    use_exact_climb = 0;
  }
  nnmax0 = nnmax;

  printf(" Done reading inpfileq \n\n");
}


void Base::structure_init(string xyzfile)
{
  printf("Reading and initializing string coordinates \n");
  printf("  -Opening structure file \n");

  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf("\n Error opening xyz file: %s \n",xyzfile.c_str());
    exit(-1);
  }

  printf("  -reading file... \n");

  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  cout <<"  -The number of atoms is: " << natoms << endl;

  success=getline(infile, line);
  vector<string> tok_line0 = StringTools::tokenize(line, " \t");
  CHARGE = 0;
  if (tok_line0.size()>0)
    CHARGE = atoi(tok_line0[0].c_str());
  if (CHARGE>5 || CHARGE<-5)
  {
    printf("   invalid charge value in initial.xyz: %2i \n",CHARGE);
    exit(1);
  }

  anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];
  frozen = new int[1+natoms];
  for (int i=0;i<natoms;i++)
    frozen[i] = 0;

  printf("  -Reading the atomic names...");
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    anumbers[i]=PTable::atom_number(anames[i]);
    amasses[i]=PTable::atom_mass(anumbers[i]);
  }

  infile.close();

  coords = new double*[1+nnmax];
  tangents = new double*[1+nnmax];
  grads = new double*[1+nnmax];
  perp_grads = new double*[1+nnmax];


  V_profile = new double[1+nnmax];
  for (int i=0;i<nnmax;i++)
    V_profile[i] = 0.;

  for (int i=0;i<nnmax;i++){
    coords[i] = new double[1+natoms*3];
    tangents[i] = new double[1+natoms*3];
    grads[i] = new double[1+natoms*3];
    perp_grads[i] = new double[1+natoms*3];
  }

  printf("  -Reading coordinates...");
  printf("Opening xyz file \n");
  infile.open(xyzfile.c_str());

  printf(" Here %i\n",isSSM);
  for (int i=0;i<2;i++)
  {
    if ((isSSM) && i==1) break;
    success=getline(infile, line);
    success=getline(infile, line);
    for (int j=0;j<natoms;j++)
    {
      if (infile.eof())
      {
        printf("   end of xyz file reached early, exiting \n");
        exit(1);
      }
      success=getline(infile, line);
      int length=StringTools::cleanstring(line);
      vector<string> tok_line = StringTools::tokenize(line, " \t");
//      cout << " i: " << i << " string: " << line << endl;
      int n;
      if (i==0) n = 0;
      else if (i==1) n = nnmax-1;
      coords[n][3*j+0]=atof(tok_line[1].c_str());
      coords[n][3*j+1]=atof(tok_line[2].c_str());
      coords[n][3*j+2]=atof(tok_line[3].c_str());
      perp_grads[i][3*j+0] = 0.0;
      perp_grads[i][3*j+1] = 0.0;
      perp_grads[i][3*j+2] = 0.0;
    }
  }

  if (isSSM)
  for (int i=0;i<3*natoms;i++)
    coords[nnmax-1][i] = coords[0][i];

  //cout << " done" << endl;
  infile.close();

  double zero = 0.;
  for (int i=0;i<3*natoms;i++)
    zero += coords[0][i]*coords[0][i];
  if (zero<0.0001) 
  {
    printf("\n ERROR: initial.xyz has NULL coordinates \n");
    exit(1);
  }

  printf("  printing frozen list:");
  for (int i=0;i<natoms;i++)
    printf(" %i",frozen[i]);
  printf("\n");

  printf("Finished reading information from structure file \n");
 
  return;
}

int Base::isomer_init(string isofilename)
{

  printf(" reading isomers \n");
  if (bondfrags == 1)
    printf("  WARNING: ignoring BONDS in ISOMERS file because BOND_FRAGMENTS == 1 \n");

  nfound = 0;

  nadd = 0;
  nbrk = 0;
  nangle = 0;
  ntors = 0;

  int maxab = 10;
  bond = new int[2*maxab];
  add = new int[2*maxab]; 
  brk = new int[2*maxab];
  angles = new int[3*maxab];
  anglet = new double[3*maxab];
  tors = new int[4*maxab];
  tort = new double[4*maxab];
  for (int i=0;i<2*maxab;i++) bond[i] = -1;
  for (int i=0;i<2*maxab;i++) add[i] = -1;
  for (int i=0;i<2*maxab;i++) brk[i] = -1;
  for (int i=0;i<3*maxab;i++) angles[i] = -1;
  for (int i=0;i<3*maxab;i++) anglet[i] = 999.;
  for (int i=0;i<4*maxab;i++) tors[i] = -1;
  for (int i=0;i<4*maxab;i++) tort[i] = 999.;

  ifstream output(isofilename.c_str(),ios::in);
  if (!output)
  {
    printf(" couldn't find ISOMERS file: %s \n",isofilename.c_str());
    return 0;
  }

  string line;
  vector<string> tok_line;
  while(!output.eof())
  {
    getline(output,line);
    //cout << " RR " << line << endl;
    if (line.find("BOND")!=string::npos && bondfrags==0)
    {
      tok_line = StringTools::tokenize(line, " \t");
      bond[2*nbond] = atoi(tok_line[1].c_str()) -1;
      bond[2*nbond+1] = atoi(tok_line[2].c_str()) -1;
      printf(" bond for coordinate system: %i %i \n",bond[2*nbond]+1,bond[2*nbond+1]+1);
      nbond++;
      if (nbond>maxab) break;
    }
    if (line.find("ADD")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      add[2*nadd] = atoi(tok_line[1].c_str()) -1;
      add[2*nadd+1] = atoi(tok_line[2].c_str()) -1;
      printf(" adding bond: %i %i \n",add[2*nadd]+1,add[2*nadd+1]+1);
      //if (!geoms[id].bond_exists(add[nfound][2*nadd],add[nfound][2*nadd+1]))
        nadd++;
      if (nadd>maxab) break;
    }
    if (line.find("BREAK")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      brk[2*nbrk] = atoi(tok_line[1].c_str()) -1;
      brk[2*nbrk+1] = atoi(tok_line[2].c_str()) -1;
      printf(" breaking bond: %i %i \n",brk[2*nbrk]+1,brk[2*nbrk+1]+1);
      //if (geoms[id].bond_exists(brk[nfound][2*nbrk],brk[nfound][2*nbrk+1]))
        nbrk++;
      if (nbrk>maxab) break;
    }
    if (line.find("ANGLE")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      angles[3*nangle+0] = atoi(tok_line[1].c_str()) -1;
      angles[3*nangle+1] = atoi(tok_line[2].c_str()) -1;
      angles[3*nangle+2] = atoi(tok_line[3].c_str()) -1;
      anglet[nangle] = atof(tok_line[4].c_str());
      printf(" angle: %i %i %i align to %4.3f \n",angles[3*nangle+0]+1,angles[3*nangle+1]+1,angles[3*nangle+2]+1,anglet[nangle]);
      nangle++;
      if (nangle>maxab) break;
    }
    if (line.find("TORSION")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      tors[4*ntors+0] = atoi(tok_line[1].c_str()) -1;
      tors[4*ntors+1] = atoi(tok_line[2].c_str()) -1;
      tors[4*ntors+2] = atoi(tok_line[3].c_str()) -1;
      tors[4*ntors+3] = atoi(tok_line[4].c_str()) -1;
      tort[ntors] = atof(tok_line[5].c_str());
      printf(" tor: %i %i %i %i align to %4.3f \n",tors[4*ntors+0]+1,tors[4*ntors+1]+1,tors[4*ntors+2]+1,tors[4*ntors+3]+1,tort[ntors]);
      ntors++;
      if (ntors>maxab) break;
    }
  }
  if (nadd > 0 || nbrk > 0 || nangle > 0 || ntors > 0) nfound++;

  printf(" found %i isomer",nfound);
  if (nfound!=1) printf("s");
  printf("\n\n");

#if DRIVE_ADD_TETRA
  //printf("\n now adding tetrahedral destinations as drivers \n");
  //Note: implemented in break_planes_ssm()
#endif


  return nfound;
}

/**
 * Adds a node between n1 and n3. In SE, n2 is nnmax-1
 * arbitrarily because the added node is not "between" any node. 
 */
int Base::addNode(int n1, int n2, int n3) //maybe need to pass ICoord objects to this function
{
  printf(" adding node: %i between %i %i \n",n2,n1,n3);

  if (n1==n2 || n2==n3 || n1==n3)
  {
    printf(" cannot add node, exiting \n"); 
    exit(1);
  }

  //TODO these need to be allocated here
  ICoord newic,intic;
  newic.alloc(natoms);
  intic.alloc(natoms);

#if USE_MOLPRO || USE_TC
  if (n2>n1) icoords[n2].grad1.seedType = -1;
  if (n2<n1) icoords[n2].grad1.seedType = -2;
#endif

  int nbonds = newic.nbonds;
  int nangles = newic.nangles;
  int ntor = newic.ntor;
  int size_ic = newic.nbonds+newic.nangles+newic.ntor;
  int len_d = newic.nicd0;
  double* ictan = new double[size_ic];
  double* ictan0 = new double[size_ic];

  double BDISTMIN = 0.01;
  double bdist = 0.;

	if (isSE_ESSM)
		BDISTMIN=0.05;

// Add a node
  int iR,iP,wR,wP,iN;
  for (int n=0;n<1;n++)
  {
    iR = n1;
    iP = n3;
    iN = n2;
    printf(" iR,iP: %i %i iN: %i ",iR,iP,iN);

    newic.reset(natoms,anames,anumbers,icoords[iR].coords);
    intic.reset(natoms,anames,anumbers,icoords[iP].coords);

    newic.update_ic();
    intic.update_ic();

    if (isSSM)
    {

      //TODO tangent_1b needs to be in base class
      bdist = tangent_1b(ictan);
			icoords[iN].bdist=bdist;
      printf(" bdist: %4.3f \n",bdist);
      if (bdist<BDISTMIN) break;
			if (bdist>icoords[iR].bdist) break; //not getting smaller
    }
    else
      tangent_1(ictan); //TODO same

#if 0
    printf(" printing ictan \n");
    for (int i=0;i<nbonds;i++)
      printf(" %1.2f",ictan[i]);
    printf("\n");
    for (int i=0;i<nangles;i++)

      printf(" %1.2f",ictan[nbonds+i]);
    printf("\n");
    if (ntors>0)
    for (int i=0;i<ntor;i++)
      printf(" %1.2f",ictan[nbonds+nangles+i]);
    printf("\n");
#endif

    double dqmag = 0.;

    newic.bmatp_create();
    newic.bmatp_to_U();
    newic.opt_constraint(ictan);

    if (isSSM)
    {
      dqmag = get_ssm_dqmag(bdist);
      if (tstype==2)
        dqmag = DQMAG_SSM_MAX;
    }
    else
    {
      for (int i=0;i<size_ic;i++) ictan0[i] = ictan[i];

      for (int j=0;j<size_ic-ntor;j++)
        dqmag += ictan0[j]*newic.Ut[newic.nicd*size_ic+j];
      for (int j=nbonds+nangles;j<size_ic;j++)
        dqmag += ictan0[j]*newic.Ut[newic.nicd*size_ic+j]; //CPMZ check
    }

    printf(" dqmag: %1.2f",dqmag);

//    newic.bmatp_create();
//    newic.bmatp_to_U();
//    newic.opt_constraint(ictan);
    newic.bmat_create();
    if (nnmax-nn!=1)
      newic.dq0[newic.nicd0-1] = -dqmag/(nnmax-nn);
    else
      newic.dq0[newic.nicd0-1] = -dqmag/2;
    if (isSSM)
      newic.dq0[newic.nicd0-1] = -dqmag; //CPMZ check

    printf(" dq0[constraint]: %1.2f \n",newic.dq0[newic.nicd0-1]);
    int success = newic.ic_to_xyz();

    if (isSSM && !success)
    {
      newic.dq0[newic.nicd0-1] = -dqmag/2.0;
      success = newic.ic_to_xyz();
      if (!success)
      {
        printf(" ERROR: couldn't add node, dqmag: %4.3f \n",dqmag);
        exit(1);
      }
      else
        printf(" add node working second time \n");
    }
    newic.update_ic();

    icoords[iN].reset(natoms,anames,anumbers,newic.coords);
    com_rotate_move(iR,iP,iN,1.0); //operates on iN via newic

    icoords[iN].bmatp_create();
    icoords[iN].bmatp_to_U();
    icoords[iN].bmat_create();

    if (!isSSM)
    {
      icoords[iN].make_Hint();
#if HESS_TANG
      //not yet implemented
      //get_eigenv_finite(iN);
#endif
      icoords[iN].newHess = 5;
    }
    else
    {
      for (int i=0;i<size_ic*size_ic;i++)
        icoords[iN].Hintp[i] = icoords[iR].Hintp[i];
      icoords[iN].newHess = 2;
    }

    active[iN] = 1;

  } //loop over interpolation

  delete [] ictan;
  delete [] ictan0;

#if CLOSE_DIST_ADD
  if (close_dist_fix(1))
    addNode(n1,n2,n3);
  else
#endif
  int success = 1;
  if (isSSM)
  {
    if (bdist>=BDISTMIN && icoords[n2].bdist < icoords[n1].bdist)
      nn++;
    else
      success = 0;
  }
  else 
    nn++;

  return success;
}