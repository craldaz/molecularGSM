#include "main.h"
void check_bonding(int nstruct, ICoord ic1, int natoms1, string* anames, int* anumbers, vector<double*> xyzalla, int* unique);
void procedure_1(string xyzfile, string targetfile);
void procedure_2(int nconf, string xyzfile, string targetfile);


#define DIST_THRESH 0.2


void opt_mopac(int charge, int natoms, string* anames, int* anumbers, vector<double*> xyzall, double* E, int type)
{
  int N3 = natoms*3;
  int nstruct = xyzall.size();

  int* frz = new int[natoms];
  for (int i=0;i<natoms;i++) frz[i] = 0;
  Mopac mopt;
  mopt.alloc(natoms);
  mopt.set_charge(charge);
  // AD - Feb 16, 2017
  // pointer called freeez which is an array of size natoms
  // what are we freezing, and why are we freezing it?
  mopt.freeze_d(frz);

  int* done = new int[nstruct];
  for (int i=0;i<nstruct;i++) done[i] = 0;

  for (int i=0;i<nstruct;i++)
  if (!done[i])
  {
    // AD - Feb 16, 2017
    // what is the difference between these three types
    // 1 = mopt = ligand only
    // 2 = moptb = ligand + metal
    // 3 = ? 
    string nstr = StringTools::int2str(i,3,"0");
    string mfilename;
    if (type==1)
      mfilename = "scratch/mopt"+nstr+".xyz";
    else if (type==2)
      mfilename = "scratch/moptb"+nstr+".xyz";
    else
      mfilename = "scratch/moptc"+nstr+".xyz";
  
//    ICoord ic1,ic2;
//    ic1.init(natoms,anumbers,anames,xyzall[i]);
//    ic1.ic_create();

    mopt.reset(natoms,anumbers,anames,xyzall[i]);
    E[i] = mopt.opt_check(mfilename);
    for (int j=0;j<N3;j++)
      xyzall[i][j] = mopt.xyz[j];

//    ic2.init(natoms,anumbers,anames,xyzall[i]);
//    ic2.ic_create();

//    int intact = compare_ic(ic1,ic2);

    done[i] = 1;
  }

#if 1
  printf("  conformer energies:");
  for (int i=0;i<nstruct;i++)
    printf(" %5.2f",E[i]);
  printf("\n");
#endif

  delete [] done;
  delete [] frz;

  return;
}

int generate_conformers_and_opt(int nconf, string filename, double* &E, vector<double*> &xyzall)
{
//  string xyzfile = filename;
  
  printf("  generating %4i conformers of %s \n",nconf,filename.c_str());
  string confile = "all.xyz";
  
  string nconfstr = StringTools::int2str(nconf,1,"0");
  string cmd = "obabel "+filename+" -O "+confile+" --confab --conf "+nconfstr;
  system(cmd.c_str());

  string cmd2 = "obabel " +filename+ " -oconfabreport -xf temgg.sdf -xr 1.0";
  system(cmd2.c_str());

  //Amanda Dewyer - Feb 20, 2017
  //appends original "test.xyz" file to the end of all.xyz so that it is sampled
  ifstream infile;
  fstream outfile;
  infile.open("test.xyz");
  outfile.open("all.xyz", ios::out | ios::app);
  while(infile.good())
  {
    string line;
    getline(infile, line);
    outfile << line << endl;
  }
  outfile.close();
  infile.close();  

  //CHARGE1 = test.xyz charge, ligand charge - structure that conformational search is done on.
  string cfilename = "CHARGE1";
  int charge = get_charge(cfilename);
  int natoms = get_natoms(filename);
  int N3 = natoms*3;
  printf("  there are %2i atoms \n",natoms);

  string* anames = new string[natoms];
  get_all_xyz(natoms,anames,xyzall,confile);
  
  int nstruct = xyzall.size();
  int* anumbers = new int[natoms];
  for (int i=0;i<natoms;i++)
    anumbers[i] = PTable::atom_number(anames[i]);


  E = new double[nstruct];
  for (int i=0;i<nstruct;i++) E[i] = 0.;
  
  opt_mopac(charge,natoms,anames,anumbers,xyzall,E,1);
  write_all_xyz(natoms,anames,E,xyzall,"all2.xyz");

#if 0
  printf("   showing all structures \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" %2i \n\n",natoms);
    for (int j=0;j<natoms;j++)
      printf(" %2s %8.5f %8.5f %8.5f \n",anames[j].c_str(),xyzall[i][3*j+0],xyzall[i][3*j+1],xyzall[i][3*j+2]);
  }
#endif

  delete [] anumbers;
  delete [] anames;

  return nstruct;
}

int main(int argc, char* argv[])
{
  printf("\n\n in main() \n");
  string xyzfile = "test.xyz";
  string targetfile = "target.xyz";
  switch (argc){
    case 2:
      //printf(" case 2 \n");
      xyzfile = argv[1];
      break;
    default:
      break;
  }
  // AD - Feb 16, 2017
  // Why is nconf hard coded?
  // We run procedure_2 because we just want to generate optimized ligand conformers?
  int nconf = 250000;
//  procedure_1(xyzfile,targetfile);
  procedure_2(nconf,xyzfile,targetfile);

  return 0;
}

void procedure_2(int nconf, string xyzfile, string targetfile)
{
  // AD - Feb 16, 2017
  // takes xyz file and runs conformer generation, and inital mopac opt
  ICoord ic1; 
  ic1.init(xyzfile);
  printf(" initial bonds \n");
  ic1.print_bonds();

  vector<double*> xyzall;
  double* E;
  // AD - Feb 16, 2017
  //why is nconf here so large in comparison to nconf for procedure 1?
  int nstruct = generate_conformers_and_opt(nconf,xyzfile,E,xyzall);

  return;
}

void procedure_1(string xyzfile, string targetfile)
{

  // AD - Feb 16, 2017
  // Runs gsm to push ligand and metal center together
  // sums up charges to get overall charge
  
  ICoord ic1; 
  ic1.init(xyzfile);
  printf(" initial bonds \n");
  ic1.print_bonds();

  vector<double*> xyzall;
  double* E;
  int nstruct = generate_conformers_and_opt(2000,xyzfile,E,xyzall);

  string cfilename = "CHARGE1";
  int charge1 = get_charge(cfilename);
  cfilename = "CHARGE2";
  int charge2 = get_charge(cfilename);
  int charget = charge1 + charge2;

  int natoms1 = get_natoms(xyzfile);
  int natoms2 = get_natoms(targetfile);
  string* anames = new string[natoms1];
  string* anamesm = new string[natoms2];
  int* anumbers = new int[natoms1];
  int* anumbersm = new int[natoms2];
  double* xyz0 = new double[3*natoms1];
  double* xyzm = new double[3*natoms2];
  xyz_read(natoms1,anames,xyz0,xyzfile);
  xyz_read(natoms2,anamesm,xyzm,targetfile);
  for (int i=0;i<natoms1;i++)
    anumbers[i] = PTable::atom_number(anames[i]);
  for (int i=0;i<natoms2;i++)
    anumbersm[i] = PTable::atom_number(anamesm[i]);

  int natomst = natoms1 + natoms2;
  int N3t = natomst*3;
  	string* anamest = new string[natomst];
  int* anumberst = new int[natomst];
  for (int i=0;i<natoms1;i++)
    anamest[i] = anames[i];
  for (int i=0;i<natoms2;i++)
    anamest[natoms1+i] = anamesm[i];
  for (int i=0;i<natomst;i++)
    anumberst[i] = PTable::atom_number(anamest[i]);

  printf("  target structure: \n");
  for (int i=0;i<natoms2;i++)
    printf(" %2s %8.5f %8.5f %8.5f \n",anamesm[i].c_str(),xyzm[3*i+0],xyzm[3*i+1],xyzm[3*i+2]);

  int* unique = new int[nstruct];

  align_and_opt(natoms1,natoms2,anames,anamesm,anamest,anumbers,anumbersm,charget,nstruct,unique,xyzall,xyzm);

 //retrieve firstnode.xyz files from GSM
  vector<double*> xyzalla;
  for (int i=0;i<nstruct;i++)
  {
    double* xyz1 = new double[3*natomst];
    xyzalla.push_back(xyz1);
  }
  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");
    string gsmfilename = "scratch/firstnode.xyz"+nstr;
    xyz_read_last(natomst,xyzalla[i],gsmfilename);
  }
  opt_mopac(charget,natomst,anamest,anumberst,xyzalla,E,2);
  write_all_xyz(natomst,anamest,E,xyzalla,"all4.xyz");

  printf(" done generating complexes \n");


  check_bonding(nstruct,ic1,natoms1,anames,anumbers,xyzalla,unique);
  printf("\n");

  // AD - Feb 16, 2017
  // Section of code chooses lowest E structures to add to all5.xyz
  // Not sure how it decides what the energy cutoff for strctures is
  // ** How does it decide?
  vector<pair<double,int> > bestc;
  for (int i=0;i<nstruct;i++)
  if (unique[i])
  {
     pair<double,int> newc(E[i],i);
     bestc.push_back(newc);
  }
  //sort(bestc.begin(),bestc.end());

  int nstructf = bestc.size();
  printf("  best structures: \n");
  for (int i=0;i<nstructf;i++)
  {
    printf("   %2i: %8.5f \n",bestc[i].second+1,bestc[i].first);
  }

  double* E1 = new double[nstructf];
  double** xyzb = new double*[nstructf];
  for (int i=0;i<nstructf;i++)
    xyzb[i] = new double[N3t];

  for (int i=0;i<nstructf;i++)
  {
    int i1 = bestc[i].second;
    E1[i] = E[i1];
    double* xyz1 = xyzalla[i1];
    for (int j=0;j<N3t;j++)
      xyzb[i][j] = xyz1[j];
  }

  write_all_xyz(natomst,anamest,nstructf,E1,xyzb,"all5.xyz");

  return;
}


void check_bonding(int nstruct, ICoord ic1, int natoms1, string* anames, int* anumbers, vector<double*> xyzalla, int* unique)
{
  // AD - Feb 16 2017
  // Is this section supposed to make sure the ligand stays intact?
  // Doesn't neccessarily work, a good number of my structures either broke or formed cycle structures
  printf("\n  now checking bonding \n");
  for (int i=0;i<nstruct;i++)
  {
    ICoord ic2;
    ic2.isOpt = 0;
    ic2.init(natoms1,anames,anumbers,xyzalla[i]);
    //ic2.print_bonds();

    if (ic2.nbonds!=ic1.nbonds)
    {
      printf("  struct %2i bonds changed, eliminating \n",i+1);
      unique[i] = 0;
    }
    else
    {
     //here check that bonding is the same
      int nf = 0;
      //for (int j=0;j<ic1.nbonds;j++)
    }
  }

  return;
}


int read_adds(int* adds, string addfile)
{
  // AD - Feb 16, 2017
  // would it make sense to read in the ADD moves separately, in other words look for
  // mono-, bi-, tri-, and tetra- dentate structures
  // rather than all ISOMERSXXXX files having all add moves present
  int nadd = 0;

  ifstream infile;
  infile.open(addfile.c_str());
  if (!infile){
    printf(" Error: couldn't open file: %s \n",addfile.c_str());
    exit(-1);
  } 
  
  string line;
  while (!infile.eof())
  {
    getline(infile,line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    if (tok_line.size()<2) break;

    adds[2*nadd+0] = atoi(tok_line[0].c_str())-1;
    adds[2*nadd+1] = atoi(tok_line[1].c_str())-1;
    //printf(" ADD found: %2i %2i \n",adds[2*nadd+0]+1,adds[2*nadd+1]+1);
    nadd++;
  }

  return nadd;
}


void align_and_opt(int natoms1, int natoms2, string* anames, string* anamesm, string* anamest, int* anumbers, int* anumbersm, int charget, int nstruct, int* unique, vector<double*> xyzall, double* xyzm)
{
  int natomst = natoms1 + natoms2;
  
  int* adds = new int[8]; //no more than tetradentate
  int nadd = read_adds(adds,"ADD");
  for (int i=0;i<nadd;i++)
  {
    adds[2*i+1] += natoms1;
    printf(" ADD found: %2i %2i \n",adds[2*i+0]+1,adds[2*i+1]+1);
  }

  int N3t = (natoms1+natoms2)*3;
  double** xyzalign = new double*[nstruct];
  for (int i=0;i<nstruct;i++)
    xyzalign[i] = new double[N3t];

  printf("  now creating initial.xyz \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" working on structure %2i \n",i+1);

    Align align1;
    align1.init(natoms1,anames,anumbers,xyzall[i],natoms2,anamesm,anumbersm,xyzm);
    align1.add_align(nadd,adds);
    for (int j=0;j<N3t;j++)
      xyzalign[i][j] = align1.xyza[j];
  }
  printf("\n");


  write_all_xyz(natomst,anamest,nstruct,NULL,xyzalign,"all3.xyz");
  write_gsm(natomst,anamest,charget,nstruct,NULL,xyzalign,nadd,adds);

#if 0
  printf(" printing aligned structures \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" %2i \n \n",natomst);
    for (int j=0;j<natomst;j++)
      printf(" %2s %8.5f %8.5f %8.5f \n",anamest[j].c_str(),xyzalign[i][3*j+0],xyzalign[i][3*j+1],xyzalign[i][3*j+2]);
  }
#endif

  do_gsm(nstruct);

  for (int i=0;i<nstruct;i++)
    delete [] xyzalign[i];
  delete [] xyzalign;

  return;
}

void do_gsm(int nstruct)
{
//  printf("  skipping GSM step \n");
//  return;

  string cmd;
  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");
    string istr = StringTools::int2str(i,1,"0");
    string filename = "scratch/firstnode.xyz"+nstr;
    struct stat sts;
    if (stat(filename.c_str(), &sts) != -1)
    {
      printf("  string %2i already done \n",i);
    }
    else
    {
      cmd = "./gfstringq.exe "+istr+" > scratch/paragsm"+nstr;
      system(cmd.c_str());
    }
  }
}

void get_all_xyz(int natoms, string* anames, vector<double*> &xyzs, string xyzfile)
{    
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file: %s \n",xyzfile.c_str());
    exit(-1);
  } 
  
  string line;
  int nf = 0;
  while (!infile.eof())
  {
    getline(infile, line);
    int natoms1 = atoi(line.c_str());
    if (natoms1==natoms)
    {
      double* coords = new double[3*natoms];
      getline(infile, line);
  
      for (int i=0;i<natoms;i++)
      {
        getline(infile, line);
        int length=StringTools::cleanstring(line);
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        if (nf==0)
          anames[i]=tok_line[0];
        coords[3*i+0]=atof(tok_line[1].c_str());
        coords[3*i+1]=atof(tok_line[2].c_str());
        coords[3*i+2]=atof(tok_line[3].c_str());
      }
      xyzs.push_back(coords);
      nf++;
    } //if adding new geom
    else
    {
      printf(" done reading after %2i structures \n",nf);
      break;
    }
  }
  infile.close();
 
  return;
}   

void xyz_read(int natoms, string* anames, double* coords, string xyzfile)
{  
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file \n");
    exit(-1);
  } 
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  printf(" natoms: %i \n",natoms);
  
  success=getline(infile, line);
  
  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    coords[3*i+0]=atof(tok_line[1].c_str());
    coords[3*i+1]=atof(tok_line[2].c_str());
    coords[3*i+2]=atof(tok_line[3].c_str());
  }
  
  infile.close();

#if 0
  printf(" XYZ: \n");
  for (int i=0;i<natoms;i++)
    printf(" %s %8.6f %8.6f %8.6f \n",anames[i].c_str(),coords[3*i+0],coords[3*i+1],coords[3*i+2]);
#endif

  printf(" done reading XYZ \n"); fflush(stdout);
 
  return;
}   

void xyz_read_last(int natoms, double* coords, string xyzfile)
{ 
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file: %s \n",xyzfile.c_str());
    exit(-1);
  } 
  
  string line;
  int nf = 0;
  while (!infile.eof())
  {
    getline(infile, line);
    int natoms1 = atoi(line.c_str());
    if (natoms1==natoms)
    {
      getline(infile, line);  
      for (int i=0;i<natoms;i++)
      {
        getline(infile, line);
        int length=StringTools::cleanstring(line);
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        coords[3*i+0]=atof(tok_line[1].c_str());
        coords[3*i+1]=atof(tok_line[2].c_str());
        coords[3*i+2]=atof(tok_line[3].c_str());
      }
      nf++;
    } //if adding new geom
    else
    {
      printf(" done reading after %2i structures \n",nf);
      break;
    }
  }
  infile.close();
 
  return;
}   

int get_charge(string filename)
{
  return get_natoms(filename);
}

int get_natoms(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf("  couldn't find file %s \n",filename.c_str());
    exit(-1);
  }

  string line;
  getline(infile, line);
  int length=StringTools::cleanstring(line);
  int natoms = atoi(line.c_str());

  infile.close();

  return natoms;
}

#if 0
int get_unique_conf(int nstruct, int* unique)
{
  int ns2 = nstruct * nstruct;
  int* similar = new int[ns2];
  for (int i=0;i<ns2;i++) similar[i] = 0;

  OBConversion obc1,obc2;
  obc1.SetInFormat("xyz");
  obc2.SetInFormat("xyz");
  OBMol mol1;
  OBMol mol2;

  int nf1 = 0;
  bool notatend1 = obc1.ReadFile(&mol1,"all2.xyz");
  while (notatend1)
  {
    //printf("  Molecular Weight: %5.3f \n",mol1.GetMolWt());

    int nf2 = 0;
    bool notatend2 = obc2.ReadFile(&mol2,"all2.xyz");
    while (notatend2)
    {
      OBAlign align1 = OBAlign(mol1,mol2);
      align1.Align();
      double dist = align1.GetRMSD();
      //printf(" difference (%2i-%2i): %10.6f \n",nf1,nf2,dist);
      if (dist<DIST_THRESH) similar[nf1*nstruct+nf2] = 1;

      mol2.Clear();
      notatend2 = obc2.Read(&mol2);
      nf2++;
    }
    mol1.Clear();
    notatend1 = obc1.Read(&mol1);
    nf1++;
  }

  for (int i=0;i<nstruct;i++) unique[i] = 1;
  for (int i=0;i<nstruct;i++)
  for (int j=0;j<i;j++)
  if (similar[i*nstruct+j])
    unique[j] = 0;

  int nf = 0;
  for (int i=0;i<nstruct;i++)
  if (unique[i])
    nf++;
  
  // AD - Feb 16, 2017
  // How does a similarity matrix work?
  printf("\n similarity matrix: \n");
  for (int i=0;i<nstruct;i++)
  {
    for (int j=0;j<nstruct;j++)
      printf(" %i",similar[i*nstruct+j]);
    printf("\n");
  }
  printf(" unique list:");
  for (int i=0;i<nstruct;i++)
    printf(" %i",unique[i]);
  printf("\n\n");

  return nf;
}
#endif

void write_all_xyz(int natoms, string* anames, double* E, vector<double*> xyzs, string xyzfile_string)
{
  int nstruct = xyzs.size();
  double** xyzs1 = &xyzs[0];
  write_all_xyz(natoms,anames,nstruct,E,xyzs1,xyzfile_string);
  return;
}

void write_all_xyz(int natoms, string* anames, int nstruct, double* E, double** xyzs, string xyzfile_string)
{
  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  char* sbuff = new char[3000];
  for (int i=0;i<nstruct;i++)
  {
    if (E==NULL)
      xyzfile << natoms << endl << endl;
    else
      xyzfile << natoms << endl << E[i] << endl;
    for (int j=0;j<natoms;j++)
    {
      sprintf(sbuff," %2s %10.6f %10.6f %10.6f \n",anames[j].c_str(),xyzs[i][3*j+0],xyzs[i][3*j+1],xyzs[i][3*j+2]);
      xyzfile << sbuff;
    }
  }

  xyzfile.close();

  delete [] sbuff;

  printf(" done writing \n");

  return;
}

void write_gsm(int natoms, string* anames, int charge, int nstruct, double* E, double** xyzs, int nadd, int* adds)
{
  char* sbuff = new char[3000];

  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");

    string xyzfile_string = "scratch/initial"+nstr+".xyz";
    ofstream xyzfile;
    xyzfile.open(xyzfile_string.c_str());
    xyzfile.setf(ios::fixed);
    xyzfile.setf(ios::left);
    xyzfile << setprecision(8);

    xyzfile << natoms << endl << charge << endl;
    for (int j=0;j<natoms;j++)
    {
      sprintf(sbuff," %2s %10.6f %10.6f %10.6f \n",anames[j].c_str(),xyzs[i][3*j+0],xyzs[i][3*j+1],xyzs[i][3*j+2]);
      xyzfile << sbuff;
    }
    xyzfile.close();

    string isofile_string = "scratch/ISOMERS"+nstr;
    ofstream isofile;
    isofile.open(isofile_string.c_str());
    isofile << "NEW" << endl;
    for (int j=0;j<nadd;j++)
      isofile << "BOND " << adds[2*j+0]+1 << " " << adds[2*j+1]+1 << endl;
    isofile.close();

    double dist1 = 2.0;
    string forcefile_string = "scratch/FORCE"+nstr;
    ofstream forcefile;
    forcefile.open(forcefile_string.c_str());
    for (int j=0;j<nadd;j++)
      forcefile << adds[2*j+0]+1 << " " << adds[2*j+1]+1 << " " << dist1 << " 0.15 " << endl;
    forcefile.close();
  }


  delete [] sbuff;

  printf(" done writing \n");

  return;
}


