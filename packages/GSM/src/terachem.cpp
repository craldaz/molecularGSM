#include "terachem.h"

using namespace std;

void TC::alloc(int natoms0)
{
  natoms = natoms0;
  anumbers = new int[natoms+1];
  anames = new string[natoms+1];

  return;
}

void TC::init(int natoms0,int nstates0, int nclosed0, int nactive0, string* anames0, double* xyz0, string basis0,string method0,int run, int rune)
{
  gradcalls = 0;
  nscffail = 0;
  firstrun = 1;

  natoms = natoms0;
  nclosed = nclosed0;
  nstates = nstates0;
  nactive = nactive0;
  basis = basis0;
  method = method0;
  ncpu = 1;

  anames = new string[natoms+1];
  xyz = new double[3*natoms];
  dvec = new double[3*natoms];
	dgrad = new double[3*natoms];
  grad = new double[3*natoms];

  for (int i=0;i<natoms0;i++)
  {
    anames[i] = anames0[i];
    xyz[3*i+0] = xyz0[3*i+0];
    xyz[3*i+1] = xyz0[3*i+1];
    xyz[3*i+2] = xyz0[3*i+2];
  }
  runNum = run;
  runend = rune;

#if TC
  printf("  TC initialized \n");
#endif

  string nstr = StringTools::int2str(runNum,4,"0");
  string runends = StringTools::int2str(runend,2,"0");
  outfile="tcout"+nstr+runends;

#if 0
  printf(" anames: ");
  for (int i=0;i<natoms;i++)
    printf(" %s",anames[i].c_str());
  printf("\n");
  printf(" anumbers: ");
  for (int i=0;i<natoms;i++)
    printf(" %i",anumbers[i]);
  printf("\n");
#endif

  return;
}


void TC::grads(int state,int runend,int runnum)
{
  int badgeom = check_array(3*natoms,xyz);
  if (badgeom)
  {
    printf(" ERROR: Geometry contains NaN, exiting! \n");
    exit(-1);
  }

  if (ncpu<1) ncpu = 1;

  for (int i=0;i<3*natoms;i++)
    grad[i] = 0.;

  int num,k,c;
  double V = -1.;

  string nstr=StringTools::int2str(runNum,4,"0");
  string runends=StringTools::int2str(runend,2,"0");
  string wstr=StringTools::int2str(state,1,"0");
  string endstr = nstr+"."+runends +"."+wstr;
  string molname = "scratch/structure"+endstr;

  ofstream geomfile(molname.c_str());
  //cout << " geomfile " << molname << endl;

  geomfile << natoms << endl << endl;
  for(int j=0;j<natoms;j++)
  {
    geomfile << setw(2) << anames[j];
    geomfile << setw(16)<< xyz[3*j+0];
    geomfile << setw(16)<< xyz[3*j+1];
    geomfile << setw(16)<< xyz[3*j+2] << endl;
  }
  geomfile.close();

  string gradfile = "scratch/GRAD"+endstr;
  string cmd;
  if (method == "FOMO")
  {
    string nactstr = StringTools::int2str(nactive,1,"0");
    string nststr = StringTools::int2str(nstates,1,"0");
    string nclstr = StringTools::int2str(nclosed,1,"0");
    cmd = "./grad.py  --method "+ method +" --fname " + endstr + " --active " + nactstr + " --nstates " + nststr + " --closed " + nclstr + " --basis " + basis + " --gradient --wstate " + wstr+ " > " + gradfile;
  }
  else if (method == "DFT")
    cmd = "./grad.py  --method "+ method +" --fname " + endstr + " --basis " + basis + " --gradient --wstate " + wstr+ " > " + gradfile;
 
  //cout << cmd << endl;
  system(cmd.c_str());
  system("wait");

  gradcalls++;

  return;
}

void TC::calc_dvec(int runend,int runnum)
{
  int badgeom = check_array(3*natoms,xyz);
  if (badgeom)
  {
    printf(" ERROR: Geometry contains NaN, exiting! \n");
    exit(-1);
  }

  if (ncpu<1) ncpu = 1;

  for (int i=0;i<3*natoms;i++)
    grad[i] = 0.;

  int num,k,c;
  double V = -1.;
  cout << " Wstates are " << wstate << " " << wstate2 << endl;
  string nstr=StringTools::int2str(runNum,4,"0");
  string runends=StringTools::int2str(runend,2,"0");
  string wstr1=StringTools::int2str(wstate,1,"0");
  string wstr2=StringTools::int2str(wstate2,1,"0");
  string endstr = nstr+"."+runends +"."+wstr1 +wstr2;
  string molname = "scratch/structure"+endstr;

  ofstream geomfile(molname.c_str());
  cout << " geomfile " << molname << endl;

  // print the molecule coordinate section
  geomfile << natoms << endl << endl;
  for(int j=0;j<natoms;j++)
  {
    geomfile << setw(2) << anames[j];
    geomfile << setw(16)<< xyz[3*j+0];
    geomfile << setw(16)<< xyz[3*j+1];
    geomfile << setw(16)<< xyz[3*j+2] << endl;
  }
  geomfile.close();

  string gradfile = "scratch/COUP"+endstr;
  string cmd;
  if (method == "FOMO")
  {
    string nactstr = StringTools::int2str(nactive,1,"0");
    string nststr = StringTools::int2str(nstates,1,"0");
    string nclstr = StringTools::int2str(nclosed,1,"0");
    cmd = "./grad.py  --method "+ method +" --fname " + endstr + " --active " + nactstr + " --nstates " + nststr + " --closed " + nclstr + " --basis " + basis + " --coupling --wstate " + wstr1 + " " + wstr2 + " > " + gradfile;
  }
  else if (method == "DFT")
	{
     printf(" No coupling in DFT!\n");
     exit(-1);
  }
 
  //cout << cmd << endl;
  system(cmd.c_str());
  system("wait");

  //printf(" done with TC grad call \n");  

  gradcalls++;

  return;
}

double TC::get_energy_grad(int choosestate,int runend,int runnum, double* grad)
{
  string nstr=StringTools::int2str(runNum,4,"0");
  string runends=StringTools::int2str(runend,2,"0");
  string wstr=StringTools::int2str(choosestate,1,"0");
  string endstr = nstr+"."+runends +"."+wstr;
  string file = "scratch/GRAD"+endstr;

  ifstream gradfile;
  gradfile.open(file.c_str());
  if (!gradfile)
  {
    printf(" Error opening gradient file! %s \n",file.c_str());
    nscffail+=1;
    return 0.0;
  }

  string line;
  bool success = true;
  success=getline(gradfile, line);
  double V;
  if (success)
  {
    V= atof(line.c_str());
    //printf(" found E: %7.5f \n",V);
  }
  else
    return 0.0;

  for (int i=0;i<natoms;i++)
  {
    if(gradfile.eof())
    {
      printf(" missing data in GRAD(1) \n"); fflush(stdout);
      nscffail+=1;
      grad[3*i+0] = grad[3*i+1] = grad[3*i+2] = 1.;
      break;
    }
    success=getline(gradfile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t[]");
    if (tok_line.size()<3)
    {
      printf(" missing data in GRAD(2) \n"); fflush(stdout);
      nscffail+=1;
      grad[3*i+0] = grad[3*i+1] = grad[3*i+2] = 1.;
    }
    else
    {
      grad[3*i+0] = atof(tok_line[0].c_str());
      grad[3*i+1] = atof(tok_line[1].c_str());
      grad[3*i+2] = atof(tok_line[2].c_str());
    }
  } //loop i over natoms

#if 0
  cout << " gradient: " << endl;
  for (int i=0;i<natoms;i++) 
    cout << grad[3*i+0] << " " << grad[3*i+1] << " " << grad[3*i+2] << endl;
#endif

  if (nscffail>25)
  {
    printf("\n\n Too many SCF failures: %i, exiting \n",nscffail);
    exit(1);
  }

  gradfile.close();

  return V*627.5; 
}

void TC::get_dvec(int runend,int runnum, double* dvec)
{
  string nstr=StringTools::int2str(runNum,4,"0");
  string runends=StringTools::int2str(runend,2,"0");
  string wstr1=StringTools::int2str(wstate,1,"0");
  string wstr2=StringTools::int2str(wstate2,1,"0");
  string endstr = nstr+"."+runends +"."+wstr1+wstr2;
  string file = "scratch/COUP"+endstr;

  ifstream coupfile;
  coupfile.open(file.c_str());
  if (!coupfile)
  {
    printf(" Error opening gradient file! %s \n",file.c_str());
    nscffail+=1;
    return;
  }

  string line;
  bool success = true;
  success=getline(coupfile, line);
  double V;
  if (success)
  {
    V= atof(line.c_str());
  }
  else
  {
    nscffail+=1;
    return;
  }

  for (int i=0;i<natoms;i++)
  {
    if(coupfile.eof())
    {
      printf(" missing data in COUP(1) \n"); fflush(stdout);
      nscffail+=1;
      dvec[3*i+0] = dvec[3*i+1] = dvec[3*i+2] = 1.;
      break;
    }
    success=getline(coupfile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t[]");
    if (tok_line.size()<3)
    {
      printf(" missing data in COUP(2) \n"); fflush(stdout);
      nscffail+=1;
      dvec[3*i+0] = dvec[3*i+1] = dvec[3*i+2] = 1.;
    }
    else
    {
      dvec[3*i+0] = atof(tok_line[0].c_str());
      dvec[3*i+1] = atof(tok_line[1].c_str());
      dvec[3*i+2] = atof(tok_line[2].c_str());
    }
  } //loop i over natoms

#if 1
  cout << " dvec: " << endl;
  for (int i=0;i<natoms;i++) 
    cout << dvec[3*i+0] << " " << dvec[3*i+1] << " " << dvec[3*i+2] << endl;
#endif

  if (nscffail>25)
  {
    printf("\n\n Too many SCF failures: %i, exiting \n",nscffail);
    exit(1);
  }

  coupfile.close();

  return; 
}

void TC::write_xyz_grad(double* coords, double* grad, string filename)
{
  ofstream xyzfile;
  string xyzfile_string = filename+".xyz";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  xyzfile << natoms << endl << energy << endl;
  for (int i=0;i<natoms;i++)
    xyzfile << anames[i] << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2] << endl;

  xyzfile.close();

  ofstream gradfile;
  string gradfile_string = filename+".grad";
  gradfile.open(gradfile_string.c_str());
  gradfile.setf(ios::fixed);
  gradfile.setf(ios::left);
  gradfile << setprecision(6);

  gradfile << natoms << endl << energy << endl;
  for (int i=0;i<natoms;i++)
    gradfile << anames[i] << " " << grad[3*i+0] << " " << grad[3*i+1] << " " << grad[3*i+2] << endl;

  gradfile.close();

  return;
}

void TC::reset(double* xyz1)
{
  for (int i=0;i<3*natoms;i++)
    xyz[i] = xyz1[i];
  return;
}
