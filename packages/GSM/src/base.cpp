#include "base.h"

void Base::init(string infilename, int run, int nprocs)
{

  string nstr=StringTools::int2str(run,4,"0");
  string xyzfile = "scratch/initial"+nstr+".xyz";
  ncpu = nprocs;
  runNum = run;
  runends = nstr;

  //parameter init from infileq
  parameter_init(infilename);

  //read xyz file
  //structure_init(xyzfile);

  //initialize ISOMER if necessary

   return;
}

void Base::driver()
{

    //preopt if turned on
    printf("Hello\n");
    //if MECP 
      //set variables, parameters, etc
      //do opt,etc

   return;
}


void Base::parameter_init(string infilename)
{

  printf(" Hello 2\n");

//  nnmax = 9;
//  hessSSM = 0; //default no Hessian given
//  tstype = 0;
//  prodelim = 1000.;
//  lastOpt = 0; //default do not optimize last SSM node 
//  initialOpt = 0; //default do not optimize first node 
//  DQMAG_SSM_MAX = 0.8; 
//  DQMAG_SSM_MIN = 0.2; 
//  QDISTMAX = 5.0;
//  PEAK4_EDIFF = 2.0;
//  isRestart = 0;
//  isSSM = 0;
//  isFSM = 0;
//	isMECI =0;
//	isPRODUCT =0;
//	isMAP_DE =0;
//	isMAP_SE=0;
//	isSE_ESSM =0; //these are for excited state methods
//  use_exact_climb = 2;
//	restart_wfn=0;  //restart molpro calc with wfu
//
//  printf("Initializing Tolerances and Parameters... \n");
//  printf("  -Opening %s \n",infilename.c_str());
//
//  ifstream infile;
//  infile.open(infilename.c_str());
//  if (!infile){
//    printf("\n Error opening %s \n",infilename.c_str());
//    exit(-1);
//  }
//
//  printf("  -reading file... \n");
//
//  // pass infile to stringtools and get the line containing tag
//  string tag="String Info";
//  bool found=StringTools::findstr(infile, tag);
//  if (!found) { cout << "Could not find tag for Default Info" << endl; exit(-1);}
//  string line, templine, tagname;
//
//  // parse the input section
//  bool stillreading=true;
//  while (stillreading)
//  {
//    stillreading=false;
//    // set to false and set back to true if we read something
//    // get filename
//    getline(infile, line);
//    vector<string> tok_line = StringTools::tokenize(line, " ,\t");
//    templine=StringTools::newCleanString(tok_line[0]);
//    tagname=StringTools::trimRight(templine);
//    // these variables are denoted by strings with same name
//    if (tagname=="MAX_OPT_ITERS") {
//      MAX_OPT_ITERS=atoi(tok_line[1].c_str());
//      stillreading=true;
//      cout <<"  -MAX_OPT_ITERS: " << MAX_OPT_ITERS << endl;
//    }
//    if (tagname=="STEP_OPT_ITERS") {
//      STEP_OPT_ITERS = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -STEP_OPT_ITERS: " << STEP_OPT_ITERS << endl;
//    }
//    if (tagname=="RESTART") {
//      isRestart = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -RESTART: " << isRestart << endl;
//    }
//    if (tagname=="TS_FINAL_TYPE") {
//      tstype = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -TS_FINAL_TYPE: " << tstype << endl;
//      if (tstype!=0 && tstype!=1 && tstype!=2)
//      {
//        printf("  TS_FINAL_TYPE must be 0, 1 or 2 \n");
//        exit(1);
//      }
//    }
//    if (tagname=="PRODUCT_LIMIT") {
//      prodelim = atof(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -PRODUCT_LIMIT: " << prodelim << endl;
//    }
//    if (tagname=="FINAL_OPT") {
//      lastOpt = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -FINAL_OPT: " << lastOpt << endl;
//    }
//    if (tagname=="INITIAL_OPT") {
//      initialOpt = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -INITIAL_OPT: " << initialOpt << endl;
//    }
//    if (tagname=="SSM_DQMAX") {
//      DQMAG_SSM_MAX = atof(tok_line[1].c_str());
//      DQMAG_SSM_MIN = DQMAG_SSM_MAX/4.;
//      stillreading = true;
//      cout <<"  -SSM_DQMAX: " << DQMAG_SSM_MAX << endl;
//      cout <<"  -SSM_DQMIN: " << DQMAG_SSM_MIN << endl;
//    }
//    if (tagname=="MIN_SPACING") {
//      QDISTMAX = atof(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -SSM_MIN_SPACING: " << QDISTMAX << endl;
//    }
//    if (tagname=="INT_THRESH") {
//      PEAK4_EDIFF = atof(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -INT_THRESH: " << PEAK4_EDIFF << endl;
//    }
//	  if (tagname=="RESTART_WFN"){
//			restart_wfn = atoi(tok_line[1].c_str());
//			stillreading=true;
//			cout <<"  -restart_wfn = " << restart_wfn << endl; 
//		}
//    if (tagname=="SM_TYPE") {
//      if (tok_line[1]=="SSM")
//      {
//        printf("  -using SSM \n");
//        isFSM = 0;
//        isSSM = 1;
//      }
//      else if (tok_line[1]=="SSMFR")
//      {
//        printf("  -using SSM with initial Hessian \n");
//        isFSM = 0;
//        isSSM = 1;
//        hessSSM = 1;
//      }
//      else if (tok_line[1]=="FSM")
//      {
//        printf("  -using FSM \n");
//        isFSM = 1;
//        isSSM = 0;
//      }
//      else if (tok_line[1]=="OPT")
//      {
//        printf("  -using OPT \n");
//        isFSM = -1;
//        isSSM = -1;
//      }
//			else if (tok_line[1]=="MECI")
//			{
//				printf("  -using MECI \n");
//				isMECI = 1;
//			}
//			else if (tok_line[1]=="PRODUCT")
//			{
//				printf("  -using PRODUCT\n");
//				isPRODUCT = 1;
//			}
//			else if (tok_line[1]=="SE-ESSM")
//			{
//				printf("  -using SE-ESSM \n");
//				isSSM = 1; //both are turned on 
//				isSE_ESSM=1;
//			}
//			else if (tok_line[1]=="DE-MAP")
//			{
//				printf("  -using DE-MAP \n");
//				isMAP_DE=1;
//			}
//			else if (tok_line[1]=="SE-MAP")
//			{
//				printf("  -using SE-MAP \n");
//				isMAP_SE=1;
//			}
//      else
//      {
//        printf("  -using GSM \n");
//        isFSM = 0;
//        isSSM = 0;
//      }
//      stillreading=true;
//    }
//    if (tagname=="CONV_TOL") {
//      CONV_TOL=atof(tok_line[1].c_str());
//      stillreading=true;
//      cout <<"  -CONV_TOL = " << CONV_TOL << endl;
//    }
//    if (tagname=="ADD_NODE_TOL"){
//      ADD_NODE_TOL=atof(tok_line[1].c_str());
//      stillreading=true;
//      cout <<"  -ADD_NODE_TOL = " << ADD_NODE_TOL << endl;
//    }
//    if (tagname=="SCALING"){
//      SCALING = atof(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -SCALING = " << SCALING << endl;
//    }
//    if (tagname=="BOND_FRAGMENTS"){
//      bondfrags = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -BOND_FRAGMENTS = " << bondfrags << endl;
//    }
//    if (tagname=="GROWTH_DIRECTION"){
//      GROWD = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -GROWTH_DIRECTION = " << GROWD << endl;
//    }
//    if (tagname=="CLIMB_TS"){
//      use_exact_climb = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -CLIMB_TS = " << use_exact_climb << endl;
//    }
//    if (tagname=="nnodes" || tagname=="NNODES"){
//      nnmax = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -NNODES = " << nnmax << endl;
//      if (nnmax < 3)
//      {
//        printf("\n\n ERROR: NNODES cannot be less than 3 \n");
//        exit(1);
//      }
//    }
//
//  } //while stillreading
//  infile.close();
// 
//  if (tstype==2)
//  {
//    printf("  TS_FINAL_TYPE == 2, turning off climbing image and TS search \n");
//    use_exact_climb = 0;
//  }
//  nnmax0 = nnmax;
//
//  printf(" Done reading inpfileq \n\n");
}
