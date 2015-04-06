#include "cctrl.hpp"
__LOGOUT

CCTRL::CCTRL() {}

ofstream& operator<< (ofstream& fout, CCTRL &ctrl)
{
	fout << ctrl.clIf;
	WRITEDATA(fout, ctrl.nDS);
	WRITEDATA(fout, ctrl.nIS);
	WRITEDATA(fout, ctrl.nMode);
	fout << ctrl.nAaAtEaDsSt;
	fout << ctrl.dsSt;
	fout << ctrl.rvDsSt;
	fout << ctrl.nSmStOfCl;
	fout << ctrl.clSz;
	
	WRITEDATA(fout, ctrl.cosChiCut);

	WRITEDATA(fout, ctrl.mnNSeq);
	WRITEDATA(fout, ctrl.facDF2);
	WRITEDATA(fout, ctrl.elimMode);
	WRITEDATA(fout, ctrl.nElimCFs);
	WRITEDATA(fout, ctrl.numSeqs);
	WRITEDATA(fout, ctrl.facDF1);
	WRITEDATA(fout, ctrl.nMScr);
	WRITEDATA(fout, ctrl.nSrchRemCFs);
	WRITEDATA(fout, ctrl.useHiarc);
	return fout;
}

ifstream& operator>> (ifstream& fin, CCTRL &ctrl)
{
	fin >> ctrl.clIf;
	READDATA(fin, ctrl.nDS);
	READDATA(fin, ctrl.nIS);
	READDATA(fin, ctrl.nMode);
	fin >> ctrl.nAaAtEaDsSt;
	fin >> ctrl.dsSt;
	fin >> ctrl.rvDsSt;
	fin >> ctrl.nSmStOfCl;
	fin >> ctrl.clSz;

	READDATA(fin, ctrl.cosChiCut);

	READDATA(fin, ctrl.mnNSeq);
	READDATA(fin, ctrl.facDF2);
	READDATA(fin, ctrl.elimMode);
	READDATA(fin, ctrl.nElimCFs);
	READDATA(fin, ctrl.numSeqs);
	READDATA(fin, ctrl.facDF1);
	READDATA(fin, ctrl.nMScr);
	READDATA(fin, ctrl.nSrchRemCFs);
	READDATA(fin, ctrl.useHiarc);
	return fin;
}

ostream& operator<< (ostream& cout, CCTRL &ctrl)
{
	char tmpstr[MAXLINE];

	cout << "### control options ###\n";
	cout << "DesignSites           "   << ctrl.fnDSites       <<"\n";
//	cout << "Exec                  "   << ctrl.fnExec         <<"\n";
	cout << "SequenceFile          "   << ctrl.fnSeqFile      <<"\n";
	cout << "ResultsFile           "   << ctrl.fnEciFile      <<"\n";
	cout << "LogFile               "   << ctrl.fnOutPut       <<"\n";
	cout << "DumpFile              "   << ctrl.fnOutPutDump   <<"\n";
	cout << "NumOfSeqs             "   << ctrl.numSeqs        <<"\n";
	cout << "MinNSeq               "   << ctrl.mnNSeq         <<"\n";
//	cout << "FacDF1                "   << ctrl.facDF1         <<"\n";
	cout << "Parameter-d           "   << ctrl.facDF2         <<"\n";
//	cout << "ElimMode              "   << ctrl.elimMode       <<"\n";
//	cout << "CutCosChi             "   << ctrl.cosChiCut      <<"\n";
	cout << "NumElimCFs            "   << ctrl.nElimCFs       <<"\n";
	cout << "Higher-orderMethod    "   << ctrl.useHiarc       <<"\n";

	if(!ctrl.prDesInfo) return cout;

	cout << "\n### Design Information ###" << endl;
	cout << "Num: index of design sites" << endl;
	cout << "nDSeq: number of amino acids for each design site" << endl;
	cout << "Number of design sites: " << ctrl.nDS << endl;
	sprintf(tmpstr," %3s %5s   %s","Num","nDSeq","List of amino acids");
	cout << tmpstr << endl;
	int n, m, m1, m2;
	for(n=0;n<ctrl.nDS;n++) {
		sprintf(tmpstr," %3d %5d  ",n+1,ctrl.nAaAtEaDsSt[n]);
		for(m=0;m<ctrl.nAaAtEaDsSt[n];m++)
			sprintf(tmpstr,"%s %2c",tmpstr,ctrl.dsSt(n,m)+'A'-1);
		cout << tmpstr << endl;
	}
	cout << "\nNum: arbitrary index of clusters" << endl;
	cout << "Size: cluster size denote single, pair, triplet, and higher-order clusters" << endl;
	cout << "nSym: number of physically-equivalent clusters which are listed in  " << endl;
	cout << "      the cluster information field" << endl;
	cout << "Cl Info: list of equivalent clusters " << endl;
	cout << "Number of total clusters: " << ctrl.nIS << endl;

	sprintf(tmpstr," %3s %4s %4s   %s","Num","Size","nSym","Cl Info");
	cout << tmpstr << endl;
	for(n=0;n<ctrl.nIS;n++) {
		sprintf(tmpstr," %3d %4d %4d  ",n+1,ctrl.clSz[n],ctrl.nSmStOfCl[n]);
		for(m=0;m<ctrl.nSmStOfCl[n];m++){
			for(m1=0;m1<ctrl.clSz[n];m1++){
				m2 = m*ctrl.clSz[n]+m1;
				sprintf(tmpstr,"%s %3d",tmpstr,ctrl.clIf(n,m2));
			}
			sprintf(tmpstr,"%s;",tmpstr);
		}
		cout << tmpstr << endl;
	}
	return cout;
}

void CCTRL::InitOpts()
{
	int ni;
	// Input list, Maximum 100 variables, If you want to increase number of input variable
	//        then change MAXINPUTITEM in general.hpp
	// ("Identifier in file", "varable", <mode>, "default value", "update")
	// <mode> 0:ordinary, 1:start from some specific position
	// update: true(if the value changes then start program from the begining), 
	//         false(no change in starting position)
	ni = -1;
	// File input and output
	lstInput[++ni](        "designsites",    &fnDSites, 0,  "tmp_DesignFile.dat", false);
	lstInput[++ni](               "exec",      &fnExec, 0,             "Command", false);
	lstInput[++ni](            "logfile",    &fnOutPut, 0,      "tmpResults.log", false);
	lstInput[++ni](           "dumpfile",&fnOutPutDump, 0,     "tmpResults.dump", false);
	lstInput[++ni](       "sequencefile",   &fnSeqFile, 0, "tmpSequenceFile.dat", false);
	lstInput[++ni](        "resultsfile",   &fnEciFile, 0,      "tmpResFile.dat", false);
	lstInput[++ni](         "saveoption",  &saveOption, 0,                     0, false);
	// Main control
	lstInput[++ni](               "redo",        &reDo, 0,                     0, false);
	lstInput[++ni](               "mode",       &nMode, 0,                     0, false);
	lstInput[++ni](            "mhighcf",     &mHighCF, 0,                     1, false);
	lstInput[++ni](       "modescrcvrms",       &nMScr, 0,                     1,  true);
	lstInput[++ni](          "numofseqs",     &numSeqs, 0,                  3000,  true);

	lstInput[++ni](            "minnseq",      &mnNSeq, 0,                     5,  true);
	lstInput[++ni](          "cutcoschi",   &cosChiCut, 1,                   1.0,  true);
	lstInput[++ni](             "facdf1",      &facDF1, 0,                   0.0,  true);
	lstInput[++ni](             "facdf2",      &facDF2, 0,                   0.0,  true);
	lstInput[++ni](           "elimmode",    &elimMode, 1,                     1,  true);
	lstInput[++ni](           "nelimcfs",    &nElimCFs, 1,                     0,  true);
	lstInput[++ni](        "nsrchremcfs", &nSrchRemCFs, 1,                   500,  true);
	lstInput[++ni](           "usehiarc",    &useHiarc, 0,                     0,  true);
	// print options
	lstInput[++ni](     "pr_design_info",   &prDesInfo, 0,                     1, false);
	lstInput[++ni](      "pr_log_energy",     &prLogEn, 0,                     1, false);
	lstInput[++ni](         "pr_log_ipr",     &prLogPn, 0,                     1, false);
	lstInput[++ni](         "pr_rem_cfs",    &prRemCFs, 0,                     1, false);
	lstInput[++ni](           "pr_cvrms",     &prCVRMS, 0,                     1, false);
	lstInput[++ni](         "pr_def_cfs",    &prDefCFs, 0,                     1, false);
	lstInput[++ni](           "pr_train",     &prTrain, 0,                     1, false);
	lstInput[++ni](        "pr_cedetail",  &prCEDetail, 0,                     1, false);

	// MCSimulation
	lstInput[++ni](              "ndata",       &nData, 0,                  1000, false);
	lstInput[++ni](           "nprestep",    &nPreStep, 0,                  1000, false);
	lstInput[++ni](          "ncoolstep",   &nCoolStep, 0,                  1000, false);
	lstInput[++ni](               "temp",       &dbl_T, 0,                 300.0, false);
	lstInput[++ni](               "code",       &codeN, 0,                  "mc", false);
	lstInput[++ni](           "nintstep",    &nIntStep, 0,                   300, false);
	lstInput[++ni](             "sttemp",     &dbl_stT, 0,                   0.0, false);
	lstInput[++ni](             "ndtemp",     &dbl_ndT, 0,                3000.0, false);
	lstInput[++ni](            "inctemp",    &dlb_incT, 0,                 100.0, false);
	lstInput[++ni](      "potentialfile",     &potFile, 0,        "Potential.ce", false);
	lstInput[++ni](         "outputfile",  &fnMCOutPut, 0,      "tmpResults.dat", false);
}

void CCTRL::ReadCtrl(FILE *fp)
{
	int nn, i, k;
	CSTRING str;

	// Connect and initialize options
	InitOpts();

	// Read control file
	nn = sizeof(lstInput)/(sizeof(CINPUT)); // total number of table items.
	while (str.ReadLine(fp)) {
		k = str.SplitLine(" ,");
		if(k<2) continue;
		for(i=0;i<nn;i++) 
			if(str(0).Lower() == lstInput[i].idName) break;
		if(i == nn) continue;
		lstInput[i].SetValue(str(1));
	}
}

bool CCTRL::CheckECIFile(char * file)
{
	ifstream fin(file,std::ios::in|std::ios::binary);
	if(!fin.is_open()) return false;

	int len;
	READDATA(fin, len);
	if(len > 10) {return false;}
	else {fin.seekg(0,std::ios::beg);}

	CSTRING tmpstr;
	fin >> tmpstr;
	if(!tmpstr.IsStrIn("ver")) return false;
	return true;
}

bool CCTRL::CheckECIFile(CSTRING & file)
{
	ifstream fin(file.GetString(),std::ios::in|std::ios::binary);
	if(!fin.is_open()) return false;

	int len;
	READDATA(fin, len);
	if(len > 10) {return false;}
	else {fin.seekg(0,std::ios::beg);}

	CSTRING tmpstr;
	fin >> tmpstr;
	if(!tmpstr.IsStrIn("ver")) return false;
	return true;
}

void CCTRL::ReadDesign(FILE *fp)
{
	int i, k;
	CSTRING str;

	// Read design sites file
	while(str.ReadLine(fp)) {
		if(str == "#design_start") break;
	}
	if(str != "#design_start") {__ERRORTHROW(ERR_DSNFORM); }
	int n=0;
	int ncol=27;
	while (str.ReadLine(fp)) {
		if(str == "#design_end") break;
		k = str.SplitLine(" ");
		if (k == 0) continue;
		n++;
		dsSt.SetDim(n,27);
		rvDsSt.SetDim(n,27);
		nAaAtEaDsSt.SetDim(n); // increase vector dimension to matrix dimension
		nAaAtEaDsSt[n-1] = k-1;
		for(i=1;i<k;i++){
			int nn = (int)(str(i)[0] - 'A'+1); // meaningful sequences start from 1
			                                     // nn =1 : A, 2: B, 3: C, 4: D
			dsSt(n-1,i-1) = nn;    // interaction range = (1 ~ nAaAtEaDsSt-1) 
			                       // Let dsSt(0,:)={1, 2, 3} then rvDsSt(0,1):0, (0,2):1, (0,3):2 
			rvDsSt(n-1,nn) = i-1;  // 0 means no interaction
		}
	}
	nDS = n;
	if(str != "#design_end") {__ERRORTHROW(ERR_DSNFORM); }

	// Read a file containing information for interaction between design sites
	rewind(fp);
	while(str.ReadLine(fp)) {
		if(str == "#cluster_start") break;
	}
	if(str != "#cluster_start") {__ERRORTHROW(ERR_DSNFORM); }
	n = 0;
	while (str.ReadLine(fp)) {
		if(str == "#cluster_end") break;
		k = str.SplitLine(";");
		int l = str(0).SplitLine(" ");
		if (k == 0) continue;
		if(k*l > clIf.Len(COLUMN)) {   // if more column is needed then increase column and add 1 row
			if(n == 0) clIf.SetDim(1,k*l);
			else clIf.SetDim(n+1,k*l);
		} else 
			clIf.SetDim(n+1,clIf.Len(COLUMN));
		nSmStOfCl.SetDim(n+1);
		clSz.SetDim(n+1);
		nSmStOfCl[n] = k;
		clSz[n] = l;
		for(i=0;i<k;i++) {  // Save each symmetric site
			CSTRING str1;
			str1 = str(i);
			int l2 = str1.SplitLine(" \t");
			int j;
			if(l != l2) {__ERRORTHROW(ERR_DSNFORM); }
			for(j=0;j<l2;j++) {  // Save each residue site of a cluster
				CSTRING str2;
				str2 = str1(j);
				int clnum = str2.CvtInt();
				if (clnum > nDS || clnum < 1) { 
					cout << "\nCheck design file" << endl;
					cout << "Error line: " << str << endl; 
					cout << "Valid site range of cluster is from 1 to " << nDS << endl;
					__ERRORTHROW(ERR_DSNFORM); 
				}
				clIf(n, i*l2+j) = clnum;
			}
		}
		n++;
	}
	nIS = n;
	if(str != "#cluster_end") {__ERRORTHROW(ERR_DSNFORM); }
}

void CCTRL::Usage(int mode) 
{
	char usagecontrol[][MAXLINE]={
		"\nEx) CONTROL FILE ",
		"--- Control file ---",
		"DesignSites   Design.dat                  // information of design sites ", 
		"Exec          CalcProtEn.exe              // name of energy model ",
		"SequenceFile  Sequence.dat ",
		"ResultsFile   TrainedResults.dat          // ECI train results  ",
		"LogFile       Results.log                 // result file     ",
		"DumpFile      Results.dat.dump            // contains information for recalculation",
		"                                                       ",
		"# determine to start from saved file (DumpFile) or not ",
		"Redo   1                                      ",
		"# mode : 1: extract eci values from sequence file  ",
		"#        2: extract eci vaules from random generation",
		"#        3: extract eci values from random generation with limitation",
		"Mode   2           ",
		"# total number of sequences to be used (working with option 0  1)     ",
		"NumOfSeqs          3000 ",
		"# cut off value of eci rmsd value  ",
		"CutPR             100.0     ",
		"CutSeqID            0.5     ",
		"nECIConv             10   ",
		"MaxNEci             500  // maximum number of updating eci values.   ",
		"                      ",
		"CutCF              0.10    ",
		"CutMaxNCF           500    // maximum number of updating cluser functions ",
		"                      ",
		"Screen_Interactions   1",
		"CutCVRMS            2.0",
		"                       ",
		"Remain_Freedom      2.0",
		"                       ",
		"SaveOption            1",
		"                       ",
		"PrSeqId               0",
		"                       ",
		"# minimum number of cluster functions included in initial sequence set ",
		"MinNCF                5      // always larger than 2 ",
		"                             ",
		"# weighting factor for degrees of freedom ",
		"FacDF                 0                   ",
		"                                          ",
		"# cluster function elimination method     ",
		"# elimMode                                ",
		"#   0: using cutCosChi (eliminate CFs with Cos(X) larger than the cut-off value), ",
		"#   1: using nElimCFs (elimininate CFs nElimCFs times).",
		"ElimMode              1      ",
		"CutCosChi           1.0      ",
		"nElimCFs              3      ",
		"                             ",
		"# Print options              ",
		"# Print design information in the header (default 0)",
		"// PR_Design_Info     1                             ",
		"# Print calculated energy in the course of evaluation (default 0)",
		"// PR_Log_Energy      1                                          ",
		"# Print participation number in the optimization routine (default 0)",
		"// PR_Log_PR          1                                    ",
		"# Print parameters of excluded CFs (default 0)             ",
		"// PR_Rem_CFs         1                                    ",
		"# Don't print CV score (default 1)                         ",
		"// PR_CVRMS           0                                    ",
		"# Don't print definition of CFs (default 1)                ",
		"// PR_Def_CFs         0                                    ",
		"# Print training sequences (default 0)                     ",
		"// PR_Train           1                                    ",
		"# Print details of CE in the individual updating CFs (default 0)  ",
		"// PR_CEDetail        1      "
	};
	char usagedesign[][MAXLINE]={
		"\nEx) DESIGN FILE ",
		"--- Design file ---",
		"this file contains dimeric coiled coil interaction and selected design sites ",
		"Possbile letters are A~Z, capital letter.                                    ",
		"                                                                             ",
		"#design_start                                                                ",
		"1   A L I V T M N R Q Y K                                                    ",
		"2   A E Q K R L N S H T D                                                    ",
		"3   A E Q K R L I H V T D                                                    ",
		"4   A L V I T N K R M C E                                                    ",
		"5   A L I V T M N R Q Y K                                                    ",
		"6   A E Q K R L N S H T D                                                    ",
		"7   A E Q K R L I H V T D                                                    ",
		"8   A L V I T N K R M C E                                                    ",
		"9   A L I V T M N R Q Y K                                                    ",
		"10  A E Q K R L N S H T D                                                    ",
		"#design_end                                                                  ",
		"                                                                             ",
		"#cluster_start                                                           ",
		"1; 5;  9;                                         ",
		"2; 6; 10;                                         ",
		"3; 7;                                         ",
		"4; 8;                                         ",
		"1 2; 5  6;  9 10;           ",
		"1 3; 5  7;           ",
		"1 4; 5  8;             ",
		"1 5; 5  9;                    ",
		"2 3; 6  7;            ",
		"#cluster_end                                                             "
	};
	char usageseqfile[][MAXLINE]={
	  "\nEx) SEQUENCE FILE ",
	  " Example of 20 design sites" 
	  " [Energy] [Sequence]",
  	  "--- Sequence file ---",
      "-197.000   S L T S K T N F A A A S I S L V K S E I",
      "-213.400   K T Q Y K T E S T S S T Q T A N F A E A",
      "-198.300   M S Q I S K T L M K A T V V D N L D V M",
      "-209.300   F L V F S K L T I N E I V Q S A S E Q T",
      "-213.200   N T S Y L S V S L Q L S T S N F T Q N I",
      "-213.600   M A E Y N N S A L Q S F V E D M F S E I",
	};
	int n, nn;
	if(mode == 0 || mode == 1) {
		nn = sizeof(usagecontrol)/(int)MAXLINE;
		for(n=0;n<nn;n++) cout << usagecontrol[n] << endl;
	}
	if(mode == 0 || mode == 2) {
		nn = sizeof(usagedesign)/(int)MAXLINE;
		for(n=0;n<nn;n++) cout << usagedesign[n] << endl;
	}
	if(mode == 0 || mode == 3) {
		nn = sizeof(usageseqfile)/(int)MAXLINE;
		for(n=0;n<nn;n++) cout << usageseqfile[n] << endl;
	}
}

CCTRL & CCTRL::operator= (CCTRL &rhs)
 {
	 InitOpts();
	         nDS = rhs.nDS;       
	         nIS = rhs.nIS;       
	    fnDSites = rhs.fnDSites;       
	      fnExec = rhs.fnExec;       
	    fnOutPut = rhs.fnOutPut;       
	fnOutPutDump = rhs.fnOutPutDump;       
	   fnSeqFile = rhs.fnSeqFile;       
	   fnEciFile = rhs.fnEciFile;       
	        reDo = rhs.reDo;       
	       nMode = rhs.nMode;       
	     mHighCF = rhs.mHighCF;       
	      mnNSeq = rhs.mnNSeq;       
	 nAaAtEaDsSt = rhs.nAaAtEaDsSt;  
	        dsSt = rhs.dsSt;  
	      rvDsSt = rhs.rvDsSt;  
	   nSmStOfCl = rhs.nSmStOfCl;  
	        clSz = rhs.clSz;  
	        clIf = rhs.clIf;  
	  saveOption = rhs.saveOption;  
	      strVer = rhs.strVer;  
	   cosChiCut = rhs.cosChiCut;  
	      facDF1 = rhs.facDF1;  
	      facDF2 = rhs.facDF2;  
	    elimMode = rhs.elimMode;  
	    nElimCFs = rhs.nElimCFs;  
	     numSeqs = rhs.numSeqs;  
	       nMScr = rhs.nMScr;  
	 nSrchRemCFs = rhs.nSrchRemCFs;  
	   prDesInfo = rhs.prDesInfo;  
	     prLogEn = rhs.prLogEn;  
	     prLogPn = rhs.prLogPn;  
	    prRemCFs = rhs.prRemCFs;  
	     prCVRMS = rhs.prCVRMS;  
	    prDefCFs = rhs.prDefCFs;  
	     prTrain = rhs.prTrain;  
	  prCEDetail = rhs.prCEDetail;  
	       nData = rhs.nData;  
	    nPreStep = rhs.nPreStep;  
	   nCoolStep = rhs.nCoolStep;  
	       dbl_T = rhs.dbl_T;  
	       codeN = rhs.codeN;  
	    nIntStep = rhs.nIntStep;  
	     dbl_stT = rhs.dbl_stT;  
	     dbl_ndT = rhs.dbl_ndT;  
	    dlb_incT = rhs.dlb_incT;  
	     potFile = rhs.potFile;  
	  fnMCOutPut = rhs.fnMCOutPut;  
	  for(int i=0; i<MAXINPUTITEM;i++) {
	    lstInput[i] = rhs.lstInput[i]; 
	  }
	  return *this;
}

