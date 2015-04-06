#include "cemethod.hpp"
#include "argparser.hpp"
#include <string>
EXTERNDEF

void PrintFileEn(CAPPCEM &ce, CARGPARSER &args);
void PrintFileEnF(CAPPCEM &ce, CARGPARSER &args);
void PrintSeqEn(CAPPCEM &ce, CARGPARSER &args, CSTRING &seq);
bool CheckFile(const char *);

std::string GetSaveFileName(std::string iname) 
{
	char fname[MAXLINE];
	int i=0;
	do {
		sprintf(fname,"%s_%d",iname.c_str(),i++);
	} while(CheckFile(fname));
	std::string value = fname;
	return value;
}

void Usage() 
{
	char usage[][MAXLINE]={
		" ",
		"CEEnergy calculates energies of sequences using an energy estimator",
		"previously trained with the cluster expansion method.              ",
		"The number of amino acids in a provided sequence should be equal to",
		"the number of design sites                                         ",
		" ",
		"Usage :                                                            ",
		"   CEEnergy [OPTIONS]                     ",      
		"    ex) CEEnergy -l output.file -r training.result -s sequence.file", 
		"        CEEnergy -l output.file -r training.result -u \"C D ...\"  ",
		" ",
        "   REQUIRED OPTIONS ",
        "     -r [training result] : results file containing energy       ",
		"                            estimator trained using CE method    ",
		"                            (or trained CE).                     ",
		"     Note: either of the options below, \"s\" or \"u\", must be ",
		"           designated",
		"     -s [sequence file]   : set of sequences for which to evaluate",
		"                            energy, Length of sequence should be ",
		"                            equal to the number of design sites ",
		"                            in trained CE. ",
 		"     -u \"single sequence\" : user can input a single sequence     ",
		"                            instead of a sequence file. Length of",
		"                            the sequence should be equal to the  ",
		"                            number of design sites in trained CE ",
		" ",
        "   OTHER OPTIONS ",
		"     -l [output file]     : file name for standard output        ",
        "     -k                   : this option executes a routine which ",
		"                            checks whether the sequence contains ",
		"                            any of the eliminated CFs. If the    ",
		"                            sequence contains any eliminated     ",
		"                            CF/CFs, then sequence energy is",
		"                            marked as \"NONE\"                   ",
		"     -i [num]             : print energy information ",
		"                            num=0: print only predicted energy   ",
		"                                   (CeEn) (default)              ",
		"                                1: print real energy (ReEn, if in",
		"                                   the sequence file the user",
		"                                   listed the structure-derived or",
		"                                   experimentally-derived energy)",
		"                                   and predicted energy (CeEn)",
		"                                2: print real energy (ReEn),     ",
		"                                   predicted energy (CeEn) energy",
		"                                   difference (EnDiff=ReEn-CeEn) ",
		"                                   , and scaled energy difference",
		"                                   (ScEnDiff)",
		"     -p                   : turn off rapid calculation mode      ",
		"                            (reducing memory requirement)       ",
		"     -c                   : print sequence information ",
		"     -d                   : print design information ",
		"     -h                   : display this help and exit                 ",
		"     -h  sequence         : display an example sequence file     "
	};
	int n, nn;
	nn = sizeof(usage)/(int)MAXLINE;
	for(n=0;n<nn;n++) cout << usage[n] << endl;
}

int AllocOption(CARGPARSER &args, CCTRL &tC) 
{
	CSTRING filename;
	int  i;

	if(args.CheckArgs("h")) { // display help
		for(i=1;i<args.GetNCommands();i++) {
			if(args.GetCommand(i) == "sequence") {tC.Usage(3); return 9;}
		}
		Usage(); return 9;
	}
	if(args.GetNCommands() > 1 ) {  // unknown arguments
		printf("Unknown arguments in command line: ");
		for(i=1;i<args.GetNCommands();i++) printf (" %s",args.GetCommand(i).GetString());
		printf("\n");
		Usage(); return 9;
	}
	if(!args.CheckArgs("r")) {
		printf("Error: the file name for trained CE is required\n");
		Usage(); return 9;
	}
	if(args.CheckArgs("su")) {
		printf("Error: the options \"s\" or \"u\" can not be written simultaneously\n");
		Usage(); return 9;
	}
	if(!(args.CheckArgs("s") || args.CheckArgs("u"))){
		printf("Error: either option of \"s\" or \"u\" is required\n");
		Usage(); return 9;
	}
	if(args.GetParam("s").GetLen() !=0) {  // Check DESIGN file
		tC.fnSeqFile = args.GetParam("s");
		if(!CheckFile(tC.fnSeqFile.GetString())) {
			ERRMSG(Could not open designated sequence file);
			Usage(); tC.Usage(3); return 9;
		}
	}
	if(args.GetParam("l").GetLen() !=0) {  // redirect output information to a log file
		tC.fnOutPut = args.GetParam("l");
	}
	return 0;
}

int main(int argc, char *argv[])
{
__TIMESTART__
__STARTEXCEPTION__
	CAPPCEM ce;
	FILE *fp;

	ofstream *fstdout;
	fstdout = new ofstream;
	streambuf* sbuf1 = cout.rdbuf();
	streambuf* sbuf2 = fstdout->rdbuf();

	// Option Analysis
	CARGPARSER args;
	args.InitArgs("hkdcp",0);       // no parameter is required
	args.InitArgs("slrui",1);       // one parameter is required
	if(!args.Parse(argc,argv)) {
		ERRMSG(Error in the command line);
		Usage(); return -1;
	}
	if (AllocOption(args, ce.tC)) return -1;
	if(args.GetParam("l").GetLen() !=0) {  // redirect output information to a log file
		fstdout->open(ce.tC.fnOutPut.GetString());
		cout.rdbuf(sbuf2);
	}
	// Read information of CE from the training result
	ce.tC.reDo = 1;
	if(ce.tC.CheckECIFile(args.GetParam("r"))) {
		ce.tC.fnOutPutDump = args.GetParam("r");
		if (!ce.Read()){
			ERRMSG(error in the training result file);
			Usage(); return -1;
		}
	} else  return -1;

	// Print energy
	if(args.CheckArgs("d")) {
		ce.PrintHeader();
		cout << ce.tC << endl;
		cout << "\n### Estimated Energy ###" << endl;
	}
	if (args.CheckArgs("u")) {     // a command line contains sequence information
		CSTRING seq = args.GetParam("u");
		int n = seq.SplitLine(" ");
		if( n != ce.tC.nDS) { // a sequence in the command line mismatch the design sites
			ERRMSG(Mismatch in sequence information);
			return -1;
		}
		PrintSeqEn(ce, args, seq);
	}else {              // a file contains sequence information
		fp = ce.tUtil.SafeOpen(ce.tC.fnSeqFile.GetString(),"r");
		ce.tSet.ReadTrainingSet(fp,ce.tC);
		fclose(fp);
		if(args.CheckArgs("p")) {
			PrintFileEn(ce, args);
		} else {
			PrintFileEnF(ce, args);
		}
	}
__DEBUGTIMEEND__
	if(args.GetParam("l").GetLen() !=0) {
		cout.rdbuf(sbuf1);
	}
	delete fstdout;
__ENDEXCEPTION__
	return 0;
}

void PrintFileEnF(CAPPCEM &ce, CARGPARSER &args) {
	CNVEC_INT tmSeq;
	CNVEC_INT tmNu;
	int nn;
	double r_en;
	char buff[MAXLINE], entxt[20];
	int imode=0;
	if(args.CheckArgs("i")) { imode = args.GetParam("i").CvtInt();}
	double cvrms = ce.tCFs.crVlSmr.GetEndVal();
	if(cvrms == 0.0 && imode == 2) imode = 1;

	switch (imode) {
		case 0:
			sprintf(buff,"%12s", "CeEn");
			cout << buff << endl;
			break;
		case 1:
			sprintf(buff,"%12s %12s", "ReEn", "CeEn");
			cout << buff << endl;
			break;
		case 2:
			sprintf(buff,"%12s %12s %12s %10s  ", "ReEn", "CeEn", "EnDiff", "ScEnDiff");
			cout << buff << endl;
			break;
	}
	// energy evaluation
	CNMAT_INT scdmat;
	CNVEC_DBL envec;
	ce.tCFs.CalcSCD(ce.tSet.dataSet,scdmat);
	envec = scdmat*ce.tCFs.lsEci;

	//print results
	nn = ce.tSet.GetNSeq(D_DATA);
	bool seq = args.CheckArgs('c');
	bool quick = !args.CheckArgs('k');
	bool ien;
	for(int i=0;i<nn;i++) {
		ce.tSet.GetSeq(D_DATA,i,tmSeq,r_en);
		ien = false;
		if(quick) {
			sprintf(entxt,"%12.4f",envec[i]);
		} else {
			if(ce.tCFs.IsItForbidden(tmSeq)){
				sprintf(entxt,"%12s","None");
				ien = true;
			} else {
				sprintf(entxt,"%12.4f",envec[i]);
			}
		}
		switch (imode) {
			case 0:
				sprintf(buff,"%12s", entxt);
				break;
			case 1:
				sprintf(buff,"%12.4f %12s", r_en, entxt);
				break;
			case 2:
				if(ien) {
					sprintf(buff,"%12.4f %12s %12s %10s", r_en, entxt,"None","None");
				}else{
					tmNu = scdmat.GetVec(ROW,i);
					double gammahetero = tmNu*ce.tCFs.htEciCo*tmNu;
					double heterosig = sqrt(cvrms*cvrms + gammahetero);
					double diff = r_en - envec[i];
					sprintf(buff,"%12.4f %12s %12.4f %10.6f", r_en, entxt, diff, diff/heterosig);
				}
				break;
		}
		if(seq) cout << buff <<  "  " <<ce.tSet.GetSeqText(ce.tC,tmSeq) << endl; 
		else cout << buff << endl;
	}
}

void PrintFileEn(CAPPCEM &ce, CARGPARSER &args) {
	CNVEC_INT tmSeq;
	CNVEC_INT tmNu;
	int nn;
	double en, r_en;
	char buff[MAXLINE], entxt[20];
	int imode=0;
	if(args.CheckArgs("i")) { imode = args.GetParam("i").CvtInt();}
	double cvrms = ce.tCFs.crVlSmr.GetEndVal();
	if(cvrms == 0.0 && imode == 2) imode = 1;

	switch (imode) {
		case 0:
			sprintf(buff,"%12s", "CeEn");
			cout << buff << endl;
			break;
		case 1:
			sprintf(buff,"%12s %12s", "ReEn", "CeEn");
			cout << buff << endl;
			break;
		case 2:
			sprintf(buff,"%12s %12s %12s %10s  ", "ReEn", "CeEn", "EnDiff", "ScEnDiff");
			cout << buff << endl;
			break;
	}
	nn = ce.tSet.GetNSeq(D_DATA);
	bool seq = args.CheckArgs('c');
	bool quick = !args.CheckArgs('k');
	bool ien;
	for(int i=0;i<nn;i++) {
		ce.tSet.GetSeq(D_DATA,i,tmSeq,r_en);
		ien = false;
		if(quick) {
			en = ce.tCFs.CalcCEEnergy(tmSeq);
			sprintf(entxt,"%12.4f",en);
		} else {
			if(ce.tCFs.IsItForbidden(tmSeq)){
				sprintf(entxt,"%12s","None");
				ien = true;
			} else {
				en = ce.tCFs.CalcCEEnergy(tmSeq);
				sprintf(entxt,"%12.4f",en);
			}
		}
		switch (imode) {
			case 0:
				sprintf(buff,"%12s", entxt);
				break;
			case 1:
				sprintf(buff,"%12.4f %12s", ce.tSet.dataEn[i], entxt);
				break;
			case 2:
				if(ien) 
					sprintf(buff,"%12.4f %12s %12s %10s", r_en, entxt,"None","None");
				else{
					ce.tCFs.CalcSCD(tmSeq, tmNu);
					double gammahetero = tmNu*ce.tCFs.htEciCo*tmNu;
					double heterosig = sqrt(cvrms*cvrms + gammahetero);
					double diff = r_en - en;
					sprintf(buff,"%12.4f %12s %12.4f %10.6f", r_en, entxt, diff, diff/heterosig);
				}
				break;
		}
		if(seq) cout << buff <<  "  " <<ce.tSet.GetSeqText(ce.tC,tmSeq) << endl; 
		else cout << buff << endl;
	}
}

void PrintSeqEn(CAPPCEM &ce, CARGPARSER &args, CSTRING &seq) {
	// Interpret sequence data
	CNVEC_INT tmSeq;
	int nn;
	double en;
	char buff[MAXLINE], entxt[20];
	int imode=0;
	if(args.CheckArgs("i")) { imode = args.GetParam("i").CvtInt();}

	tmSeq.SetDim(ce.tC.nDS);
	for(nn=0;nn<ce.tC.nDS;nn++) {
		tmSeq[nn] = ce.tC.rvDsSt(nn,seq(nn)[0] - 'A' + 1);
		//rvDsSt convert sequence information to internal format for energy evalution
	}
	bool bseq = args.CheckArgs('c');
	bool quick = !args.CheckArgs("k");

	if(quick) {
		en = ce.tCFs.CalcCEEnergy(tmSeq);
		sprintf(entxt,"%12.4f",en);
	} else {
		if(ce.tCFs.IsItForbidden(tmSeq)){
			sprintf(entxt,"%12s","None");
		} else {
			en = ce.tCFs.CalcCEEnergy(tmSeq);
			sprintf(entxt,"%12.4f",en);
		}
	}
	sprintf(buff,"%12s", entxt);
	if(bseq) cout << buff <<  "  " <<ce.tSet.GetSeqText(ce.tC,tmSeq) << endl; 
	else cout << buff << endl;
}

bool CheckFile(const char * file)
{
	bool ok;
	ifstream f_tmp_in(file);
	ok = f_tmp_in.is_open();
	return ok;
}

