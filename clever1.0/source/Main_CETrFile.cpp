#include "cemethod.hpp"
#include "argparser.hpp"
#include <string>
EXTERNDEF

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
        "CETrFile trains an energy estimator using the cluster expansion method ",
		"starting from both a design file and a sequence file that includes ",
		"associated energies (or other sequence-related property). To remove   ",
		"parts of sequence space not well described by the CE method, such",
		"sequence groups may be identified and eliminated based on their having",
		"cluster functions with large normalized heterogeneity ratios(NH-ratio).",
		" ",
		"Usage :                                                                      ",
		"   CETrFile  [OPTIONS]                  ",
		"     ex) CETrFile -d design.file -s sequence.file -r training.result",
		"         CETrFile -d design.file -s sequence.file -r training.result",
		"                  -n 3000 -b 5 -l output.file                        ",
		" ",
		"   REQUIRED OPTIONS                                                         ",
		"     -d [design file]         : design file containing design information, ",
		"                                (number of design sites and list of possible",
		"                                amino acids at each site) and structural",
		"                                information (clusters)",
		"     -s [sequence file]       : training set containing sequences and ",
		"                                their associated energies  ",
		" ",
		"   OTHER OPTIONS                                                             ",
        "     -r [training result]     : binary results file containing energy",
		"                                estimator trained using CE method",
		"                                (default : trainresults.eci)",
		"     -l [output file]         : log file containing various information",
		"                                generated during training process. If",
		"                                this option is omitted then the log results",
		"                                print on standard output device",
		"     -n [number of sequences] : number of sequences from training set to ",
		"                                use, starting from beginning of sequence ",
		"                                file (default : number of sequences in", 
		"                                a sequence file) ",
		"     -b [d parameter]         : a paramter to adjust the stringency of CF",
		"                                selection. Increasing the value of ",
		"                                this option lowers the aceptance ratio for",
		"                                selecting new CFs, so that the total number",
		"                                of selected CFs decreases (default 0.0)    ",
		"     -e [number]              : number of sequence groups to be eliminated  ",
		"                                by eliminating sequences containing the ",
		"                                cluster functions with maximum NH-ratio ",
		"                                (default : 0) ",
		"     -m [number]              : minimum number of sequences for a CF. If a ",
		"                                training set does not contain this number",
		"                                of sequences for a given CF, then the CF and ",
		"                                its associated sequences are excluded from ",
		"                                being selected. This value must be 3 or larger",
		"                                (default : 5)", 
		"     -t [mode]                : determine how to deal with higher-order CFs ",
		"                                (triplet) (default : 1)",
		"                                mode=0 use all triplets listed in design file ",
		"                                mode=1 use all triplets listed in design file ",
		"                                and use triplets made from selected pair CFs.",
		"                                (if all three possible pair CFs exist between",
		"                                the three sites, include the triplet CF for",
		"                                those three sites",
		"     -v                       : print version information                  ",
		"     -h                       : display this help and exit                  ",
		"     -h  design               : display an example design file        ",
		"     -h  sequence             : display an example sequence file      "
	};
	int n, nn;
	nn = sizeof(usage)/(int)MAXLINE;
	for(n=0;n<nn;n++) cout << usage[n] << endl;
}

bool AllocOption(CARGPARSER &args, CCTRL &tC) 
{
	CSTRING filename;
	int i;

	if(args.CheckArgs("v")) {
		cout << "CLEVER  " << VERSION << endl; 
		return false;
	}

	if(args.CheckArgs("h")) { // display help
		for(i=1;i<args.GetNCommands();i++) {
			if(args.GetCommand(i) == "design") { tC.Usage(2); exit(-1);	}
			if(args.GetCommand(i) == "sequence") {tC.Usage(3); exit(-1);}
		}
		Usage(); return false;
	}
	if(args.GetParam("d").GetLen() == 0) { // check required options
		cout << "Error: missing design file" << endl;
		Usage(); return false;
	} else {
		tC.fnDSites = args.GetParam("d");
		if(!CheckFile(tC.fnDSites.GetString())) {
			cout << "Error: incorrect design file" << endl;
			Usage(); tC.Usage(2); return false;
		}
	}
	if( args.GetParam("s").GetLen() ==0){
		cout << "Error: missing sequence file" << endl;
		Usage(); return false;
	} else {
		tC.fnSeqFile = args.GetParam("s");
		if(!CheckFile(tC.fnSeqFile.GetString())) {
			cout << "Error: incorrect sequence file" << endl;
			Usage(); tC.Usage(3); return false;
		}
	}
	if(args.GetNCommands() > 1 ) {  // unknown arguments
		cout << "Unknown arguments in command line: " << endl;
		for(i=1;i<args.GetNCommands();i++) printf (" %s",args.GetCommand(i).GetString());
		printf("\n");
		Usage(); return false;
	}
	if(args.GetParam("r").GetLen() !=0) tC.fnEciFile = args.GetParam("r");
	else tC.fnEciFile = "trainresults.eci";

	filename = tC.fnEciFile + "_dump";
	tC.fnOutPutDump = GetSaveFileName(filename);

	if(args.GetParam("n").GetLen() !=0) tC.numSeqs = args.GetParam("n").CvtInt();
	else {
		CUTIL util;
		FILE *fp;
		fp = util.SafeOpen(tC.fnSeqFile.GetString(),"r");
		char buff[1000];
		int count=0;
		while(fgets(buff,1000,fp)) {
			if (strlen(buff)> 2) count++;
			else break;
		}
		fclose(fp);
		cout << "Setting number of training sequences to " << count << endl;
		tC.numSeqs = count;
	}

	if(args.GetParam("b").GetLen() !=0) tC.facDF2 = args.GetParam("b").CvtDbl();
	else tC.facDF2 = 0.0;
	 
	if(args.GetParam("m").GetLen() !=0) {
		i = args.GetParam("m").CvtInt();
		if(i > 2) tC.mnNSeq = i;
		else {
			cout << "Error: number for option \"m\" must be 3 or higher" << endl;
			Usage(); return false;
		}
	} 
	if(args.GetParam("l").GetLen() !=0) tC.fnOutPut = args.GetParam("l");

	if(args.GetParam("t").GetLen() !=0) tC.useHiarc = args.GetParam("t").CvtInt();
	else tC.useHiarc = 1;

	if(args.GetParam("e").GetLen() !=0) tC.nElimCFs = (int)args.GetParam("e").CvtInt();
	else tC.nElimCFs = 0;
	
	return true;
}

int main(int argc, char *argv[])
{
__TIMESTART__
__STARTEXCEPTION__

	CAPPCEM ce;

	// List of Log file redirection
	ofstream *fstdout;
	fstdout = new ofstream;
	streambuf* sbuf1 = cout.rdbuf();
	streambuf* sbuf2 = fstdout->rdbuf();

	// Initialize Options to carry out CE for the data from file
	ce.tC.InitOpts(); // Initialize options
	ce.tC.nMode = 1;
	ce.tC.nSrchRemCFs=0;

	// Option Analysis
	CARGPARSER args;
	args.InitArgs("hv",0);      // no parameter is required
	args.InitArgs("dsrnbltme",1); // one parameter is required
	if(!args.Parse(argc,argv)) {
		ERRMSG(Error in the command line);
		Usage(); return -1;
	}
	if(!AllocOption(args, ce.tC)) return -1;
	if(args.GetParam("l").GetLen() !=0) {  // redirect output information to a log file
		fstdout->open(ce.tC.fnOutPut.GetString());
		cout.rdbuf(sbuf2);
	}
	// Read Design file
	ce.ReadDesign(ce.tC.fnDSites);
	// Print header file
	ce.PrintHeader();
	// print Options
	cout << ce.tC << endl;
	cout << "\n### Log ###" << endl;
	ce.CalcCE_File();
	cout << "\nTotal running time : ";
	APPENDTIME;
	cout << endl;

	ce.PrintResults();
	ce.SaveTrainResults();

__PRINTTIME__
__TIMEEND__

	cout.rdbuf(sbuf1);
	delete fstdout;
	// delete a temporary dump file
	remove(ce.tC.fnOutPutDump.GetString());
	
__ENDEXCEPTION__
	return 0;
}

bool CheckFile(const char * file)
{
	ifstream fin(file,std::ios::in|std::ios::binary);
	if(!fin.is_open()) return false;
	return true;
}

