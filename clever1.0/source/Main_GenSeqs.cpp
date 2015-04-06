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
		"GenSeqs uses a given design file to randomly generate a set of sequences ",
		"by independently drawing residues from the list of allowed amino acids at",
		"each design site. This can be used to either generate a training set or a",
		"test set. If generating a test set, a useful option is -o which allows the",
		"user to generate a test set that does not share sequences with the training",
		"set. Note, the same seed is used for each random generation of sequences.",
		" ",
		"Usage :                                                                    ",
		"   GenSeqs [OPTIONS]                     ",
		"     ex) GenSeqs -n 3000 -d design.file -s sequence.file                   ",
		"         GenSeqs -n 3000 -d design.file -s sequence.file -o sequence.old   ",
		" ",
		"   REQUIRED OPTIONS                                                        ",
		"     -n [number]      : number of sequences to be randomly generated",
		"     -d [design file] : design file containing design information",
		"                        (number of design sites and list of possible",
		"                        amino acids at each site) and structural",
		"                        information (clusters)",
		" ",
		"   OTHER OPTIONS                                                            ",
		"     -o [old sequence]: old sequence file, additional sequences different   ",
		"                        from these old sequences are generated and can be   ",
		"                        used as a test set containing no training set       ",
		"                        sequences",
		"     -s [new sequence]: result file containing additional sequences. If this",
		"                        option is omitted then the generated sequences     ",
		"                        print on standard output device",
		"     -h               : display this help and exit                  ",
		"     -h design        : display an example for a design file        "
	};
	int n, nn;
	nn = sizeof(usage)/(int)MAXLINE;
	for(n=0;n<nn;n++) cout << usage[n] << endl;
}

bool AllocOption(CARGPARSER &args, CCTRL &tC, CMCSIMUL &tMC) 
{
	CSTRING filename;
	int i;

	if(args.CheckArgs("h")) { // display help
		for(i=1;i<args.GetNCommands();i++) {
			if(args.GetCommand(i) == "design") { tC.Usage(2); return false;	}
		}
		Usage(); return false;
	}
	if(args.GetParam("n").GetLen()==0) { // check required options
		cout << "Error: option for number of sequences is required" << endl;
		Usage(); return false;
	} else
		tMC.nData = args.GetParam("n").CvtInt();

	if(args.GetParam("d").GetLen()!=0)
		tC.fnDSites = args.GetParam("d");
	else {
		cout << "Error: design file is required" << endl;
		Usage(); return false;
	}

	if( args.GetParam("o").GetLen() !=0){
		tC.fnSeqFile = args.GetParam("o");
		if(!CheckFile(tC.fnSeqFile.GetString())) {
			cout << "Error: missing old sequence file" << endl;
			Usage(); tC.Usage(3); return false;
		}
	}

	if(args.GetParam("s").GetLen() !=0) {  // Check DESIGN file
		tC.fnOutPut = args.GetParam("s");
	} 

	return true;
}

void Print(CCTRL &tC, CTSET &tSet) {
	int nn, idx;
	char buff[MAXLINE];
	double en;

	for(nn=0;nn<tSet.GetNSeq(D_ITS);nn++) {
		idx = nn;
		CNVEC_INT & tmpV = tSet.GetSeq(D_ITS,idx);
		if(idx < tSet.itsEn.GetItsDim()) en = tSet.itsEn[idx];
		else en=0.0;
		sprintf(buff,"%12.5f %s",en,tSet.GetSeqText(tC,tmpV).GetString());
		cout << buff << endl;
	}
}

int main(int argc, char *argv[])
{
	ofstream *fstdout;
	fstdout = new ofstream;
	streambuf* sbuf1 = cout.rdbuf();
	streambuf* sbuf2;

__TIMESTART__
	CAPPCEM ce;
	FILE *fp;


	// Option Analysis
	CARGPARSER args;
	args.InitArgs("h",0);           // no parameter is required
	args.InitArgs("ndos",1);       // one parameter is required
	if(!args.Parse(argc,argv)) {
		cout << "Error: command line error" << endl;
		Usage(); return -1;
	}
	if(!AllocOption(args, ce.tC, ce.tMC)) {
		return -1;
	}
	ce.tMC.bSeq = true; // true: skip randomly generated sequence appeared in old sequences

	if(args.GetParam("s").GetLen() != 0) {  // redirect output information to a log file
		fstdout->open(ce.tC.fnOutPut.GetString());
		sbuf2 = fstdout->rdbuf();
		cout.rdbuf(sbuf2);
	}

	// Get desgin information
	ce.ReadDesign(ce.tC.fnDSites);
//	ce.tCFs.GenInitCFs(ce.tC);

	if (args.GetParam("o").GetLen()!=0) {
		fp = ce.tUtil.SafeOpen(ce.tC.fnSeqFile.GetString(),"r");
		ce.tSet.ReadTrainingSet(fp,ce.tC);
		fclose(fp);
		ce.tSet.itsEn  = ce.tSet.dataEn;
		ce.tSet.itsSet = ce.tSet.dataSet;
		ce.tSet.UpSeqOrdIdx(ce.tSet.itsSet, ce.tSet.itsOrd);
	}

	int oldnn = ce.tSet.GetNSeq(D_ITS);
	ce.tMC.GenRandSeqs(ce.tCFs,ce.tSet);
	if (oldnn != 0) {
		ce.tSet.itsSet.Remove(COLUMN,0,oldnn-1);
	}
	Print(ce.tC, ce.tSet);

	if(args.GetParam("s").GetLen() != 0) {
		cout.rdbuf(sbuf1);
	}
	delete fstdout;
	return 0;
}

bool CheckFile(const char * file)
{
	ifstream fin(file,std::ios::in|std::ios::binary);
	if(!fin.is_open()) return false;
	return true;
}

