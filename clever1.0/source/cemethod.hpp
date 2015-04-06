#pragma once
#include "cutil.hpp"
#include "cctrl.hpp"
#include "ccfs.hpp"
#include "ctset.hpp"
#include "mcsimul.hpp"

class CAPPCEM {
public:
	void CalcCE_File();    // carry out CE method from a sequence file
	void CalcCE_Rand();    // 

	void GenRandSeqs();    // generate random sequences for printing to standard output
	void TrainResults();
	void PrintResults();
	void PrintCEFromSeq();
	void PrintTestResults();
	void ReadInput(CSTRING fname);  // read both a control file and design file (fname is control file name)
	void ReadDesign(CSTRING fname); // read a design files (fname means design file name)
	void Save(int np);             // save the data to a file for future evaluation
	void Save();
	bool Read();                    // read a saving file for recalculation
	void PrintECIs(CCFS &);         // print ECI values
	void PrintHeader();             // print header of a log file
	void PrintRandSeqs();           // print te randome sequences to standard output
	void SaveTrainResults();        // save an eci file, this function reduces file size 
	CAPPCEM();
	CCTRL         tC;          // contain control file data
	CTSET       tSet;          // contain training set
	CCFS        tCFs;          // contain cluster functions
	CMCSIMUL     tMC;          // MC simulator/ Random sequence generator
	CUTIL      tUtil;          // General utility
private:
	int       cStep;          // current step
	int       sStep;          // saved step
	int    nUpCFsOk;          // control parameter for managing save position
	int      rmsEci;
	int      countN;
	int    nStep_cs;          // number of step to go out second loop
	int       nStep;          // number of step to go out first loop
	int      itsNum;
	int    itsCount;
	int    nECIStep;
	CSTRING   strVer;          // software version
	CSTRING rFileVer;          // reading file version
};

