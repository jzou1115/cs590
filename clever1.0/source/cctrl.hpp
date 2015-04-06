#pragma once
/***********************************************************************************
	class CCTRL
	The purpose of CCTRL class is to manage user input and to control CE routine.
************************************************************************************/

#include <math.h>
#include <time.h>
#include "general.hpp"
#include "cstring.hpp"
#include "cnvec.hpp"
#include "cnmat.hpp"
#include "cinput.hpp"

class CCTRL {
public:
	int             nDS;       // number of design sites
	int             nIS;       // number of clusters without a reference cluster
	CSTRING     fnDSites;       // name of a file containing design sites/cluster information
	CSTRING       fnExec;       // name of external program for energy calcuation of sequences
	CSTRING     fnOutPut;       // file name to record stdandard output   
	CSTRING fnOutPutDump;       // a file name to record recalculation point
	CSTRING    fnSeqFile;       // sequence file name containing {energy, sequence} information
	CSTRING    fnEciFile;       // a file name for saving ECI train results
	int             reDo;       // 1 : start from a saving position, 0: Calculate from an initial condition 
	int            nMode;       // selection of calculation routine,
	                            // mode 1: extract CE from sequence file, mode 2: extract CE from random generation
	                            // mode 3: calculate CE by eliminating bad sequences from a file
	int          mHighCF;       // value = 1 : consider all Clusters included in a design file
	                            // value = 2 : consider only combination of CFs which contains all lower order CFs
	                            // when value = 2 : if higher order CFs (larger than pair) are described then CFs are calculated in the pool
	                            // of the described CFs
	int           mnNSeq;       // selected CFs should have sequences more than this number of sequences 
	CINPUT       lstInput[MAXINPUTITEM];  // mapping control information from file to variables in CCFS class
	CNVEC_INT    nAaAtEaDsSt;    // number of amino acids at each design site
/*****************************************************************************************
  Design        |    Sequence         External index              Internal index
     0 1 2 3 4  |    1 2 3 ...        1  2  3 ...                 1 2 3 ...
    ----------- |    ---------        -----------       "rvDsSt"  ---------
  1  C A D E G  |    C F R ...        3  6 18 ...      -------->  0 1 2 ...
  2  A F G H K  |    D H Q ...  <==>  4  8 17 ...                 2 3 3 ...
  3  F T R Q K  |    A A T ...        1  1 20 ...       "dsSt"    1 0 1 ...
  ...           |    E F K ...        5  6 11 ...      <--------  3 1 4 ...
  The external index is used to describe sequence information with numbers in an amino-acid space.
  The internal index is used to describe sequence information with numbers in a design space.
  "rvDsSt" function is used to convert the external index to the internal index
     rvDsSt(Design site, External index)
     ex) rvDsSt(0,{3,4,1,5}) ={0,2,1,3}
  "dsSt" function is used to convert the internal index to the external index
     dsSt(Design site, Internal index)
     ex) dsSt(0,{0,2,1,3}) = {3,4,1,5}
  Here, some amino acids which do not belong to the given design space are expressed as zero 
  in the design space.
******************************************************************************************/
	CNMAT_INT         dsSt;    // information of design sites, (all amino acids at each design site
	CNMAT_INT       rvDsSt;    // reverse mapping for dsSt
	CNVEC_INT    nSmStOfCl;    // number of symmetric sites of each cluster [0...nIS-1]
	CNVEC_INT         clSz;    // cluster size  [0...nIS-1]
	CNMAT_INT         clIf;    // cluster information [0...nIS-1, interaction sites and symmetric sites of a cluster]
	int          saveOption;    // 0 : overwrite, 1 : make a saving file for each updating point
	CSTRING          strVer;    // version control for compatibility to read a saving file of old version
	double        cosChiCut;    // bad interaction cutoff criterion (error cosine)
	double           facDF1;    // a weighting factor for degrees of freedom used in screening CFs for point CFs
	double           facDF2;    // a weighting factor for degrees of freedom used in screening CFs for pair and higher order CFs
	int           elimMode;    // method of elimination of CFs, 0: using cutCosChi, 1: using nElimCFs
	int           nElimCFs;    // used when elimMode=1; number of sequence groups to be eliminated which contain bad CFs
	int            numSeqs;    // number of used sequences for training
	int              nMScr;    // mode for screening CFs, mode=0: both point and pair CFs screened
	                            // mode=1: only pair CFs screened
	int        nSrchRemCFs;    // number of searched number of CFs which are not selected for lsSeCf
	// Print options
	int          prDesInfo;    // Print design information to standard output
	int            prLogEn;    // Print energy calculation process
	int            prLogPn;    // Print participation number of CFs which has less sequences than criterion
	int           prRemCFs;    // Print information of some other CFs in lsOtCFs
	int            prCVRMS;    // Print CVRMS data evaluated in the course of calculation
	int           prDefCFs;    // Print definition of all CFs
	int            prTrain;    // Print the training set
	int         prCEDetail;    // Print all information (PR, Gamma, ECI, CFs, error cosign) appeared in the process of calculation
	int           useHiarc;    // Determine whether using hierarchy of CFs or not. 1=> use, 0=> not use
	// MC simulation
	int              nData;    // total number of predicted sequences
	int           nPreStep;    // number of pre-equilibrium step
	int          nCoolStep;    // number of annealing step
	double            dbl_T;    // simulation temperature
	CSTRING           codeN;    // simulation method, {e.g. MC, Simulated annealing, temperature scan..}
	int           nIntStep;    // step interval between data points
	double          dbl_stT;    // starting temperature
	double          dbl_ndT;    // ending temperature
	double         dlb_incT;    // temperature increment for a step 
	CSTRING         potFile;    // Energy calculation file, eci file from ECI trainer
	CSTRING      fnMCOutPut;    // Results file for saving MC results

	CCTRL();
	// for version control to read old version saving file
	void SetVer(CSTRING str){strVer = str; };
	// read a control file
	void ReadCtrl(FILE *fp);
	// read a design file
	void ReadDesign(FILE *fp);
	// print usage of this program
	void Usage(int mode);
	// check whether the file is a control file or a CE-training-results file (eci file) 
	bool CheckECIFile(char *);
	bool CheckECIFile(CSTRING & file);
	// Initialize options
	void InitOpts();
	CCTRL & operator= (CCTRL &rhs);
};   // Control data;

ofstream& operator<< (ofstream& fout, CCTRL &ctrl);
ifstream& operator>> (ifstream& fin, CCTRL &ctrl);
ostream& operator<< (ostream& cout, CCTRL &ctrl);
istream& operator>> (istream& cin, CCTRL &ctrl);
