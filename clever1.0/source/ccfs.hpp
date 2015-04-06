#pragma once
/***********************************************************************************
	class CCFS
	The purpose of CCFS class is to manipulate and evaluate cluster functions and
	their effective cluster interactions. 
	CCFS class contains the definition and evaluation of cluster functions.
************************************************************************************/

#include <math.h>
#include <time.h>
#include "general.hpp"
#include "cstring.hpp"
#include "cnvec.hpp"
#include "cnmat.hpp"
#include "cctrl.hpp"
#include "ctset.hpp"
#include "ccflist.hpp"

class CCTRL;
class CTSET;

class CCFS {
public:
	/********************************************************************************
	 Preparation step of cluster expansion
	     1 A D C E F G   <-- definition of design sites
		 2 D E F G           nAaAtEaDsStRn[0] = 5; {'A' is a reference}
		 3 A C D E F G
		 4 D E F G
		 -------
		 1; 3            <-- definition of clusters
		 2; 4                clIfRn(0,:)={0, 2}; nSmStOfCl[0]=2
		 1 2; 3 4            clSz[2] = 2;
		 -------
		 0  constant     <-- definition of CFs
		 1  D 1              here 11, 12 mean index of CFs, {D E, D F} decoration
		 8  G 2              
		 9  D E 1 2          pairWtPnt(1,:)={9, 10, 11}
		10  D F 1 2         
		 -------    
		 nTtCf = 24; clRgInCf[0,1,2,3]={0,1,9,24}
    ********************************************************************************/
	CCFLIST         cfDefs;       // Definition of CFs
	CNVEC_INT    nAaAtEaDsStRn;  // number of amino acdis at each design site (renumbered, -1 shifted)
	CNMAT_INT      clIfRn;       // cluster (design sites, zero basis) information
	CNVEC_INT   nSmStOfCl;       // number of symmetric sites in each cluster
	CNVEC_INT        clSz;       // cluster size: point(1), pair(2), triplet (3)
	int             nTtCf;       // number of total CFs
	CNVEC_INT    clRgInCf;       // cluster range in index of CFs
	CNMAT_INT    pmFoAlCf;       // (Future use) matrix variable to save a probe matrix for all cluster functions

	/* parameters used in algorithms to check whether enough sequences are in a training set*/
	int            mnNSeq;       // A training set contains CFs more than mnNSeq times

	/********************************************************************************
	 Evaluation step of cluster expansion
		 If we assume a list of selected CFs as lsSeCf={0,1,8,9} then
		 information of training sequences is converted to information of CFs, 
		    E  Seq               Ep probe matrix (M)          
		 -123  C D A G         -125 1 1 1 0        sePm(0,:)={1,1,1,0}
		 -235  A E D F     =>  -234 1 0 0 0
		 -345  C G C E         -345 1 2 1 1
		 Ep = [M*(MT*M)^-1*MT]*E = [M*eciCo*MT]*E = Pj*E;
		 dfInEn = E-Ep;
		 htEciCo = eciCo*MT*[diagonal matrix of n*(CVRMS vector)^2_i]*M*eciCo;
    ********************************************************************************/
	CNMAT_INT        sePm;       // selected probe matrix (the element is a return value of CF operation on a sequence) 
	CNVEC_DBL        lsEci;       // effective cluster interaction energies
	CNMAT_DBL        eciCo;       // eciCo = (MT*M)^-1, scaled covariance matrix of ECI 
	CNVEC_DBL         dgPj;       // dgPj = [M*(MT*M)^-1*MT]_ii, diagonal part of a projection matrix
	CNVEC_DBL       dfInEn;       // dfInEn = E - E(expected), difference between energies
	CNMAT_DBL      htEciCo;       // heterogeneous covariance matrix of ECI
	double          magErr;       // magnitude of error vector |n*Err_unit^2 -1|= |Dcv|
	CNVEC_DBL         pnVl;       // (Future use) participation number value

	/* variables which are used in cluster function manipulation */
	CNVEC_INT      lsTtCf;       // list of total CFs   (total possible CFs >= lsSeCF+lsOtCf+lsElCf)
	CNVEC_INT      lsSeCf;       // list of selected cluster functions
	CNVEC_INT      lsOtCf;       // temperary eliminated from lsTtCf ( used when updating CFs)
	CNVEC_INT      lsElCf;       // a list of Eliminated CFs (if sequence decreased then not used for updating CFs
	                              // if sequence increased then  used for updating)
	CNVEC_INT elimCFsList;       // a list of eliminated CFs, sequences are also eliminated (used in determining forbidden sequences)
	                              // CFs connected to the elements of elimCFsList is permanently removed in the possible CF list

	/* control parameters which are used in fitting algorithm */
	int             nMScr;       // choose a method for screening CFs (screen point and pair or screen pair only)
	double          facDF1;       // a weighting factor for screening point cluster functions
	double          facDF2;       // a weighting factor for screening higher-order cluster functions
	int       nSrchRemCFs;       // Number of CFs for looking up in the lsOtCf list for elimination (for efficiency of code)
	CSTRING         strVer;       // version control to read old saved file
	int          useHiarc;       // determine whether using hierarchy of CFs

	/* Variables for saving results of cluster expansion*/
	int            nUpCFs;       // count number of UpdateCFs calling
	CNVEC_DBL         crVl;       // CVRMS values
	CNVEC_DBL      crVlSmr;       // collection of final CVRMS value
	CNVEC_DBL       htNrEr;       // normalized error by applying heterogenous error distribution
	CNMAT_DBL         svRs;       // IPR, Gamma, Cos, ECI of each step are saved
	CNMAT_INT     svRsIdx;       // Index of CFs corresponding to row of "svRs"
	CNMAT_DBL      svRsSmr;       // summary of results of each step

    /* determine quality of prediction (future use) */
	CNVEC_INT      ctgCFSeqs;    // Category of CFs of sequences
	CNMAT_DBL       pctSeqsCFs;   // Percentage of sequences containing given set of CFs
	CNVEC_INT      ctgCFList;    // Indicate where the sequences is involved in
	CNMAT_DBL       ctgCFRange;   // Save training-set rmsd value and sequence range

	int mHighCF;                 // Including higher-order CFs
	int maxClOrd;                // Maximum order of CF

	CCFS();
	//Calculate change in CVRMS by comparing old CVRMS value to current value
	//this function was used to determine to finish updating ECI value when including more sequences which 
	//are predicted to have lower eneriges.
	double       CalcChangeInCVRMS();
	// This function is used to version control
	// To read a backup file generated from an old version program
	void         SetVer(CSTRING str){strVer = str; };
	// Calculate difference in energy between authentic energy and CE-derived energy
	void         CalcEnDiff(CNVEC_DBL &env, CNMAT_INT &seqm);
	// Calculate a vector of probe matrix from sequence matrix for the CF "ncf"
	void         CalcSCD(CTSET &tSet, int ncf, CNVEC_INT &cfvec); 
	void         CalcSCD(CNMAT_INT &setSeq, CNMAT_INT &setSCD); 
	// Calculate a probe matrix from a training set for the vector list "vec"
	void         CalcSCD(CNVEC_INT &vec, CNVEC_INT &cfvec);
	void         CalcSCD(CNVEC_INT &seq, CNVEC_INT &cfvec, CNVEC_INT &listcfs); 
	// Initialize CFs
	void         GenInitCFs(CCTRL &);
	// Determine order of for screening point CFs
	// return an ordered list of point CFs
	void DtOrdOfCFsLast(CTSET &tSet, int dimcf, CNVEC_INT & ord);
	// Determine order for screeening pair CFs
	// return an ordered list of pair CFs
	void DtOrdOfCFs1st(CTSET &tSet, int dimcf, CNVEC_INT & ord);
	// Effective update of CFs and ECIs when sequences are increased
	void         UpdateECIs(CCTRL &, CTSET &);
	// Update CFs and ECI values from a training set
	int         UpdateCFs(CCTRL &, CTSET &);
	// Calculate ECI values from a list of CFs designated by "listCFs"
	void         UpdateCFs(CCTRL &tC, CTSET &tSet, CNVEC_INT listCFs);
	// Calculate CE-derived energy from a sequence "seq"
	double       CalcCEEnergy(CNVEC_INT & seq);
	// Select worst CF and return true when the elimination of CF-including sequences are valid or else false 
	bool         ChooseErrInteract(CCTRL & tC, CTSET &tSet);
	// Eliminate sequences which contain the worst CF and choose a worst-fitting sequence
	void         ElimSeqs(CCTRL &tC, CTSET &tSet);
	// return the sequence contains a CF in elimCFsList when modifying "nSite" to the amino acid "nSeq" in the sequence  
	bool         IsItForbidden(int nSite, int nSeq, CNVEC_INT &seqV);
	// return the sequence contains a CF in elimCFsList
	bool         IsItForbidden(CNVEC_INT &seqV);
	// Calculate "htEciCo" 
	int          CalcVarECIs();
	// Calculate "htNrEr"
	void         CalcStdErrData();
	// Calculate ECI values, participation number, and so on for the CFs in the list "lsOtCf"
	void         CalcRemCFsECIs(CCTRL &tC, CTSET &tSet);
	// Categorize sequences into several groups according to energy prediction accuracy
	void         AnalSeqErr(CCTRL &tC, CTSET &tSet);
	// Move a CF "ncfs" and their higher-order-derivated CFs in the list "lsTtCf" to "lsSeCf", "lsOtCf" and "lsElCf"
	// mode 0: move to lsSeCf
	// mode 1: move to lsOtCf
	// mode 2: move to lsElCf
	// mode 3: eliminate from the list
	void         MoveCFs(int ncfs, int mode);
	// point CFs which are appeared in training sequences less than "mnNSeq" times are moved to "lsElCf"
	void         SequenceFilter(CTSET &tSet);
};   // Control data;

ofstream& operator<< (ofstream& fout, CCFS &ccfs);
ifstream& operator>> (ifstream& fin, CCFS &ccfs);
