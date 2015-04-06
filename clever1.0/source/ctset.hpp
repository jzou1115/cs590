#pragma once

#include <math.h>
#include <time.h>
#include "general.hpp"
#include "cstring.hpp"
#include "cctrl.hpp"
#include "ccfs.hpp"
#include "cnvec.hpp"
#include "cnmat.hpp"

class CCTRL;
class CCFS;

enum DATA_SET { D_ITS, D_ADDED, D_ELIM, D_DATA};
class CTSET {
public:
	CNMAT_INT itsSet;      // a training set
	CNVEC_DBL  itsEn;      // energy of a training set
	CNVEC_INT itsOrd;      // order of sequences
	CNMAT_INT addedSet;    // the sequences randomly generated are allocated to this variable 
	CNVEC_DBL  addedEn;    // energy of additional set
	CNVEC_INT addedOrd;    // order of sequences
	CNMAT_INT elimSet;     // eliminated sequences are saved to this variable
	CNVEC_DBL  elimEn;     // eliminated energy set
	CNMAT_INT dataSet;     // sequences from a file are saved to this variable
	CNVEC_DBL  dataEn;     // energies of dataSet
	CNMAT_INT sel_BadSeqs; // worst-fitting sequences for worst-fitting CFs
	int          aN;       // this variable used for efficient converting to probe matrix 
	int  nAdd;             // number of added sequences for a loop
	CSTRING strVer;        // used for compatibility to read old version saving file

	CTSET();
	// compatibility for reading a saved file from old version
	void SetVer(CSTRING str){strVer = str; };
	// Read a training set from a sequence file
	void ReadTrainingSet(FILE *fp, CCTRL &tC);
	// Energy evalution by using a program which is supplied by user
	void CalcEnergy(CCTRL &, CCFS &);
	// This routine is used to update energies after increasing backbone space
	void ReCalcEnergy(CCTRL & tC);
	// add addedSet data to itsSet
	void Union();
	// print a training set information
	void PrintSeqs(CCTRL &tC, CCFS &tCFs);
	// Print sequences in the "addedSet" variable
	void PrintAddSeqs(CCTRL &tC);
	// append a sequence which is not in the training set to the addedSet
	//bool AddSeq(CNMAT_INT &addV);
	void AddSeq(DATA_SET set, CNVEC_INT &addV);
	void AddSeq(DATA_SET set, CNVEC_INT &addV, double en);
	void CopySeq(DATA_SET orig, DATA_SET dest);
	bool CheckSeq(CNVEC_INT &addV);
	// FindSeqPos : If there is same sequence then return -1, 
	// or otherwise return offset position for inserting
	int  FindSeqPos(CNMAT_INT &seqM, CNVEC_INT &ordV, CNVEC_INT &seq);
	// Check whether the given sequence is previously appeared
	// return true when sequence addition is complete
	// return false when sequence addition is fail
	bool CheckAddSeq(DATA_SET set, CNVEC_INT &seq);
	// Update sequence order vector
	void UpSeqOrdIdx(CNMAT_INT &seqM, CNVEC_INT &ordV);
	// Remove sequences
	bool RemSeq(DATA_SET set, int offset);
	void RemSeq(DATA_SET set, CNVEC_INT & list);
	void ClearSeq(DATA_SET set);
	// Move sequences
	void MoveSeq(DATA_SET from_set, CNVEC_INT &list, DATA_SET to_set);
	// Get a sequence vector and its energy from a data set, the sequence vector is saved in "seq"
	bool GetSeq(DATA_SET set, int offset, CNVEC_INT &seq, double &en);
	CNVEC_INT & GetSeq(DATA_SET set, int offset);
	int GetNSeq(DATA_SET set);
	// Get a sequence line from a int vector information
	CSTRING GetSeqText(CCTRL &tC, CNVEC_INT &vec);
	// command line executor
	double Execution( const CSTRING & command);
};
ofstream& operator<< (ofstream& fout, CTSET &ctset);
ifstream& operator>> (ifstream& fin, CTSET &ctset);

