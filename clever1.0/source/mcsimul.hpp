#pragma once
#include "cutil.hpp"
#include "cctrl.hpp"
#include "ccfs.hpp"
#include "ctset.hpp"

class CMCSIMUL
{
public:
	CMCSIMUL(void);
	~CMCSIMUL(void);
	// read a control file for MC
	void ReadInput(char *filename, CCTRL &tC);
	// Print usage of this program
	void Usage(int mode);
	// generate sequences
	void GenTSet(CCTRL &tC, CCFS &tCFs, CTSET &tSet);
	// Change element "num" times based on old sequence "vec"
	void ChangeElements(CCFS &tCFs, CNVEC_INT &vec, int num);
	// Obtain three sequences from random, quenched and middle temperature from MC simulation 
	void MCLowMedHigh(CCFS & tCFs, CTSET &tSet);
	// Obtain sequences from temperature scan
	void MCTempScan(CCFS &tCFs, CTSET &tSet, CNVEC_DBL &tList, CNMAT_INT &listM);
	// MC simulation
	void MCSimul(CCFS &tCFs, CTSET &tSet);
	// Updating parameters from tC
	void UpdateParams(CCTRL &tC);
	// Quenching routine
	void MCAnneal(CCFS & tCFs, CNVEC_INT &initseq, CNVEC_INT &finalseq);
	void SimAnneal(CCFS &tCFs, CTSET &tSet);
	// Obtain a sequence which has an energy close to "en"
	void MCSpecEn(CCFS &tCFs, double en, CNVEC_INT & seq);
	void GenRandSeqs(CCFS &tCFs, CTSET &tSet);
	double Rand2();
	// Random sequence generation routine
	void RandSeq(CCFS &tCFs, CNVEC_INT &seq);

public:
	CUTIL     tUtil;        // utility for general purpose
	int       nData;        // total number of step
	CSTRING fnOutPut;        // name of output file
	CSTRING  potFile;        // CE force field file
	int    nPreStep;        // number of pre-equilibrium step
	int   nCoolStep;        // number of annealing step
	double     dbl_T;        // simulation temperature
	CSTRING    codeN;        // simulation code
	int    nIntStep;        // step interval between data points
	double   dbl_stT;        // starting temperature
	double   dbl_ndT;        // ending temperature
	double  dlb_incT;        // temperature interval 
	int        nAdd;
	bool          bQ;        // skip invalid sequences
	bool        bSeq;        // skip duplicate sequences 

	// Variables to read CE energy
	int              nDS;   // Number of design sites
	CNVEC_INT nAaAtEaDsSt;   // Number of amino acids at each site
};

