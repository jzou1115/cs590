#pragma once
#include "general.hpp"
#include "cctrl.hpp"
#include "cstring.hpp"
#include "cnvec.hpp"

class CCFDEF
{
public:
	CCFDEF(void);
	CCFDEF(CCFDEF &cf);
	~CCFDEF(void);
	CSTRING GetCFText(CCTRL &tC);
	CSTRING GetCFDetails(CCTRL &tC);
	CSTRING GetCFInfo(CCTRL &tC, int mode);
	void PutClust(int clust);
	void SetDecor(int decor);
	void PutDecor(int decor);
	void Clear();
	CCFDEF & operator= (CCFDEF &cf);

public:
	CNVEC_INT lClust;      // cluster of CF
	CNVEC_INT lDecor;      // decoration of CF
	CNVEC_INT lHoAddSt;    // list of additional point CF to build higher-order CF based on this CF
	CNVEC_INT lHoCf;       // list of higher-order CFs which include this CF
	CNVEC_INT lHoSmStCf;   // Symmetric site for this CF 
	CNVEC_INT lLoCf;       // list of lower-order CFs in this CF
	CNVEC_INT lSt;         // list of sites for point CFs involved in this CF
	int nSym;              // a number of symmetric sites in this CF
	int nOrd;              // CF size
	int cfIdx;             // Cluster function label;
};

ofstream & operator<< (ofstream& fout, CCFDEF &rhs);
ifstream & operator>> (ifstream&  fin, CCFDEF &rhs);

