/*******************************************************************************************
	Manage a list of CF definitions

	Example)
	   CCFLIST cfls;
	   CCFDEF  cf;
	   cf.nOrd = 2;  cf.nSym = 1;
	   cfls.Add(cf);
	   cf.nOrd = 2; cf.nSym = 1;
	   cfls.Add(cf);
	   cfls.Print(ce.tC);
	   cfls.Remove(1);
	   cfls.Print(ce.tC);
********************************************************************************************/

#pragma once
#include "general.hpp"
#include "cctrl.hpp"
#include "cstring.hpp"
#include <map>
#include "general.hpp"
#include "ccfdef.hpp"
#include "cnvec.hpp"

class CCFLIST
{
public:
	CCFLIST(void);
	~CCFLIST(void);
	void Add(CCFDEF & cf);
	void RemCFLst(int offset);
	int Length();
	void PrintDefCFs(CCTRL &tC);
	void Clear();
	void SetStIdxMap(int decor, int clust, int ncf); // clust : external index
	void SetIntClIdx(int extidx, int intidx);
	int GetCfIdx(int decor, int clust);              // clust : external index
	int GetCfIdx(int site);
	int GetStIdx(int decor, int clust);              // clust : external index
	CCFDEF  & operator[] (int ncf);
	CCFLIST & operator= (CCFLIST &cflst);
public:
	CNVEC_INT stIdx;        // index of sites for point CFs, stIdx[site] = index of a corresponding point CF
	CNVEC_INT extToInt;     // convert external index for cluster to internal index, clToIdx[ex idx] = int idx
	CNVEC_INT intToExt;     // convert internal index for cluster to external index, idxToCl[int idx] = ex idx
	std::map<int,CCFDEF*> cfLst;
};

ofstream & operator<< (ofstream& fout, CCFLIST &rhs);
ifstream & operator>> (ifstream&  fin, CCFLIST &rhs);

