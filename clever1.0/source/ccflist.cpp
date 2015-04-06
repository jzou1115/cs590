#include "ccflist.hpp"

ofstream & operator<< (ofstream& fout, CCFLIST &rhs)
{
	CSTRING str;
	str = "CFLIST";
	fout << str;
	int nmax = rhs.Length();
	WRITEDATA(fout,nmax);
	std::map<int, CCFDEF*>::const_iterator iter;
    for (iter=rhs.cfLst.begin(); iter != rhs.cfLst.end(); ++iter) {
		fout << *(iter->second);
	}
	fout << rhs.stIdx;
	fout << rhs.extToInt;
	fout << rhs.intToExt;
	return fout;
}
ifstream & operator>> (ifstream&  fin, CCFLIST &rhs)
{
	CSTRING str;
	fin >> str;
	if(str != "CFLIST") {__ERRORTHROW(ERR_DATAFORM);}
	int nmax;
	READDATA(fin,nmax);
	CCFDEF cf;
	rhs.Clear();
	for(int n=0;n<nmax;n++) {
		fin >> cf;
		rhs.Add(cf);
	}
	fin >> rhs.stIdx;
	fin >> rhs.extToInt;
	fin >> rhs.intToExt;
	return fin;
}

CCFLIST::CCFLIST(void)
{
}

CCFLIST::~CCFLIST(void)
{
	Clear();
}

void CCFLIST::Clear()
{
	std::map<int, CCFDEF*>::const_iterator iter;
    for (iter=cfLst.begin(); iter != cfLst.end(); ++iter) {
		delete iter->second;
	}
	cfLst.clear();
	extToInt.SetDim(0);
	intToExt.SetDim(0);
	stIdx.SetDim(0);
}
void CCFLIST::SetStIdxMap(int decor, int clust, int ncf)
{
	const int MAXDECOR = 30;
	int site;
	site = extToInt[clust]*MAXDECOR + decor;
	if(site >= stIdx.GetItsDim()) stIdx.SetDim(site+1);
	stIdx[site] = ncf;
}
void CCFLIST::SetIntClIdx(int extidx, int intidx)
{
	if(extToInt.GetItsDim()<=extidx) extToInt.SetDim(extidx+1);
	if(intToExt.GetItsDim()<=intidx) intToExt.SetDim(intidx+1);
	extToInt[extidx]=intidx;
	intToExt[intidx]=extidx;
}

int CCFLIST::GetCfIdx(int decor, int clust)
{
	const int MAXDECOR = 30;
	int site;
	site = extToInt[clust]*MAXDECOR + decor;
	return stIdx[site];
}

int CCFLIST::GetCfIdx(int site)
{
	return stIdx[site];
}

int CCFLIST::GetStIdx(int decor, int clust)
{
	const int MAXDECOR = 30;
	int site;
	site = extToInt[clust]*MAXDECOR + decor;
	return site;
}

void CCFLIST::Add(CCFDEF &cf) 
{
	CCFDEF *pcf;
	SAFEALLOC(CCFDEF,pcf);
	*pcf = cf;
	int ncf = cf.cfIdx;
	if(cfLst[ncf] != NULL) delete cfLst[ncf];
	cfLst[ncf]=pcf;
}

void CCFLIST::RemCFLst(int ncf)
{
	if(cfLst[ncf] != NULL) delete cfLst[ncf];
	cfLst.erase(ncf);
}

int CCFLIST::Length()
{
	return (int)cfLst.size();
}

void CCFLIST::PrintDefCFs(CCTRL &tC)
{
	if(tC.prDefCFs == 0) return;
	cout << "Definition of cluster functions" << endl;
	std::map<int, CCFDEF*>::const_iterator iter;
    for (iter=cfLst.begin(); iter != cfLst.end(); ++iter) {
		if(tC.prDefCFs == 1) {
			cout << iter->second->GetCFText(tC) << endl;
		} else {
			cout << iter->second->GetCFDetails(tC) << endl;
		}
	}
}

CCFDEF & CCFLIST::operator[] (int ncf)
{
	if(cfLst[ncf]==NULL) {
		cout << "No CF defintion for CF(" << ncf << ")" << endl;
		__ERRORTHROW(ERR_NOCFDEFS); 
	}
	return *(cfLst[ncf]);
}

CCFLIST & CCFLIST::operator= (CCFLIST &cflst)
{
	Clear();
	extToInt = cflst.extToInt;
	intToExt = cflst.intToExt;
	stIdx = cflst.stIdx;

	std::map<int, CCFDEF*>::const_iterator iter;
    for (iter=cflst.cfLst.begin(); iter != cflst.cfLst.end(); ++iter) {
		Add(*(iter->second));
	}
	return *this;
}
