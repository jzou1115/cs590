#include "ccfdef.hpp"

ofstream & operator<< (ofstream& fout, CCFDEF &rhs)
{
	fout << rhs.lClust;
	fout << rhs.lDecor;
	fout << rhs.lLoCf;
	fout << rhs.lHoCf;
	fout << rhs.lHoAddSt;
	fout << rhs.lSt;
	WRITEDATA(fout,rhs.nOrd);
	WRITEDATA(fout,rhs.nSym);
	WRITEDATA(fout,rhs.cfIdx);
	return fout;
}

ifstream & operator>> (ifstream& fin, CCFDEF &rhs)
{
	fin >> rhs.lClust;
	fin >> rhs.lDecor;
	fin >> rhs.lLoCf;
	fin >> rhs.lHoCf;
	fin >> rhs.lHoAddSt;
	fin >> rhs.lSt;
	READDATA(fin,rhs.nOrd);
	READDATA(fin,rhs.nSym);
	READDATA(fin,rhs.cfIdx);
	return fin;
}

CCFDEF::CCFDEF(void) {}

CCFDEF::~CCFDEF(void) {}

CCFDEF::CCFDEF(CCFDEF &cf)
{
	lClust = cf.lClust;
	lDecor = cf.lDecor;
	lLoCf  = cf.lLoCf;
	lHoCf  = cf.lHoCf;
	lHoAddSt = cf.lHoAddSt;
	lSt  = cf.lSt;
	nOrd = cf.nOrd;
	nSym = cf.nSym;
	cfIdx = cf.cfIdx;
}
CCFDEF & CCFDEF::operator= (CCFDEF &cf)
{
	lClust = cf.lClust;
	lDecor = cf.lDecor;
	lLoCf  = cf.lLoCf;
	lHoCf  = cf.lHoCf;
	lHoAddSt = cf.lHoAddSt;
	lSt  = cf.lSt;
	nOrd = cf.nOrd;
	nSym = cf.nSym;
	cfIdx = cf.cfIdx;
	return *this;
}

void CCFDEF::Clear()
{
	lClust.SetDim(0);
	lDecor.SetDim(0);
	lLoCf.SetDim(0);
	lHoCf.SetDim(0);
	lHoAddSt.SetDim(0);
	lSt.SetDim(0);
	nOrd = 0;
	nSym = 0;
	cfIdx = 0;
}

void CCFDEF::PutClust(int clust)
{
	lClust.Add(clust);
}

void CCFDEF::SetDecor(int decor)
{
	lDecor.SetDim(0);
	lDecor.Add(decor);
}

void CCFDEF::PutDecor(int decor)
{
	lDecor.Add(decor);
}

CSTRING CCFDEF::GetCFText(CCTRL &tC)
{
	char line[MAXLINE];
	CSTRING tmpS;

	// constant cluster function
	if(nOrd == 0) {
		sprintf(line,"%6d   %s",0,"Constant CF");
		tmpS = line;
		return tmpS;
	}

	// other cluster functions
	int l;
	sprintf(line,"%6d   ",cfIdx);
	for(l=0;l<nOrd;l++){
		sprintf(line,"%s %c",line,(char)('A' + tC.dsSt(lClust[l],lDecor[l])-1));
	}
	// print sites of interaction
	for(l=0;l<nOrd;l++){
		sprintf(line,"%s %3d",line,lClust[l]+1);
	}
	tmpS = line;
	return tmpS;
}

CSTRING CCFDEF::GetCFDetails(CCTRL &tC)
{
	char line[MAXLINE];
	CSTRING tmpS;

	sprintf(line,"%6d   \n",cfIdx);
	// constant cluster function
	if(nOrd == 0) {
		sprintf(line,"%s  %-18s:",line,"Constant CF");
		tmpS = line;
		return tmpS;
	}

	// other cluster functions
	int l,k;
	sprintf(line,"%s  %-18s:",line,"Decor");
	for(l=0;l<nOrd;l++){
		sprintf(line,"%s %c",line,(char)('A' + tC.dsSt(lClust[l],lDecor[l])-1));
	}
	tmpS = line;
	sprintf(line,"\n  %-18s:","Clust");
	tmpS += line;
	for(k=0;k<nSym;k++) {
		strcpy(line,"");
		for(l=0;l<nOrd;l++){
			sprintf(line,"%s %3d",line,lClust[k*nOrd+l]+1);
		}
		if(k != nSym-1) sprintf(line,"%s;",line);
		tmpS += line;
	}
	// print sites of interaction
	for(l=0;l<nOrd;l++){
		sprintf(line,"%s %3d",line,lClust[l]+1);
	}
	sprintf(line,"\n  %-18s:", "Sites");
	tmpS += line;
	for(l=0;l<lSt.GetItsDim();l++){
		sprintf(line," %5d",lSt[l]);
		tmpS += line;
	}
	sprintf(line,"\n  %-18s:", "Lower-order CFs");
	tmpS += line;
	for(l=0;l<lLoCf.GetItsDim();l++){
		sprintf(line," %5d",lLoCf[l]);
		tmpS += line;
	}
	sprintf(line,"\n  %-18s:", "Higher-order CFs");
	tmpS += line;
	for(l=0;l<lHoCf.GetItsDim();l++){
		sprintf(line," %5d",lHoCf[l]);
		tmpS += line;
	}
	sprintf(line,"\n  %-18s:", "Sym Level for CFs");
	tmpS += line;
	for(l=0;l<lHoSmStCf.GetItsDim();l++){
		sprintf(line," %5d",lHoSmStCf[l]);
		tmpS += line;
	}
	sprintf(line,"\n  %-18s:", "Add Sites for HOCF");
	tmpS += line;
	for(l=0;l<lHoAddSt.GetItsDim();l++){
		sprintf(line," %5d",lHoAddSt[l]);
		tmpS += line;
	}
	return tmpS;
}

// mode 0: cluster size, 1: interaction amino acids, 2: interaction sites
CSTRING CCFDEF::GetCFInfo(CCTRL &tC, int mode)
{
	char line[MAXLINE];
	CSTRING tmpS;

	// constant cluster function
	if(nOrd == 0) {
		if(mode == 0) {
			sprintf(line,"%d",0);
		}else {
			sprintf(line,"%s","Constant CF");
		}
		tmpS = line;
		return tmpS;
	}

	// other cluster functions
	if(mode == 0) {
		sprintf(line,"%d",nOrd);
		tmpS = line;
		return tmpS;
	}
	// other cluster functions
	int l;
	sprintf(line,"");
	if(mode == 1) {
		for(l=0;l<nOrd;l++){
			sprintf(line,"%s %c",line,(char)('A' + tC.dsSt(lClust[l],lDecor[l])-1));
		}
	}
	// print sites of interaction
	if(mode == 2) {
		for(l=0;l<nOrd;l++){
			sprintf(line,"%s %3d",line,lClust[l]+1);
		}
	}
	tmpS = line;
	return tmpS;
}
