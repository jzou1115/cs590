#include "ccfs.hpp"

__LOGOUT
EXTERN

ofstream& operator<< (ofstream& fout, CCFS &ccfs)
{
	fout << ccfs.nAaAtEaDsStRn;   
	fout << ccfs.clIfRn;  
	fout << ccfs.nSmStOfCl;  
	fout << ccfs.clSz;   
	fout << ccfs.pmFoAlCf;    
	WRITEDATA(fout, ccfs.nTtCf);  
	fout << ccfs.clRgInCf;  
	fout << ccfs.lsSeCf;   
	fout << ccfs.sePm;   
	fout << ccfs.crVlSmr;  
	fout << ccfs.lsEci; 
	fout << ccfs.eciCo;   
	fout << ccfs.dgPj;    
	fout << ccfs.dfInEn;    
	fout << ccfs.svRs;
	fout << ccfs.elimCFsList;
	WRITEDATA(fout, ccfs.nUpCFs);
	fout << ccfs.svRsSmr;

	fout << ccfs.crVl;  
	fout << ccfs.htEciCo;
	fout << ccfs.htNrEr;
	fout << ccfs.magErr;

	fout << ccfs.ctgCFSeqs;    
	fout << ccfs.pctSeqsCFs;   
	fout << ccfs.ctgCFList;    
	fout << ccfs.ctgCFRange;   

	WRITEDATA(fout, ccfs.mnNSeq);
	fout << ccfs.lsElCf;
	fout << ccfs.lsOtCf;
	fout << ccfs.lsTtCf;

	WRITEDATA(fout, ccfs.nMScr);
	WRITEDATA(fout, ccfs.facDF1);
	WRITEDATA(fout, ccfs.facDF2);
	WRITEDATA(fout, ccfs.nSrchRemCFs);
	fout << ccfs.svRsIdx;
	fout << ccfs.cfDefs;
	fout << ccfs.useHiarc;

	return fout;
}

ifstream& operator>> (ifstream& fin, CCFS &ccfs)
{
	fin >> ccfs.nAaAtEaDsStRn;   
	fin >> ccfs.clIfRn;  
	fin >> ccfs.nSmStOfCl;  
	fin >> ccfs.clSz;   
	fin >> ccfs.pmFoAlCf;    
	READDATA(fin, ccfs.nTtCf);  
	fin >> ccfs.clRgInCf;  
	fin >> ccfs.lsSeCf;   
	fin >> ccfs.sePm;   
	fin >> ccfs.crVlSmr;  
	fin >> ccfs.lsEci; 
	fin >> ccfs.eciCo;   
	fin >> ccfs.dgPj;    
	fin >> ccfs.dfInEn;   
	fin >> ccfs.svRs;
	fin >> ccfs.elimCFsList;
	READDATA(fin, ccfs.nUpCFs);
	fin >> ccfs.svRsSmr;

	fin >> ccfs.crVl;  
	fin >> ccfs.htEciCo;
	fin >> ccfs.htNrEr;
	fin >> ccfs.magErr;

	fin >> ccfs.ctgCFSeqs;    
	fin >> ccfs.pctSeqsCFs;   
	fin >> ccfs.ctgCFList;    
	fin >> ccfs.ctgCFRange;

	READDATA(fin, ccfs.mnNSeq);
	fin >> ccfs.lsElCf;
	fin >> ccfs.lsOtCf;
	fin >> ccfs.lsTtCf;

	READDATA(fin, ccfs.nMScr);
	READDATA(fin, ccfs.facDF1);
	READDATA(fin, ccfs.facDF2);
	READDATA(fin, ccfs.nSrchRemCFs);
	fin >> ccfs.svRsIdx;
	fin >> ccfs.cfDefs;
	fin >> ccfs.useHiarc;

	return fin;
}

CCFS::CCFS() {
	nUpCFs = 0;
}

void CCFS::GenInitCFs(CCTRL &tC)
{
	// PRINTLOG(<RUN> Generate cluster function candidate);
	// PRINTMIDREMVAR();
	int nn=0;
	int i, j, k;
	int idx, tot;
	int site1, site2;
	int size1;
	CNVEC_INT cfv;

	// Evaluate initial cluster functions
	clIfRn.SetDim(tC.clIf.GetItsDim(ROW), tC.clIf.GetItsDim(COLUMN));
	nAaAtEaDsStRn.SetDim(tC.nAaAtEaDsSt.GetItsDim());
	clSz       = tC.clSz;
	nSmStOfCl  = tC.nSmStOfCl;
	mnNSeq     = tC.mnNSeq;
	mnNSeq     = tC.mnNSeq;
	nMScr      = tC.nMScr;
	facDF1     = tC.facDF1;
	facDF2     = tC.facDF2;
	nSrchRemCFs= tC.nSrchRemCFs;
	useHiarc   = tC.useHiarc;

	for(i=0;i<tC.nDS;i++) 
		nAaAtEaDsStRn[i] = tC.nAaAtEaDsSt[i]-1;          // The first residue is regarded as an amino-acid reference (renumbered)
	for(i=0;i<tC.nIS;i++) {
		nn=tC.clSz[i]*tC.nSmStOfCl[i];
		for(k=0;k<nn;k++)	clIfRn(i,k) = tC.clIf(i,k)-1;  // interaction index start from 0;
	}
	cfv.SetDim(nAaAtEaDsStRn.GetItsDim());
	tot = 0;
	for(i=0;i<clSz.GetItsDim();i++) {
		if(clSz[i] != 1) break;
		tot += nAaAtEaDsStRn[i];
		for(j=0;j<tC.nSmStOfCl[i];j++) {
			idx = clIfRn(i,j);
			cfv[idx]=i;
		}
	}

	// Find the maximum order of Clusters
	int nmax=0;
	for(i=0;i<clSz.GetItsDim();i++) {
		if(nmax < clSz[i]) nmax = clSz[i];
	}
	maxClOrd = nmax;

	// Set definitions of CFs
	cfDefs.Clear();
	// Set definition of a costant CF
	clRgInCf.SetDim(0);       // Range of CFs
	clRgInCf.Add(0);          // 0 <= constant CF < 1
	clRgInCf.Add(1);          // 1 <= point CF
	CCFDEF cf;
	cf.Clear();
	cf.nOrd=0; cf.nSym=1; cf.cfIdx=0;
	cfDefs.Add(cf);
	int decor, cfidx, clust;
	// Find maximum symmetry
	int maxsym=1;
	for(i=0;i<tC.nIS;i++) {
		if(clSz[i] != 1) continue; // only point CFs pass
		if(maxsym < tC.nSmStOfCl[i]) maxsym = tC.nSmStOfCl[i];
	}
	// Set internal cluster index 
	nn=0;
	for(k=0;k<maxsym;k++) {
		for(i=0;i<tC.nIS;i++) {
			if(clSz[i] != 1) continue; // only point CFs pass
			if(k < tC.nSmStOfCl[i]) {
				cfDefs.SetIntClIdx(clIfRn(i,k),nn);
				nn++;
			}
		}
	}
	// Set definitions for point CFs
	for(i=0;i<tC.nIS;i++) {
		if(clSz[i] != 1) continue;
		cf.Clear();
		cf.nOrd = 1;
		cf.nSym = tC.nSmStOfCl[i];
		// Set cluster information
		for(k=0;k<cf.nSym;k++) {cf.PutClust(clIfRn(i,k));}
		// Set information of decoration
		size1 = nAaAtEaDsStRn[clIfRn(i,0)];
		if(size1 == 0) continue;
		cf.lDecor.SetDim(1);
		for(j=1;j<=size1;j++) {
			cf.lSt.SetDim(0);
			decor = j;
			cf.lDecor[0] = decor;
			cfidx = cfDefs.Length();
			for(k=0;k<cf.nSym;k++) {
				clust = clIfRn(i,k);
				cfDefs.SetStIdxMap(decor,clust,cfidx);
				cf.lSt.Add(cfDefs.GetStIdx(decor,clust));
			}
			cf.cfIdx = cfidx;
			cfDefs.Add(cf);
		}
	}
	clRgInCf.Add(cfDefs.Length());       // point CF < cfDefs.Length()
	// Set definitions for remaining CFs
	CNVEC_INT cfsize;
	int ordcf, ncf, n1, n2, site, kmax, ordsm, ncf2, ordsm2;
	for(ordcf=2;ordcf <= maxClOrd; ordcf++) {
		for(i=0;i<tC.nIS;i++) {
			if(clSz[i] != ordcf) continue;
			cf.Clear();
			cf.nOrd = ordcf;
			cf.nSym = tC.nSmStOfCl[i];
			// Set cluster information
			for(k=0;k<cf.nOrd*cf.nSym;k++) {cf.PutClust(clIfRn(i,k));}
			// Set information of decoration
			cfsize.SetDim(0);
			tot=1;
			for(k=0;k<cf.nOrd;k++) {
				cfsize.Add(nAaAtEaDsStRn[clIfRn(i,k)]);
				tot *= cfsize.GetEndVal();
			}
			if(tot== 0) continue;
			cf.lDecor.SetDim(cf.nOrd);
			cf.lDecor.Fill(1);
			cf.lDecor[cf.nOrd-1]=0;
			while(1) {
				for(k=cf.nOrd;k>=1;k--) {
					if(cfsize[k-1] == cf.lDecor[k-1]) {
						cf.lDecor[k-1]=1;
					} else {
						++cf.lDecor[k-1];
						break;
					}
				}
				if(k == 0) break;
				cfidx = cfDefs.Length();
				// Set point CF list involved in this CF
				cf.lSt.SetDim(0);
				for(j=0;j<cf.nOrd;j++) {
					clust = cf.lClust[j];
					decor = cf.lDecor[j];
					ncf = cfDefs.GetStIdx(decor,clust);
					cf.lSt.Add(ncf);
				}
				cf.lSt.Sort(ASCEND);
				// Find lower-order CFs or Build hierarchy of CFs
				cf.lLoCf.SetDim(0);
				for(j=0;j<cf.nOrd;j++) {
					site = cf.lSt[j];
					for(n1=0;n1<cf.nOrd;n1++) { if(n1 != j) break;}
					site1 = cf.lSt[n1];
					ncf = cfDefs.GetCfIdx(site1);
					ordsm = cfDefs[ncf].lSt.FindElem(site1);
					if(ordsm == -1) continue;            // no such point
					for(n2=0;n2<cf.nOrd;n2++) {
						if(n2 == n1 || n2 == j) continue;
						site2 = cf.lSt[n2];
						if(ordsm !=0) {
							ncf2 = cfDefs.GetCfIdx(site2);
							ordsm2 = cfDefs[ncf2].lSt.FindElem(site2);
							site2 = cfDefs[ncf2].lSt[ordsm2-ordsm];
						}
						kmax = cfDefs[ncf].lHoAddSt.GetItsDim();
						for(k=0;k<kmax;k++) {
							if (site2 == cfDefs[ncf].lHoAddSt[k] && cfDefs[ncf].lHoSmStCf[k] == 0) break;
						}
						if(k==kmax) break;
						ncf = cfDefs[ncf].lHoCf[k];
					}
					if(cfDefs[ncf].nOrd + 1 != cf.nOrd) continue;  // No lower-order CF
					cf.lLoCf.Add(ncf);
					cfDefs[ncf].lHoCf.Add(cfidx);
					cfDefs[ncf].lHoAddSt.Add(site);
					cfDefs[ncf].lHoSmStCf.Add(ordsm);
				}
				cf.lLoCf.Sort(ASCEND);
				cf.cfIdx = cfidx;
				cfDefs.Add(cf);
			}
		}
		clRgInCf.Add(cfDefs.Length());       // Higer-order CF < cfDefs.Length()
	}
	nTtCf = clRgInCf.GetEndVal();
	// Cluster function list
	lsTtCf.SetDim(nTtCf);
	lsTtCf.FillInc();
}
double CCFS::CalcChangeInCVRMS()
{
	int n;
	n = crVlSmr.GetItsDim();
	double value = fabs(crVlSmr[n-1]-crVlSmr[n-2])/crVlSmr[n-1];
	char tmpstr[MAXLINE];
	sprintf(tmpstr,"CVRMS variation by changing CFs %9.4f %9.4f %9.4f",crVlSmr[n-2],crVlSmr[n-1],value);
	cout << tmpstr << endl;
	return value;
}

//determine order of cluster functions for screening
//based on decrease in CVRMS value by applying linear regression with each CF
void CCFS::DtOrdOfCFs1st(CTSET &tSet, int dimcf, CNVEC_INT& ordv)
{
	// Calculate screening order of cluster functions
	CNVEC_DBL tmDlist;
	CNVEC_INT tmNu, tmDCFs, tmElList;
	CNVEC_DBL tmdblQO, tmQN, tmEdN, tmVec1, tmVec2, tmEDO, tmQ;
	double tmMu, tmT, tmD, tmCvrms, tm_dblPO, dn;
	int ncf, n, ntotcf;

	int nset = tSet.GetNSeq(D_ITS);
	ntotcf = lsTtCf.GetItsDim();
	if (nset == 0 || ntotcf==0) { ordv.SetDim(0); return;}

	dn = (double) nset;
	tm_dblPO = 1.0/dn;
	tmdblQO.SetDim(nset);
	tmdblQO = 1.0/dn;
	tmEDO = tSet.itsEn - tSet.itsEn.Average();
	//tmVec1 = tmEDO.DirectDot(tmEDO);
	//cvrms = dn*sqrt(tmVec1.Average())/(dn-1.0);
	tmQ.SetDim(nset);
	tmDlist.SetDim(0);
	tmDCFs.SetDim(0);
	for(n=0; n<ntotcf; n++) {
		ncf = lsTtCf[n];
		CCFDEF &cf = cfDefs[ncf];
		if ( cf.nOrd != dimcf) continue;
		CalcSCD(tSet,ncf,tmNu);
		if (tmNu.CountElemExcept(0) <= mnNSeq) { 
			tmElList.Add(ncf); continue; 
		} // screening by using sequences
		tmMu = tmNu.Sum();
		tmT = tm_dblPO*tmMu;
		tmD  = tmNu*tmNu - tmMu*tmT;
		if(tmD < INVCUTOFF) {
			tmElList.Add(ncf);
		} else {
			tmQ = tmT;
			tmVec1 = tmQ-tmNu;
			tmQN   = tmdblQO + tmVec1.DirectDot(tmVec1)/tmD;
			tmEdN  = tmEDO + tmVec1*(tmNu*tmEDO)/tmD;
			tmVec2 = tmEdN/(-tmQN+1.0);       // direct divide
			tmVec2 = tmVec2.DirectDot(tmVec2);
			tmCvrms = sqrt(tmVec2.Average());
			tmDlist.Add(tmCvrms);
			tmDCFs.Add(ncf);
		}
	}
	int uhier_old = useHiarc;
	useHiarc = 1;
	for(n=0; n<tmElList.GetItsDim(); n++) {
		ncf = tmElList[n];
		MoveCFs(ncf, 2);
	}
	useHiarc = uhier_old;
	// Calculate screening order of cluster functions
	if(tmDlist.GetItsDim() == 0) ordv.SetDim(0);
	else ordv = tmDCFs.GetSelect(tmDlist.GetSortIdx(ASCEND));
}

//determine order of cluster functions for screening
//based on decrease in CVRMS value by applying multiple linear regression with an additional CF
void CCFS::DtOrdOfCFsLast(CTSET &tSet, int dimcf, CNVEC_INT &ordv)
{
	// Calculate screening order of CFs
	CNVEC_INT tmDCFs, tmNu, tmMu, tmElList;
	CNVEC_DBL tmDlist, tmQ, tmQN, tmEdN, tmVec1, tmVec2, tmT;
	int ncf, n;
	double tmD, tmCvrms, cvrms;

	int nset = tSet.GetNSeq(D_ITS);
	if (nset == 0 ) {ordv.SetDim(0); return ;}

	tmQ.SetDim(nset);
	tmDlist.SetDim(0);
	tmDCFs.SetDim(0);
	cvrms = crVl.GetEndVal();
	for(n=0; n<lsTtCf.GetItsDim(); n++) {
		ncf = lsTtCf[n];
		CCFDEF &cf = cfDefs[ncf];
		if ( cf.nOrd != dimcf) continue; 
		CalcSCD(tSet,ncf,tmNu);
		if (tmNu.CountElemExcept(0) <= mnNSeq) {
			tmElList.Add(ncf); continue; 
		}
		tmMu = tmNu*sePm;
		tmT  = eciCo*tmMu;
		tmD  = tmNu*tmNu - tmMu*tmT;
		if(tmD < INVCUTOFF) {
			tmElList.Add(ncf);
		} else {
			tmQ = sePm*tmT;
			tmVec1 = tmQ-tmNu;
			tmQN = dgPj + tmVec1.DirectDot(tmVec1)/tmD;
			tmEdN = dfInEn + tmVec1*(tmNu*dfInEn)/tmD;
			tmVec2 = tmEdN/(-tmQN+1.0);       // direct divide
			tmVec2 = tmVec2.DirectDot(tmVec2);
			tmCvrms = sqrt(tmVec2.Average());
			tmDlist.Add(tmCvrms);
			tmDCFs.Add(ncf);
		}
	}
	for(n=0; n<tmElList.GetItsDim(); n++) {
		ncf = tmElList[n];
		MoveCFs(ncf, 2);
	}
	// Calculate screening order of cluster functions
	if(tmDlist.GetItsDim() == 0) ordv.SetDim(0);
	else ordv = tmDCFs.GetSelect(tmDlist.GetSortIdx(ASCEND));
}

void CCFS::MoveCFs(int ncf, int mode)
{
	// mode 0: move relevant CFs to lsSeCf
	// mode 1: move relevant CFs to lsOtCf
	// mode 2: move relevant CFs to lsElCf
	// mode 3: permanently remove CFs
	CNVEC_INT tmV;
	int n;
	// removed already
	if(lsSeCf.FindElem(ncf) != -1 || lsElCf.FindElem(ncf) != -1 || lsOtCf.FindElem(ncf) != -1) return;
	CCFDEF &cf = cfDefs[ncf];
	switch (mode) {
		case 0:
			lsSeCf.Add(ncf);
			break;
		case 1:
			tmV = lsTtCf.RemElem(ncf);
			lsOtCf.Add(tmV);
			// if point CF is eliminated then corresponding pair CFs are also eliminated
			if (useHiarc == 1) { 
				for(n=0;n<cf.lHoCf.GetItsDim();n++) MoveCFs(cf.lHoCf[n],1);
			}
			break;
		case 2:
			tmV = lsTtCf.RemElem(ncf);
			lsElCf.Add(tmV);
			// if point CF is eliminated then corresponding pair CFs are also eliminated
			if (useHiarc == 1) { 
				for(n=0;n<cf.lHoCf.GetItsDim();n++) MoveCFs(cf.lHoCf[n],2);
			}
			break;
		case 3:
			// if point CF is eliminated then corresponding pair CFs are also eliminated
			if (cf.nOrd == 1) { 
				tmV = lsTtCf.RemElem(cf.lHoCf);
			}
			break;
	}
}

// Calculate probe matrix from a test set.
double CCFS::CalcCEEnergy(CNVEC_INT &vec)
{
	double energy;
	int n, ncf;
	int add,j,k,decor,clust;
	energy = lsEci[0];
	for(n=1;n<lsSeCf.GetItsDim();n++) {
		ncf = lsSeCf[n];
		CCFDEF &cf = cfDefs[ncf];
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				decor = cf.lDecor[k];
				clust = cf.lClust[j*cf.nOrd+k];
				if(vec[clust] != decor) {add = 0; break;}
			}
			if(add == 1) energy += lsEci[n];
		}	
	}
	return energy;
}

void CCFS::CalcEnDiff(CNVEC_DBL &env, CNMAT_INT &seqm)
{
	CNVEC_INT tmpV;
	CNVEC_DBL dblV;
	int nn = seqm.GetItsDim(ROW);
	dblV.SetDim(nn);
	char str[MAXLINE];
	sprintf(str,"\n %12s %12s %12s","Calc","CEExpect","Diff");
	cout << str << endl;  
	for(int i=0;i<nn;i++) {
		tmpV = seqm.GetVec(COLUMN,i);
		double en = CalcCEEnergy(tmpV);
		double diff = env[i] - en;
		sprintf(str," %12.3f %12.3f %12.5f",env[i],en,diff);
		cout << str << endl;  
	}
}

 
/* 
   sePm contains selected characteristic data of system(CDS) (constant, points, pairs )
   lsSeCf contains selected list of CDS (points, pairs)
*/

int CCFS::UpdateCFs(CCTRL &tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Update cluster functions);
	PRINTMIDREMVAR();

	int nset = tSet.GetNSeq(D_ITS);
	if (nset ==0) { 
		cout << "Error: No remaining sequences in the function UpdateCFs" << endl;
		exit(0);
	}
	CNVEC_INT tmNu, tmMu, tmTmp;
	CNVEC_INT ordV;
	CNVEC_DBL  tmT, tmQ, tmQN, tmEdN, tmVec1, tmVec2, tmVec3;
	double tmD, tmcvrms, dn, cvrms;
	int ncf, n, ntot, nord;
	bool bsel, ok;

	
	// For updating CF list, add old CF list to CF pool
	if (lsSeCf.GetItsDim() != 0) {
		lsTtCf.Add(lsSeCf);
		lsTtCf.Add(lsOtCf);
		lsTtCf.Sort(ASCEND);
	}
	lsSeCf.SetDim(0);
	lsOtCf.SetDim(0);

	// Build a constant cluster functions
	dn = (double) nset;
	eciCo.SetDim(1,1);
	eciCo(0,0) = 1.0/dn;
	dgPj.SetDim(nset);
	dgPj = 1.0/dn;
	dfInEn = tSet.itsEn - tSet.itsEn.Average();
	CNVEC_DBL tmVec = dfInEn.DirectDot(dfInEn);
	crVl.Add(dn*sqrt(tmVec.Average())/(dn-1.0));
	MoveCFs(0,0); // move CF "0" to lsSeCf
	tmNu.SetDim(nset);
	tmNu.Fill(1);
	sePm.SetDim(0,0);
	sePm.Add(COLUMN,tmNu);

	// Build point CFs
	DtOrdOfCFs1st(tSet,1, ordV);  // Get a list of ordered point CFs
	ntot = ordV.GetItsDim();
	for(n=0; n<ntot;n++) {
		ncf = ordV[n];
		CalcSCD(tSet, ncf, tmNu);
		tmMu = tmNu*sePm;
		tmT = eciCo*tmMu;
		tmD = tmNu*tmNu - tmMu*tmT;
		if(tmD < INVCUTOFF) {
			MoveCFs(ncf, 1);
		} else {
			// q = M(:, 1:end-1)*t
			tmQ = sePm*tmT;
			// QO = QO + ((q-nu).^2)/d;
			tmVec1 = tmQ-tmNu;
			tmQN = dgPj + tmVec1.DirectDot(tmVec1)/tmD;
			// edO = edO + (q-nu)*(nu'*edO)/d;
			tmEdN = dfInEn + tmVec1*(tmNu*dfInEn)/tmD;  
			// cv_rms_p = sqrt(mean((edN./(1-QN)).^2));
			// cv_rms = [cv_rms cv_rms_p]
			tmVec2 = tmEdN/(-tmQN+1.0);       // direct divide
			tmVec2 = tmVec2.DirectDot(tmVec2);
			tmcvrms = sqrt(tmVec2.Average());
			if (nMScr == 0) {  // determine whether applying screening procedure with CVRMS
				double nufree = (double)(sePm.GetItsDim(ROW)-sePm.GetItsDim(COLUMN))*2.0;
				cvrms = crVl.GetEndVal();
				if ((cvrms - tmcvrms)*nufree > cvrms*facDF1) bsel = true;
				else bsel =false;
			} else bsel = true;
			if (bsel) {
				crVl.Add(tmcvrms);
				dgPj = tmQN;
				dfInEn = tmEdN;
				// PO = [ PO + t*t'/d   - t/d; -t'/d  1/d]
				eciCo += tmT.Dyadic(tmT)/tmD;
				tmVec3 = tmT/(-tmD);
				eciCo.Add(COLUMN, tmVec3);
				tmVec3.Add(1.0/tmD);
				eciCo.Add(ROW, tmVec3);
				// ac = [ac i]
				MoveCFs(ncf, 0);
				sePm.Add(COLUMN,tmNu);
			} else {
				MoveCFs(ncf, 1);
			}
		}
	}
	//	Build higher order CFs
	int i;
	for(nord=2;nord<=maxClOrd;nord++) {
		DtOrdOfCFs1st(tSet,nord, ordV);
		ntot = ordV.GetItsDim();
		for(n=0; n<ntot;n++) {
			int ncf = ordV[n];
			if(useHiarc == 1) {        // Check all lower order CFs were selected
				CCFDEF & cf = cfDefs[ncf];
				ok = true;
				if(cf.lLoCf.GetItsDim() != nord) continue;
				for(i=0;i<cf.lLoCf.GetItsDim();i++) {
					if(lsSeCf.FindElem(cf.lLoCf[i]) == -1){ 
						ok = false; break;
					}
				}
				if(!ok) continue;
			}
			CalcSCD(tSet, ncf, tmNu);
			// M(:,1:end-1)'*nu;
			tmMu = tmNu*sePm;
			tmT = eciCo*tmMu;
			tmD = tmNu*tmNu - tmMu*tmT;
			if(tmD < INVCUTOFF) {
				MoveCFs(ncf, 1);
			} else {
				// q = M(:, 1:end-1)*t
				tmQ = sePm*tmT;   
				// QN = QO + ((q-nu).^2)/d;
				tmVec1 = tmQ - tmNu;
				tmQN = dgPj + tmVec1.DirectDot(tmVec1)/tmD;
				// edN = edO + (q-nu)*(nu'*edO)/d;
				tmEdN = dfInEn + tmVec1*(tmNu*dfInEn)/tmD;  
				// cv_rms_p = sqrt(mean((edN./(1-QN)).^2));
				// cv_rms = [cv_rms cv_rms_p]
				tmVec2 = tmEdN/(-tmQN+1.0);       // direct divide
				tmVec2 = tmVec2.DirectDot(tmVec2);
				tmcvrms = sqrt(tmVec2.Average());
				cvrms = crVl.GetEndVal();
				double nufree = (double)(sePm.GetItsDim(ROW)-sePm.GetItsDim(COLUMN))*2.0;
				if ((cvrms - tmcvrms)*nufree > cvrms*facDF2) {
					crVl.Add(tmcvrms);
					dgPj = tmQN;
					dfInEn = tmEdN;
					// PO = [ PO + t*t'/d   - t/d; -t'/d  1/d]
					eciCo += tmT.Dyadic(tmT)/tmD;
					tmVec3 = tmT/(-tmD);
					eciCo.Add(COLUMN, tmVec3);
					tmVec3.Add(1.0/tmD);
					eciCo.Add(ROW, tmVec3);
					// ac = [ac i]
					MoveCFs(ncf, 0);
					sePm.Add(COLUMN,tmNu);
				} else {
					MoveCFs(ncf, 1);
				}
			}
		}
	}

	lsEci = eciCo*(tSet.itsEn*sePm);
	// calculate fluctuation of eci values
	double ratio = (double)(tSet.GetNSeq(D_ITS) - sePm.GetItsDim(COLUMN))/(double)sePm.GetItsDim(COLUMN);
	++nUpCFs;
	// print CVRMS value
	cout << "\t"<< "Number of Step               : " << nUpCFs << endl;
	cout << "\t"<< "CVRMS value                  : " << crVl.GetEndVal() << endl;
	cout << "\t"<< "Number of sequences          : " << tSet.itsEn.GetItsDim() << endl; 
	cout << "\t"<< "Number of selected CFs       : " << lsSeCf.GetItsDim() << endl; 
	cout << "\t"<< "Remaining degrees of freedom : " << ratio << endl;
	// save data
	crVlSmr.Add(crVl.GetEndVal());
	int nsave = svRsSmr.GetItsDim(ROW);
	svRsSmr.SetDim(nsave+1,14);
	svRsSmr(nsave,0) = (double) nUpCFs;
	svRsSmr(nsave,1) = (double) tSet.itsEn.GetItsDim();
	svRsSmr(nsave,2) = (double) lsSeCf.GetItsDim();
	svRsSmr(nsave,9) = ratio;
	svRsSmr(nsave,13) = sqrt(dfInEn*dfInEn/(double)dfInEn.GetItsDim());
	svRsSmr(nsave,3) = crVl.GetEndVal();
	// Calc magnitude of error distribution vector
	if (dfInEn.Magnitude() == 0.0) {
		svRsSmr(nsave,4)=0.0;
	}else {
		tmVec1 = dfInEn/(1.0-dgPj);
		tmVec1 = tmVec1.UnitVec();
		tmVec1 = tmVec1.DirectDot(tmVec1);
		tmVec1 = tmVec1*dfInEn.GetItsDim()-1.0;
		svRsSmr(nsave,4) = tmVec1.Magnitude();
	}

	return 0;                                           // adding more sequnces
}

void CCFS::UpdateCFs(CCTRL &tC, CTSET &tSet, CNVEC_INT listCFs)
{
	PRINTLOG(<RUN> Update cluster functions);
	PRINTMIDREMVAR();

	// build single cluster functions
	eciCo.SetDim(1,1);
	int nset = tSet.GetNSeq(D_ITS);
	double dn = (double) nset;
	eciCo(0,0) = 1.0/dn;
	dgPj.SetDim(nset);
	dgPj = 1.0/dn;
	dfInEn = tSet.itsEn - tSet.itsEn.Average();
	CNVEC_DBL tmVec = dfInEn.DirectDot(dfInEn);
	crVl.Add(dn*sqrt(tmVec.Average())/(dn-1.0));
	lsSeCf.SetDim(1);
	lsSeCf[0] = 0;

	CNVEC_INT tmNu, tmMu;
	tmNu.SetDim(nset);
	tmNu.Fill(1);
	sePm.SetDim(0,0);
	sePm.Add(COLUMN,tmNu);

	CNVEC_DBL  tmT, tmQ, tmQN, tmEdN, tmVec1, tmVec2, tmVec3;
	double tmD, tmcvrms;

	int ncf, n;
	for(n=1; n<listCFs.GetItsDim();n++) {
		ncf = listCFs[n];
		CalcSCD(tSet, ncf, tmNu);
		tmMu = tmNu*sePm;
		tmT = eciCo*tmMu;
		tmD = tmNu*tmNu - tmMu*tmT;
		if(tmD > INVCUTOFF) {
			// q = M(:, 1:end-1)*t
			tmQ = sePm*tmT;
			// QO = QO + ((q-nu).^2)/d;
			tmVec1 = tmQ-tmNu;
			tmQN = dgPj + tmVec1.DirectDot(tmVec1)/tmD;
			dgPj = tmQN;
			// edO = edO + (q-nu)*(nu'*edO)/d;
			tmEdN = dfInEn + tmVec1*(tmNu*dfInEn)/tmD;  
			dfInEn = tmEdN;
			// cv_rms_p = sqrt(mean((edN./(1-QN)).^2));
			// cv_rms = [cv_rms cv_rms_p]
			tmVec2 = tmEdN/(-tmQN+1.0);       // direct divide
			tmVec2 = tmVec2.DirectDot(tmVec2);
			tmcvrms = sqrt(tmVec2.Average());
			crVl.Add(tmcvrms);
			// PO = [ PO + t*t'/d   - t/d; -t'/d  1/d]
			eciCo += tmT.Dyadic(tmT)/tmD;
			tmVec3 = tmT/(-tmD);
			eciCo.Add(COLUMN, tmVec3);
			tmVec3.Add(1.0/tmD);
			eciCo.Add(ROW, tmVec3);
			// ac = [ac i]
			lsSeCf.Add(ncf);
			sePm.Add(COLUMN,tmNu);
		}
	}
	lsEci = eciCo*(tSet.itsEn*sePm);
}

// If choose a bad interaction then return true, else return false
bool CCFS::ChooseErrInteract(CCTRL & tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Choose a cluster function with highest NH-ratio);
	PRINTMIDREMVAR();
	// if there are no sequences then return
	if(tSet.itsEn.GetItsDim() == 0 ) return false;

	if(maxClOrd < 2) return false;

	// Calculate error of interaction for selected cluster functions
	CNVEC_DBL  tmvCos;
	CNVEC_INT tmV2, tmlCFs;
	int ncf, ni, nc;
	double cosval;

	// Choose interaction for elimination in pair cluster functions
	nc = svRs.GetItsDim(COLUMN)-1;
	tmvCos = svRs.GetVec(COLUMN,nc);
	nc = svRsIdx.GetItsDim(COLUMN)-1;
	tmlCFs = svRsIdx.GetVec(COLUMN,nc);
	tmV2 = tmvCos.GetSortIdx(DESCEND);
	for(ni=0;ni<tmV2.GetItsDim();ni++) {
		ncf = tmlCFs[tmV2[ni]];
		CCFDEF &cf = cfDefs[ncf];
		if(cf.nOrd > 1) break;
	}
	if (ni >= tmV2.GetItsDim()) {
		cout << "\t" << "No additional bad interaction is selected" << endl;
		return false;
	}
	cosval = tmvCos[tmV2[ni]];

	// eliminating sequence cutoff value
	switch(tC.elimMode) {
		case 0:
			if(cosval < tC.cosChiCut) {
				cout << "\t" << "No additional bad interaction is selected" << endl;
				cout << "  " << cosval << "/" << tC.cosChiCut << "  " << cfDefs[ncf].GetCFText(tC) << endl;
				return false;
			} else {
				cout << "\t"<< "Selected bad interaction is" << endl;
			}
			break;
		case 1:
			if(elimCFsList.GetItsDim() >= tC.nElimCFs) {
				cout << "\t" << "Largest NH-ratio in pair or higher-order CFs" << endl;
				CSTRING tmpS;
				tmpS = cfDefs[ncf].GetCFText(tC);
				cout << "  " << cosval << "  " << tmpS  << endl;
				return false;
			} else {
				cout << "\t" << "selected CF for elimination" << endl;
			}
			break;
	}
	int nsave = svRsSmr.GetItsDim(ROW)-1;
	elimCFsList.Add(ncf);
	MoveCFs(ncf, 3);
	int nni = lsSeCf.FindElem(ncf);
	if(nni != -1) nni = 1;
	svRsSmr(nsave,11) = (double) nni;
	cout << "  " << cosval << "  " << cfDefs[ncf].GetCFText(tC) << endl;
	cout << "\t"<< "Number of eliminated interactions : " << elimCFsList.GetItsDim() << endl;

	return true;
}

// calculate error angles of CFs 
void CCFS::CalcRemCFsECIs(CCTRL &tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Calculation of NH-ratios);
	PRINTMIDREMVAR();
	// if there are no sequences then return
	if(tSet.itsEn.GetItsDim() == 0 ) return;
	// all ECI values of other CFs are zero then return;
	if(dfInEn.Magnitude() == 0.0) return;

	CalcVarECIs();

	// Calculate error of interaction for selected cluster functions
	int n, ncf, ntot, row, col, col1;
	CNVEC_DBL  tmvT, tmvG, tmvECI, tmvDGG, tmvED, tmvEE, tmvDA, tmvQ, tmvPR, tmvCos;
	CNVEC_DBL  tmList, tmPR1, tmDA1, tmCos1;
	CNVEC_INT tmvMu, tmvNu, tmvOrd, tmvIdx;
	double hratio, errcos, cvrms, dblnseq, pr, tmD;


	// Calculate ECI fluctuation and error cosine in lsSeCf
	cvrms = crVlSmr.GetEndVal();
	tmvIdx = lsSeCf;
	ntot = lsSeCf.GetItsDim();
	for(n=0; n<ntot;n++) {
		ncf = lsSeCf[n];
		tmvG = sePm*eciCo.GetVec(COLUMN,n);
		tmvG = tmvG/sqrt(eciCo(n,n));
		tmvDGG = tmvG.DirectDot(tmvG);
		pr = 1.0/(tmvDGG*tmvDGG);
	    hratio = htEciCo(n,n)/(cvrms*cvrms*eciCo(n,n));
	    errcos = (hratio-1.0)*sqrt(pr)/magErr;
		tmvECI.Add(lsEci[n]);
		tmvDA.Add(hratio);
		tmvPR.Add(pr);
		tmvCos.Add(errcos);
	}
	tmPR1 = tmvPR;
	tmDA1 = tmvDA;
	tmCos1 = tmvCos;
	// Making list of CFs in lsOtCf
	if(nSrchRemCFs >0) {
		ntot = lsOtCf.GetItsDim();
		tmList.SetDim(ntot);
		tmvDGG = dfInEn.DirectDot(dfInEn);
		for(n=0; n<ntot;n++) {
			ncf = lsOtCf[n];
			CalcSCD(tSet, ncf, tmvNu);
			if (tmvNu.CountElemExcept(0) < mnNSeq ) {
				tmList[n] = -1.0;
				continue;
			}
			tmvT = tmvNu.CvtDbl();
			tmvT = tmvT.UnitVec();
			tmvT = tmvT.DirectDot(tmvT);
			tmList[n] = tmvDGG*tmvT;
		}
		tmvOrd = tmList.GetSortIdx(DESCEND);

		// Calculate error of interaction for not selected cluster functions
		dblnseq = (double)sePm.GetItsDim(ROW);

		ntot = tmvOrd.GetItsDim();
		if (nSrchRemCFs > ntot) nSrchRemCFs = ntot;
		for(n=0; n<nSrchRemCFs;n++) {
			ncf = lsOtCf[tmvOrd[n]];
			tmvIdx.Add(ncf);
			CalcSCD(tSet, ncf, tmvNu);
			tmvMu = tmvNu*sePm;
			tmvT = eciCo*tmvMu;
			tmvG = tmvNu.CvtDbl() - sePm*tmvT;
			tmD = tmvG*tmvG;
			if(tmD > INVCUTOFF) {
				// Calculate ECI value
				tmvECI.Add(tmvG*tSet.itsEn/tmD);
				// Calculate error in energy
				tmvG = tmvG/sqrt(tmD);
				tmvED = dfInEn - (tmvG*dfInEn)*tmvG;
				tmvDGG = tmvG.DirectDot(tmvG);
				pr = tmvDGG.Magnitude();
				pr = 1.0/(pr*pr);
				tmvQ = 1.0 - dgPj - tmvDGG;
				tmvEE = tmvED/tmvQ;
				tmvEE = tmvEE.UnitVec();
				tmvEE = tmvEE.DirectDot(tmvEE)*dblnseq - 1.0;
				hratio = 1.0 + tmvDGG*tmvEE;
				errcos = (hratio-1.0)*sqrt(pr)/tmvEE.Magnitude();
				tmvDA.Add(hratio);
				tmvPR.Add(pr);
				tmvCos.Add(errcos);
			} else { 
				tmvECI.Add( -100.0);
				tmvDA.Add(-1.0);
				tmvPR.Add(0.0);
				tmvCos.Add(-2.0);
			}
		}
	}
	// Save results
	tmvOrd = tmvIdx.GetSortIdx(ASCEND);
	tmvIdx = tmvIdx.GetSelect(tmvOrd);
	tmvECI = tmvECI.GetSelect(tmvOrd);
	tmvDA  = tmvDA.GetSelect(tmvOrd);
	tmvPR  = tmvPR.GetSelect(tmvOrd);
	tmvCos = tmvCos.GetSelect(tmvOrd);
	row = svRs.GetItsDim(ROW);        // save results on the variable "svRs"
	col = svRs.GetItsDim(COLUMN);
	col1 = svRsIdx.GetItsDim(COLUMN);
	n = tmvOrd.GetItsDim();
	if( row > n) {
		tmvIdx.SetDim(row);
		tmvECI.SetDim(row);
		tmvDA.SetDim(row);
		tmvPR.SetDim(row);
		tmvCos.SetDim(row);
	} else if (row < n) {
		svRs.SetDim(n, col);
		svRsIdx.SetDim(n, col1);
	}
	svRsIdx.Add(COLUMN,tmvIdx);
	svRs.Add(COLUMN,tmvECI);
	svRs.Add(COLUMN,tmvDA);
	svRs.Add(COLUMN,tmvPR);
	svRs.Add(COLUMN,tmvCos);

	int nsave = svRsSmr.GetItsDim(ROW)-1;
	svRsSmr(nsave,5) = tmPR1.Average();
	svRsSmr(nsave,6) = tmCos1.Max();
	svRsSmr(nsave,7) = tmDA1.Max();
	svRsSmr(nsave,8) = tmDA1.Skewness();

}

void CCFS::CalcStdErrData()
{
	double cvrms = crVlSmr.GetEndVal();
	double gammahetero, heterosig;

	htNrEr = dfInEn/(1. - dgPj);
	CNVEC_INT sortV, tmNu;
	for(int i=0;i<dfInEn.GetItsDim();i++) {
		tmNu = sePm.GetVec(ROW,i);
		gammahetero = tmNu*htEciCo*tmNu;
		heterosig = sqrt(cvrms*cvrms + gammahetero);
		htNrEr[i] = fabs(htNrEr[i])/heterosig;
	}
	sortV = htNrEr.GetSortIdx(ASCEND);
	cout << htNrEr << endl;
	cout << sortV << endl;
	htNrEr = htNrEr.GetSelect(sortV);
	cout << htNrEr << endl;
}

// Calculate covariance matrix of ECI values
int CCFS::CalcVarECIs()
{
	// if there are no sequences then return
	if(sePm.GetItsDim(ROW) == 0 ) return 0;
	double cvrms = crVlSmr.GetEndVal();

	CNVEC_DBL tmpV, tmpVUnit;
	if(dfInEn.Magnitude() == 0.0) {
		htEciCo.SetDim(eciCo.GetItsDim(ROW),eciCo.GetItsDim(COLUMN));
		htEciCo.SetZero();
	}else {
		tmpV = dfInEn/(1.0 - dgPj);
		tmpV = tmpV.DirectDot(tmpV);
		tmpVUnit = tmpV/(cvrms*cvrms) - 1.0;
		magErr = tmpVUnit.Magnitude();
		htEciCo = eciCo*sePm.SMTransMat(tmpV)*eciCo;
	}
	return 1;
}

// Check each sequence whether the sequnce contains eliminated cluster function (CF)
// Make a list of sequences containing last bad CF
// Move such sequences from itsSet to elimSet
void CCFS::ElimSeqs(CCTRL &tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Eliminate sequences containing erroneous interactions);
	PRINTMIDREMVAR();
	// if there are no sequences then return
	if(tSet.itsEn.GetItsDim() == 0 ) return;

	// Make a list for eliminated sequences
	int nseq;
	CNVEC_INT tmvSeq, tmvList, tmvBad;
    // Check only last selected CFs
	CalcSCD(tSet,elimCFsList.GetEndVal(), tmvBad);
	for(nseq=0;nseq<tmvBad.GetItsDim();nseq++) {
		if(tmvBad[nseq]!=0) tmvList.Add(nseq);
	}
	if(tmvList.GetItsDim() == 0) {__ERRORTHROW(ERR_ELIMSEQS); }
	// Choose most erroneous sequence
	CNVEC_INT tmvErrList;
	CNVEC_DBL tmvErr;
	tmvErr = dfInEn.GetSelect(tmvList);
	tmvErr = tmvErr.DirectDot(tmvErr);
	tmvErrList = tmvErr.GetSortIdx(DESCEND);
	double en;
	tSet.GetSeq(D_ITS,tmvList[tmvErrList[0]],tmvSeq, en);   // a sequence with worst error in energy
	tSet.sel_BadSeqs.Add(COLUMN,tmvSeq);
	cout << "\t"<< "Bad sequence\n";
	cout << "\t"<< tSet.GetSeqText(tC,tmvSeq) << endl;


	// copy selected sequences to elimSet
	tSet.MoveSeq(D_ITS,tmvList,D_ELIM);

	// eliminate sequences containing the bad interaction
	cout << "\t"<< "Total number of eliminated sequences containing the removed interaction : ";
	cout << tmvList.GetItsDim() << endl;
	int nsave = svRsSmr.GetItsDim(ROW)-1;
	svRsSmr(nsave,10) = (double) elimCFsList.GetEndVal();
	svRsSmr(nsave,12) = (double) tmvList.GetItsDim();

}

void CCFS::SequenceFilter(CTSET &tSet)
{
	PRINTLOG(<RUN> Exclude CFs which are appeared in a traing set less than criterion);
	PRINTMIDREMVAR();

	if(tSet.itsEn.GetItsDim() == 0 ) return;

	int n, nord, ncf;
	CNVEC_INT tmNu, tmvList;
	int uhier_old;

	uhier_old = useHiarc;
	useHiarc = 1;

	// Check all point CFs are involved in training sequences more than mnNSeq 
	for(nord=1;nord<=maxClOrd;nord++) {
		tmvList.SetDim(0);
		for(n=1; n<lsTtCf.GetItsDim(); n++) {
			ncf = lsTtCf[n];
			CCFDEF &cf = cfDefs[ncf];
			if(cf.nOrd != nord) continue;
			CalcSCD(tSet,ncf, tmNu);
			if ( mnNSeq >= tmNu.CountElemExcept(0)) {
				tmvList.Add(ncf);
			}
		}
		for(n=0; n < tmvList.GetItsDim();n++) {
			ncf = tmvList[n];
			MoveCFs(ncf, 1);
		}
	}
	useHiarc = uhier_old;
}

void CCFS::AnalSeqErr(CCTRL &tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Error in energy prediction of sequences);
	PRINTMIDREMVAR();

	int  ncf, ni, nc, nj, nl, ns, mi;

	// if there are no sequences then return
	ns = tSet.itsEn.GetItsDim();
	if( ns == 0 ) return;

	// Calculate error of interaction for selected cluster functions
	CNVEC_DBL  tmvCos;
	CNVEC_INT   tmV1;
	// Make a list for eliminated sequences
	double di, en;
	double cosval;
	int nni, nseq;
	CNVEC_INT tmvSeq, tmvList, tmvOrd;
	int nsave;

	// Choose a CF for elimination in pair cluster functions
	nc = svRs.GetItsDim(COLUMN)-1;
	tmvCos = svRs.GetVec(COLUMN,nc);
	tmvOrd = tmvCos.GetSortIdx(DESCEND);
	// Eliminate constant and point CFs in the list
	nl = tmvOrd.GetItsDim();
	mi = 0;
	for(ni=0;ni<nl;ni++) {
		ncf = tmvOrd[mi];
		CCFDEF &cf = cfDefs[ncf];
		if(cf.nOrd <= 1) {
			tmvOrd.Remove(mi);
		} else mi++;
	}
	// no higher order CFs than point CFs
	nl = tmvOrd.GetItsDim();
	if(nl == 0 ) return;
	ncf =tmvOrd[0];
	nsave = svRsSmr.GetItsDim(ROW)-1;
//	elimCFsList.Add(ncf);
	nni = lsSeCf.FindElem(ncf);
	if(nni != -1) nni = 1;
	cosval = tmvCos[ncf];
	svRsSmr(nsave,11) = (double) nni;
	svRsSmr(nsave,13) = sqrt(dfInEn*dfInEn/(double)ns);
	cout << "\t"<< "Selected CF not fitting to the structure : " << endl;
	cout << "ERRCOS: " << cosval << ", CF : " << cfDefs[ncf].GetCFText(tC) << endl;

    // Categorize each sequence as the order of Error cosine.
	ctgCFSeqs.SetDim(ns);
	pctSeqsCFs.SetDim(nl,3);
	pctSeqsCFs.SetZero();
	int nscf = clIfRn.GetItsDim(ROW);
	int add, decor, clust,k,j;
	int *rv;
	for(nseq = 0; nseq< ns;nseq++) {
		tSet.GetSeq(D_ITS,nseq,tmvSeq, en);
		rv = tmvSeq.GetItsVec();
		for(ni=0;ni<nl;ni++){
			nni = tmvOrd[ni];
			CCFDEF &cf = cfDefs[nni];
			for(j=0;j<cf.nSym;j++) {
				add=1;
				for(k=0;k<cf.nOrd;k++) {
					decor = cf.lDecor[k];
					clust = cf.lClust[j*cf.nOrd+k];
					if(rv[clust] != decor) {add = 0; break;}
					//if(tSet.itsSet(nseq,clust) != decor) {add = 0; break;}
				}
				if(add==1) {
					break;
				}
			}
			if(add==1) {
				ctgCFSeqs[nseq] = ni;
				if(ni==0) {tmvList.Add(nseq);};
				break;
			}
		}
		ni = ctgCFSeqs[nseq];
		pctSeqsCFs(ni,0) += 1.0;
		pctSeqsCFs(ni,1) += dfInEn[nseq]*dfInEn[nseq];
	}
	double ntset = (double)ns/100.0;
	pctSeqsCFs(0,2)=pctSeqsCFs(0,0);
	tmV1.SetDim(0);
	tmV1.Add(0);
	nl = pctSeqsCFs.GetItsDim(ROW);
	for(ni=1;ni<nl;ni++){
		pctSeqsCFs(ni,2) = pctSeqsCFs(ni-1,2)+pctSeqsCFs(ni,0);
		if(pctSeqsCFs(ni,0) >0) tmV1.Add(ni);
		pctSeqsCFs(ni-1,2) /= ntset;
	}
	pctSeqsCFs(ni-1,2) /= ntset;
	// Determine cluster function index close to 10, 20, ---, 100% sequence range
	ctgCFRange.SetDim(10,2); ctgCFRange.SetZero();
	ctgCFList.SetDim(10); ctgCFList.SetZero();
	nl = tmV1.GetItsDim();
	for(ni=0;ni<nl;ni++) {
		nj = tmV1[ni];
		for(mi=0; mi<10; mi++) {
			di = (double)(mi*10.0+10.0);
			if(ni == 0 || (fabs(ctgCFRange(mi,0)-di)>fabs(pctSeqsCFs(nj,2) - di))) {
				ctgCFRange(mi,0) = pctSeqsCFs(nj,2);
				ctgCFList[mi] = nj;
			}
		}
	}
	// Eliminate index which indicates the same cluster function
	mi=1; nl=ctgCFList.GetItsDim();
	for(ni=1;ni<nl;ni++) {
		if(ctgCFList[mi]-ctgCFList[mi-1] ==0) {
			ctgCFList.Remove(mi-1);
			ctgCFRange.Remove(ROW,mi-1);
		}else mi++;
	}
	// Calculate training-set RMSD value of each range
	int nst, nend;
	double sum, count;
	nl = ctgCFList.GetItsDim();
	for(ni=0;ni<nl;ni++) {
		if (ni==0) nst =0;
		else nst = ctgCFList[ni-1]+1;
		nend = ctgCFList[ni];
		sum = 0.0; count =0.0;
		for(mi=nst;mi<=nend;mi++) {
			count += pctSeqsCFs(mi,0);
			sum   += pctSeqsCFs(mi,1);
		}
		if (count == 0.0) ctgCFRange(ni,1) = 0.0;
		else ctgCFRange(ni,1) = sqrt(sum/count);
	}
	if(tmvList.GetItsDim() == 0) {__ERRORTHROW(ERR_ELIMSEQS); }
	// Choose a sequence with the largest prediction error
	CNVEC_INT tmvErrList;
	CNVEC_DBL tmvErr;
	tmvErr = dfInEn.GetSelect(tmvList);
	tmvErr = tmvErr.DirectDot(tmvErr);
	tmvErrList = tmvErr.GetSortIdx(DESCEND);
	tSet.GetSeq(D_ITS,tmvList[tmvErrList[0]],tmvSeq, en); // a sequence with worst error in energy
	
	tSet.sel_BadSeqs.Add(COLUMN,tmvSeq);
	cout << "\t"<< "Selected sequence not fitting the the structure : " << endl;
	cout << "\t"<< tSet.GetSeqText(tC,tmvSeq) << endl;
}

bool CCFS::IsItForbidden(int residue, int aa, CNVEC_INT &seqV)
{
	int ncf, n, j, val, k, add;
	int decor, clust, resaa;
	for(n=0;n<elimCFsList.GetItsDim();n++) {
		ncf = elimCFsList[n];
		CCFDEF & cf = cfDefs[ncf];
		val = 0;
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				decor = cf.lDecor[k];
				clust = cf.lClust[j*cf.nOrd+k];
				if(clust == residue) resaa = aa;
				else resaa =seqV[clust];
				if(resaa != decor) {add = 0; break;}
			}
			val += add;
		}
		if(val!=0) return true;
	}
	return false;
}
bool CCFS::IsItForbidden(CNVEC_INT &seqV)
{
	int ncf, n, j, val, k, add;
	int decor, clust;
	for(n=0;n<elimCFsList.GetItsDim();n++) {
		ncf = elimCFsList[n];
		CCFDEF & cf = cfDefs[ncf];
		val = 0;
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				decor = cf.lDecor[k];
				clust = cf.lClust[j*cf.nOrd+k];
				if(seqV[clust] != decor) {add = 0; break;}
			}
			val += add;
		}
		if(val!=0) return true;
	}
	return false;
}
// obtain system characteristic data
void CCFS::CalcSCD(CNVEC_INT &seq, CNVEC_INT &cfvec) 
{
	int dim = lsSeCf.GetItsDim();
	cfvec.SetDim(dim);
	if(dim == 0) {return;}
	cfvec.SetZero();

	int *rp = cfvec.GetItsVec();

	rp[0] = 1;
	int n, ncf;
	int add,j,k,clust, val;
	int *rvseq = seq.GetItsVec();
	for(n=1;n<dim;n++) {
		ncf = lsSeCf[n];
		CCFDEF &cf = cfDefs[ncf];
		val = 0;
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				clust = cf.lClust[j*cf.nOrd+k];
				if(rvseq[clust] != cf.lDecor[k]) {add = 0; break;}
			}
			val += add;
		}
		rp[n] = val;
	}
}

// obtain system characteristic data
void CCFS::CalcSCD(CNVEC_INT &seq, CNVEC_INT &cfvec, CNVEC_INT &listcfs) 
{
	int dim = listcfs.GetItsDim();
	cfvec.SetDim(dim);
	if(dim == 0) {return;}
	cfvec.SetZero();

	int *rp = cfvec.GetItsVec();

	int n, ncf;
	int add,j,k,clust, val;
	int *rvseq = seq.GetItsVec();
	for(n=0;n<dim;n++) {
		ncf = listcfs[n];
		CCFDEF &cf = cfDefs[ncf];
		val = 0;
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				clust = cf.lClust[j*cf.nOrd+k];
				if(rvseq[clust] != cf.lDecor[k]) {add = 0; break;}
			}
			val += add;
		}
		rp[n] = val;
	}
}

// return system characteristic data for given cluster function(ncf)
void CCFS::CalcSCD(CTSET &tSet, int ncf, CNVEC_INT &cfvec) 
{
	if(ncf >= nTtCf) {__ERRORTHROW(ERR_MEMBDR); }
	int dim = tSet.GetNSeq(D_ITS);
	cfvec.SetDim(dim);
	if(dim == 0) {return;}
	if(ncf == 0) { cfvec = 1; return; } // return constant cluster function.
	cfvec.SetZero();

	int *rp = cfvec.GetItsVec();
	CCFDEF & cf = cfDefs[ncf];
	int *cfdecor = cf.lDecor.GetItsVec();
	int *cfclust = cf.lClust.GetItsVec();
	int add,i,j,k,clust;
	int *rv;
	for(i=0; i<dim; i++) {
		CNVEC_INT &seqV = tSet.GetSeq(D_ITS,i);
		rv = seqV.GetItsVec();
		for(j=0;j<cf.nSym;j++) {
			add=1;
			for(k=0;k<cf.nOrd;k++) {
				clust = cfclust[j*cf.nOrd+k];
				if(rv[clust] != cfdecor[k]) {add = 0; break;}
			}
			rp[i] += add;
		}	
	}
}

// convert sequence information (setSeq) into sequence characteristic information (setSCD)
void CCFS::CalcSCD(CNMAT_INT &setSeq, CNMAT_INT &setSCD) 
{
	int nrow = setSeq.GetItsDim(COLUMN);
	int ncol = lsEci.GetItsDim();
	if(nrow == 0 || ncol == 0) {return;}
	setSCD.SetDim(nrow,1);
	setSCD.Fill(1);
	setSCD.SetDim(nrow,ncol);
	int *rp, *rv;
	int add,j,k,clust, ncf;
	for(int nc=1;nc<ncol;nc++) {
		rp = setSCD.GetVec(COLUMN,nc).GetItsVec();
		ncf = lsSeCf[nc];
		CCFDEF & cf = cfDefs[ncf];
		int *cfdecor = cf.lDecor.GetItsVec();
		int *cfclust = cf.lClust.GetItsVec();
		for(int nr=0; nr<nrow; nr++) {
			CNVEC_INT &seqV = setSeq.GetVec(COLUMN,nr);
			rv = seqV.GetItsVec();
			for(j=0;j<cf.nSym;j++) {
				add=1;
				for(k=0;k<cf.nOrd;k++) {
					clust = cfclust[j*cf.nOrd+k];
					if(rv[clust] != cfdecor[k]) {add = 0; break;}
				}
				rp[nr] += add;
			}	
		}
	}
}

void CCFS::UpdateECIs(CCTRL &tC, CTSET &tSet)
{
	PRINTLOG(<RUN> Update ECI values);
	PRINTMIDREMVAR();
	// Update ECI values by using new traing set and returning rmsd value beween old and new eci values
	int nn = sePm.GetItsDim(ROW);
	CNVEC_INT tmSeq, tmNu;
	CNVEC_DBL  tmT, tmMu, tmvECIs, tmpVec;
	tmvECIs = lsEci;
	double tmD, r_en;
	for(int i = nn; i< tSet.GetNSeq(D_ITS); i++) {
		tSet.GetSeq(D_ITS,i,tmSeq, r_en);
		CalcSCD(tmSeq, tmNu);
		// tmMu = P(n)*nu
		tmMu = eciCo*tmNu;
		// tmD = 1+ nu*P(n)*nu;
		tmD = 1.0 + tmNu*tmMu;
		if(tmD > INVCUTOFF) {
			// P(n+1) = P(n) - P(n)*nu*nu^T*P(n)/d
			eciCo -= tmMu.Dyadic(tmMu)/tmD;
			// b(n+1) = b(n) + (e - nu*b(n))*P(n)*nu/d
			double fac1 = (r_en - tmNu*lsEci)/tmD;
			lsEci += fac1*tmMu;
			// Q(n+1) = [Q(n) - M(n)*P(n)*nu*(M(n)*P(n)*nu)^T/d]_11+[1- 1/d]_22
			tmT = sePm*tmMu; // tmT = M(n)*P(n)*nu;
			dgPj -= tmT.DirectDot(tmT)/tmD;
			dgPj.Add(1.0 - 1.0/tmD);
			// dE(N+1) = [dE(N) - (e - nu*b(n))*tmT/d, (e - nu*b(n))/d]
			dfInEn -= fac1*tmT;
			dfInEn.Add(fac1);
			sePm.Add(ROW,tmNu);
		} else {
			tSet.RemSeq(D_ITS,i);     // eliminate the sequence
			i--;                           // adjust according to reduced sequence
		}
	}

	// cv_rms_p = sqrt(mean((edN./(1-QN)).^2));
	// cv_rms = [cv_rms cv_rms_p]
	tmpVec = dfInEn/(1.0-dgPj);       // direct divide
	tmpVec = tmpVec.DirectDot(tmpVec);
	double cvrms = sqrt(tmpVec.Average());
	crVl.Add(cvrms);
	crVlSmr.Add(cvrms);
}

