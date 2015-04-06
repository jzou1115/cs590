#include <stdio.h>
#include <stdlib.h>
#include "mcsimul.hpp"
__LOGOUT
EXTERN

CMCSIMUL::CMCSIMUL(void)
{
	nCoolStep = 1000;
	srand((unsigned) time (NULL));
}

CMCSIMUL::~CMCSIMUL(void)
{
}

void CMCSIMUL::ReadInput(char *filename, CCTRL &tC)
{
	FILE *fp;

	fp = tUtil.SafeOpen(filename,"r");
		tC.ReadCtrl(fp);
	fclose(fp);
	// Update MC parameters
	UpdateParams(tC);
}

void CMCSIMUL::UpdateParams(CCTRL &tC)
{
	nDS = tC.nDS;
	nAaAtEaDsSt = tC.nAaAtEaDsSt;
	codeN=tC.codeN;        
}

void CMCSIMUL::GenTSet(CCTRL &tC, CCFS &tCFs, CTSET &tSet)
{
	// PRINTLOG(<RUN> Generate sequences);
	// PRINTMIDREMVAR();
	CNVEC_INT tmpV;
	int nn;
	// Generate list of amino acid sequences
	if(tC.numSeqs <= tSet.GetNSeq(D_ADDED)) return;
	nAdd = tC.numSeqs - tSet.GetNSeq(D_ADDED);
	nn = 0;
	while(nn < nAdd) {
		RandSeq(tCFs, tmpV);
		if(tSet.CheckAddSeq(D_ADDED, tmpV)) ++nn;
	}
}

void CMCSIMUL::MCAnneal(CCFS &tCFs, CNVEC_INT &initseq, CNVEC_INT &finalseq)
{
	double en, cur;
	CNVEC_INT seqcV, seqpV;
	double RT;
	RT = 0.591*(dbl_ndT/298);

	seqpV = initseq;
	seqcV = initseq;
	cur = tCFs.CalcCEEnergy(seqpV);
	int ns = seqpV.GetItsDim();
	for(int i=0;i<nCoolStep;i++) {
		ChangeElements(tCFs, seqpV, 1);
		en = tCFs.CalcCEEnergy(seqpV);
		if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
		else { seqpV = seqcV;}
	}
	finalseq = seqcV;
}

void CMCSIMUL::SimAnneal(CCFS &tCFs, CTSET &tSet)
{
	double en, cur;
	int i;
	double RT;
	CNVEC_INT seqcV, seqpV, seqannV;
	RT = 0.591*(dbl_T/298);
	tSet.addedEn.SetDim(0);
	tSet.addedSet.SetDim(0,0);

	RandSeq(tCFs, seqpV);
	seqcV = seqpV;
	cur = tCFs.CalcCEEnergy(seqpV);
	for(i=0;i < nPreStep; i++) { 
		ChangeElements(tCFs, seqpV, 1);
		en = tCFs.CalcCEEnergy(seqpV);
		if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
		else { seqpV = seqcV;}
	}
	bool ok=false;
	for(i=0;; i++) { 
		ChangeElements(tCFs, seqpV, 1);
		en = tCFs.CalcCEEnergy(seqpV);
		if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
		else { seqpV = seqcV;}
		// Save data
		if(i%nIntStep == 0) ok=true;
		if(ok) {
			MCAnneal(tCFs, seqcV, seqannV);
			if(bQ) {if(tCFs.IsItForbidden(seqannV)) continue;}
			if(bSeq) {if(!tSet.CheckAddSeq(D_ADDED, seqannV)) continue;}
			else {
				tSet.AddSeq(D_ADDED, seqannV);
			}
			en = tCFs.CalcCEEnergy(seqannV);
			tSet.addedEn.Add(en);
			if(tSet.addedEn.GetItsDim() >= nData) break;
			ok = false; i=0;
		}
	}
}

void CMCSIMUL::MCSimul(CCFS &tCFs, CTSET &tSet)
{
	double en, cur;
	int i;
	double RT;
	CNVEC_INT seqcV, seqpV;
	RT = 0.591*(dbl_T/298);
	tSet.addedEn.SetDim(0);
	tSet.addedSet.SetDim(0,0);

	RandSeq(tCFs, seqpV);
	seqcV = seqpV;
	cur = tCFs.CalcCEEnergy(seqpV);
	for(i=0;i < nPreStep; i++) { 
		ChangeElements(tCFs, seqpV, 1);
		en = tCFs.CalcCEEnergy(seqpV);
		if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
		else { seqpV = seqcV;}
	}
	bool ok=false;
	for(i=0;; i++) { 
		ChangeElements(tCFs, seqpV, 1);
		en = tCFs.CalcCEEnergy(seqpV);
		if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
		else { seqpV = seqcV;}
		// Save data
		if(i%nIntStep == 0) ok=true;
		if(ok) {
			if(bQ) {if(tCFs.IsItForbidden(seqcV)) continue; }
			if(bSeq) {if(!tSet.CheckAddSeq(D_ADDED, seqcV)) continue;}
			else {
				tSet.AddSeq(D_ADDED, seqcV);
			}
			tSet.addedEn.Add(cur);
			if(tSet.addedEn.GetItsDim() >= nData) break;
			ok = false; i = 0;
		}
	}
}

void CMCSIMUL::MCSpecEn(CCFS &tCFs, double en, CNVEC_INT &seq)
{
	double newen, cur, medgap;
	int i;
	CNVEC_INT seqcV, seqpV, seqMed;

	RandSeq(tCFs, seqpV);
	seqcV = seqMed =seqpV;
	cur = tCFs.CalcCEEnergy(seqpV);
	medgap = fabs(cur - en);
	for(i=0;i < nCoolStep; i++) { 
		ChangeElements(tCFs, seqpV, 1);
		newen = tCFs.CalcCEEnergy(seqpV);
		if(cur >= en) {
			if (newen < cur) { cur = newen; seqcV = seqpV; }
			else { seqpV = seqcV;}
		} else {
			if (newen > cur) { cur = newen; seqcV = seqpV;}
			else { seqpV = seqcV;}
		}
		if(medgap > fabs(cur - en)) {
			medgap = fabs(cur - en);
			seqMed = seqcV;
		}
	}
	seq = seqMed;
}

void CMCSIMUL::MCLowMedHigh(CCFS &tCFs, CTSET &tSet)
{
	CNVEC_INT tmvSeq;
	double enlow, enhigh, enmed;
	// Add a sequence with lower temperature
	do {
		RandSeq(tCFs, tmvSeq);
		MCAnneal(tCFs,tmvSeq,tmvSeq);
		if(tSet.CheckAddSeq(D_ADDED,tmvSeq)) break;
	} while (1);
	enlow = tCFs.CalcCEEnergy(tmvSeq);

	// Add a random sequence preventing systematic convergence with bias
	do {
		RandSeq(tCFs,tmvSeq);
		if(tCFs.IsItForbidden(tmvSeq)) continue;
		if(tSet.CheckAddSeq(D_ADDED,tmvSeq)) break;
	} while (1);
	enhigh = tCFs.CalcCEEnergy(tmvSeq);

	// Add a sequence with intermediate energy between lowest and random energy 
	enmed = enlow + (enhigh - enlow)/2.0;
	do {
		MCSpecEn(tCFs, enmed, tmvSeq);
		if(tSet.CheckAddSeq(D_ADDED,tmvSeq)) break;
	} while (1);
}

// Temperature scan
void CMCSIMUL::MCTempScan(CCFS &tCFs, CTSET &tSet, CNVEC_DBL &tList, CNMAT_INT &listM)
{
	double en, cur;
	int i, nn, it;
	double RT;
	CNVEC_INT seqcV, seqpV;
	CNVEC_DBL tVec;
	tVec = 0.591*tList/298.0;

	tSet.addedEn.SetDim(0);
	tSet.addedSet.SetDim(0,0);
	int dim = tVec.GetItsDim();

	listM.SetDim(nData,dim);
	listM.SetZero();
	for(nn=0;nn<nData;nn++) {
		RandSeq(tCFs,seqpV);
		seqcV = seqpV;
		cur = tCFs.CalcCEEnergy(seqpV);
		RT = tVec[0];
		for(i=0;i < nPreStep; i++) { 
			ChangeElements(tCFs, seqpV, 1);
			en = tCFs.CalcCEEnergy(seqpV);
			if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
			else { seqpV = seqcV;}
		}
		for(it=0;it<dim;it++) {
			RT = tVec[it];
			bool ok=false;
			for(i=1;; i++) { 
				ChangeElements(tCFs, seqpV, 1);
				en = tCFs.CalcCEEnergy(seqpV);
				if(Rand2() < exp(-(en-cur)/RT)) { cur = en; seqcV = seqpV; }
				else { seqpV = seqcV;}
				// Save data
				if(i%nIntStep == 0) ok=true;
				if(ok) {
					if(bQ) {if(tCFs.IsItForbidden(seqcV)) continue; }
					if(bSeq) {if(!tSet.CheckAddSeq(D_ADDED, seqcV)) continue;}
					else {
						tSet.AddSeq(D_ADDED, seqcV);
					}
					tSet.addedEn.Add(cur);
					break;
				}
			}
			listM(nn,it)=tSet.addedEn.GetItsDim()-1;
		}
	}
}

// Temperature scan
void CMCSIMUL::GenRandSeqs(CCFS &tCFs, CTSET &tSet)
{
	int nn;
	CNVEC_INT seqV;

	for(nn=0;nn<nData;nn++) {
		do {
			RandSeq(tCFs,seqV);
			if(bQ) {if(tCFs.IsItForbidden(seqV)) continue; }
			if(bSeq) {if(tSet.CheckAddSeq(D_ITS, seqV)) break;}
			else {
				tSet.AddSeq(D_ITS, seqV);
				break;
			}
		} while (1);
	}
}

// evaluate random sequence without forbidden interaction
void CMCSIMUL::RandSeq(CCFS &tCFs, CNVEC_INT &seq)
{
	int i, j;
	int nsite, nseq, indseq, indsite;

	CNVEC_INT tmpV, tmpSite, tmpSeq;
	tmpV.SetDim(nDS);
	tmpSite.SetDim(nDS);
	tmpSite.FillInc();
	for(i=0;i<nDS;i++) {
		indsite = (int)((nDS-i)*Rand2());
		nsite = tmpSite[indsite];
		int na = nAaAtEaDsSt[nsite];
		tmpSeq.SetDim(na);
		tmpSeq.FillInc();
		for(j=0;j<na;j++) {
			indseq = (int)((na-j)*Rand2());
			nseq = tmpSeq[indseq];
			if(!tCFs.IsItForbidden(nsite,nseq,tmpV)) {
				tmpV[nsite] = nseq;
				break;
			}
			tmpSeq[indseq] = tmpSeq[na-j-1];
		}
		tmpSite[indsite] = tmpSite[nDS-i-1];
	}
	seq = tmpV;
}

void CMCSIMUL::ChangeElements(CCFS &tCFs, CNVEC_INT &vec, int num)
{
	int i, j, indseq, nseq;
	if(num == 0) return;
	i=0;
	CNVEC_INT tmpSeq;
	while(1) {
		int n = (int)(nDS*Rand2());
		int oldele = vec[n];
		int na = nAaAtEaDsSt[n];
		tmpSeq.SetDim(na);
		tmpSeq.FillInc();
		for(j=0;j<na;j++) {
			indseq = (int)((na-j)*Rand2());
			nseq = tmpSeq[indseq];
			if(!tCFs.IsItForbidden(n,nseq,vec)) {
				vec[n] = nseq;
				break;
			}
			tmpSeq[indseq] = tmpSeq[na-j-1];
		}
		if(vec[n] == oldele && j != na-1) continue;
		i++;
		if(i >= num) break;
	}
}

double CMCSIMUL::Rand2()
{
	return (double)rand()/(RAND_MAX+1.0);

}
