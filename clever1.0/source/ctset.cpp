#include <stdio.h>
#include <stdlib.h>
#include "ctset.hpp"
__LOGOUT
EXTERN

ofstream& operator<< (ofstream& fout, CTSET &ctset)
{
	fout << ctset.itsSet;     
	fout << ctset.itsEn;     
	fout << ctset.itsOrd;     
	fout << ctset.addedSet;   
	fout << ctset.addedEn;    
	fout << ctset.addedOrd;    
	fout << ctset.sel_BadSeqs;
	fout << ctset.elimSet;
	fout << ctset.elimEn;
	WRITEDATA(fout, ctset.aN);    
	WRITEDATA(fout, ctset.nAdd);

	fout << ctset.dataSet;   
	fout << ctset.dataEn;    

	return fout;
}

ifstream& operator>> (ifstream& fin, CTSET &ctset)
{
	fin >> ctset.itsSet;     
	fin >> ctset.itsEn;     
	fin >> ctset.itsOrd;     
	fin >> ctset.addedSet;   
	fin >> ctset.addedEn;    
	fin >> ctset.addedOrd;    
	fin >> ctset.sel_BadSeqs;
	fin >> ctset.elimSet;
	fin >> ctset.elimEn;
	READDATA(fin, ctset.aN);    
	READDATA(fin, ctset.nAdd);

	fin >> ctset.dataSet;   
	fin >> ctset.dataEn;    

	return fin;
}

CTSET::CTSET(){}

CSTRING CTSET::GetSeqText(CCTRL &tC, CNVEC_INT &vec)
{
	char seq[MAXLINE];
	CSTRING tmpS;

	strcpy(seq,"");
	for(int k=0;k<vec.GetItsDim();k++) {
		sprintf(seq, "%s %c",seq, tC.dsSt(k,vec[k]) +'A'-1);
	}
	tmpS.SetString(seq);
	if (vec.IsTmp()) delete &vec;
	return tmpS;
}

void CTSET::UpSeqOrdIdx(CNMAT_INT &seqM, CNVEC_INT &ordV)
{
	int nn, dim, pos;
	ordV.SetDim(0);
	dim = seqM.GetItsDim(COLUMN);
	for(nn=0;nn<dim;nn++) {
		CNVEC_INT &seq=seqM[nn];
		pos = FindSeqPos(seqM,ordV,seq);
		if(pos == -1) continue;
		ordV.Insert(pos,nn);
	}
}

int CTSET::FindSeqPos(CNMAT_INT &seqM, CNVEC_INT &ordV, CNVEC_INT &seq)
{
	int min, max, ind, gap;
	int sn;
	int dim = ordV.GetItsDim();
	if (dim == 0) return 0;
	int *ord = ordV.GetItsVec();
	min = 0; max = dim-1;
	gap = max-min;
	if(gap == 0)  {
		sn = seq.CompVec(seqM[ord[min]]);
		if(sn == -1) return min;
		if(sn ==  1) return min+1;
		if(sn == 0) return -1;
	} else{
		sn = seq.CompVec(seqM[ord[min]]);
		if(sn == -1) return min;
		if(sn == 0) return -1;
		sn = seq.CompVec(seqM[ord[max]]);
		if(sn == 1) return max+1;
		if(sn == 0) return -1;
		if (gap == 1) return max;
	} 
	do {
		ind = gap/2 + min;
		sn = seq.CompVec(seqM[ord[ind]]);
		if(sn == -1) max = ind;
		else if (sn == 1)min = ind;
		else return -1;
		gap = max-min;
	} while (gap != 1);
	return max;
}

bool CTSET::CheckAddSeq(DATA_SET set, CNVEC_INT &seq)
{
	int pos;
	switch(set) {
		case D_ITS:
			pos = FindSeqPos(itsSet,itsOrd,seq);
			if(pos == -1) return false;
			itsOrd.Insert(pos,itsSet.GetItsDim(COLUMN));
			itsSet.Add(COLUMN,seq);
			break;
		case D_ADDED:
			pos = FindSeqPos(addedSet,addedOrd,seq);
			if(pos == -1) return false;
			addedOrd.Insert(pos,addedSet.GetItsDim(COLUMN));
			addedSet.Add(COLUMN,seq);
			break;
	}
	return true;
}

// Check validity of a given sequence
bool CTSET::CheckSeq(CNVEC_INT &addV)
{
	int i;
	for(i=0;i<GetNSeq(D_ITS);i++) {
		if(addV == GetSeq(D_ITS,i) ) return false;
	}
	for(i=0;i<GetNSeq(D_ADDED);i++) {
		if(addV == GetSeq(D_ADDED,i)) return false;
	}
	return true;
}

// Calculate probe matrix from a test set.
void CTSET::AddSeq(DATA_SET set, CNVEC_INT &addV) 
{
	switch(set) {
		case D_ITS:
			itsSet.Add(COLUMN,addV);
			break;
		case D_ADDED:
			addedSet.Add(COLUMN,addV);
			break;
		case D_ELIM:
			elimSet.Add(COLUMN,addV);
			break;
		case D_DATA:
			dataSet.Add(COLUMN,addV);
			break;
	}
}

void CTSET::CopySeq(DATA_SET orig, DATA_SET dest) 
{
	CNMAT_INT *tmMat;
	CNVEC_DBL *tmVec;
	switch(orig) {
		case D_ITS:
			tmMat = &itsSet;
			tmVec = &itsEn;
			break;
		case D_ADDED:
			tmMat = &addedSet;
			tmVec = &addedEn;
			break;
		case D_ELIM:
			tmMat = &elimSet;
			tmVec = &elimEn;
			break;
		case D_DATA:
			tmMat = &dataSet;
			tmVec = &dataEn;
			break;
	}
	switch(dest) {
		case D_ITS:
			itsSet = *tmMat;
			itsEn = *tmVec;
			break;
		case D_ADDED:
			addedSet = *tmMat;
			addedEn = *tmVec;
			break;
		case D_ELIM:
			elimSet = *tmMat;
			elimEn = *tmVec;
			break;
		case D_DATA:
			dataSet = *tmMat;
			dataEn = *tmVec;
			break;
	}
}

void CTSET::AddSeq(DATA_SET set, CNVEC_INT &addV, double en) 
{
	switch(set) {
		case D_ITS:
			itsSet.Add(COLUMN,addV);
			itsEn.Add(en);
			break;
		case D_ADDED:
			addedSet.Add(COLUMN,addV);
			addedEn.Add(en);
			break;
		case D_ELIM:
			elimSet.Add(COLUMN,addV);
			elimEn.Add(en);
			break;
		case D_DATA:
			dataSet.Add(COLUMN,addV);
			dataEn.Add(en);
			break;
	}
}

bool CTSET::RemSeq(DATA_SET set, int offset)
{
	bool ok = false;
	switch(set) {
		case D_ITS:
			if(offset < GetNSeq(D_ITS)) {
				itsSet.Remove(COLUMN,offset);
				itsEn.Remove(offset);
				ok = true;
			}
			break;
		case D_ADDED:
			if(offset < GetNSeq(D_ADDED)) {
				addedSet.Remove(COLUMN,offset);
				addedEn.Remove(offset);
				ok = true;
			}
			break;
		case D_ELIM:
			if(offset < GetNSeq(D_ELIM)) {
				elimSet.Remove(COLUMN,offset);
				elimEn.Remove(offset);
				ok = true;
			}
			break;
		case D_DATA:
			if(offset < GetNSeq(D_DATA)) {
				dataSet.Remove(COLUMN,offset);
				dataEn.Remove(offset);
				ok = true;
			}
			break;
	}
	return ok;
}

void CTSET::ClearSeq(DATA_SET set)
{
	switch(set) {
		case D_ITS:
			itsSet.SetDim(0,0);
			itsEn.SetDim(0);
			break;
		case D_ADDED:
			addedSet.SetDim(0,0);
			addedEn.SetDim(0);
			break;
		case D_ELIM:
			elimSet.SetDim(0,0);
			elimEn.SetDim(0);
			break;
		case D_DATA:
			dataSet.SetDim(0,0);
			dataEn.SetDim(0);
			break;
	}
}

void CTSET::RemSeq(DATA_SET set, CNVEC_INT & list)
{
	switch(set) {
		case D_ITS:
			itsSet.Remove(COLUMN,list);
			itsEn.Remove(list);
			break;
		case D_ADDED:
			addedSet.Remove(COLUMN,list);
			addedEn.Remove(list);
			break;
		case D_ELIM:
			elimSet.Remove(COLUMN,list);
			elimEn.Remove(list);
			break;
		case D_DATA:
			dataSet.Remove(COLUMN,list);
			dataEn.Remove(list);
			break;
	}
}

void CTSET::MoveSeq(DATA_SET from_set, CNVEC_INT &list, DATA_SET to_set)
{
	CNVEC_DBL tmV;
	CNMAT_INT tmM;
	switch(from_set) {
		case D_ITS:
			tmM = itsSet.GetSelect(COLUMN,list);
			tmV = itsEn.GetSelect(list);
			itsSet.Remove(COLUMN,list);
			itsEn.Remove(list);
			break;
		case D_ADDED:
			tmM = addedSet.GetSelect(COLUMN,list);
			tmV = addedEn.GetSelect(list);
			addedSet.Remove(COLUMN,list);
			addedEn.Remove(list);
			break;
		case D_ELIM:
			tmM = elimSet.GetSelect(COLUMN,list);
			tmV = elimEn.GetSelect(list);
			elimSet.Remove(COLUMN,list);
			elimEn.Remove(list);
			break;
		case D_DATA:
			tmM = dataSet.GetSelect(COLUMN,list);
			tmV = dataEn.GetSelect(list);
			dataSet.Remove(COLUMN,list);
			dataEn.Remove(list);
			break;
	}
	switch(to_set) {
		case D_ITS:
			itsSet.Add(COLUMN,tmM);
			itsEn.Add(tmV);
			break;
		case D_ADDED:
			addedSet.Add(COLUMN,tmM);
			addedEn.Add(tmV);
			break;
		case D_ELIM:
			elimSet.Add(COLUMN,tmM);
			elimEn.Add(tmV);
			break;
		case D_DATA:
			dataSet.Add(COLUMN,tmM);
			dataEn.Add(tmV);
			break;
	}
}

bool CTSET::GetSeq(DATA_SET set, int offset, CNVEC_INT &seq, double &en)
{
	bool ok = false;
	switch(set) {
		case D_ITS:
			if(offset < GetNSeq(D_ITS)) {
				seq = itsSet.GetVec(COLUMN,offset);
				en  = itsEn[offset];
				ok = true;
			}
			break;
		case D_ADDED:
			if(offset < GetNSeq(D_ADDED)) {
				seq = addedSet.GetVec(COLUMN,offset);
				en  = addedEn[offset];
				ok = true;
			}
			break;
		case D_ELIM:
			if(offset < GetNSeq(D_ELIM)) {
				seq = elimSet.GetVec(COLUMN,offset);
				en  = elimEn[offset];
				ok = true;
			}
			break;
		case D_DATA:
			if(offset < GetNSeq(D_DATA)) {
				seq = dataSet.GetVec(COLUMN,offset);
				en  = dataEn[offset];
				ok = true;
			}
			break;
	}
	return ok;
}

CNVEC_INT & CTSET::GetSeq(DATA_SET set, int offset)
{
	if(set == D_ITS) return itsSet.GetVec(COLUMN,offset);
	else if(set == D_ADDED) return addedSet.GetVec(COLUMN,offset);
	else if(set == D_ELIM) return elimSet.GetVec(COLUMN,offset);
	else return dataSet.GetVec(COLUMN,offset);
}

int CTSET::GetNSeq(DATA_SET set)
{
	if(set == D_ITS) return itsSet.GetItsDim(COLUMN);
	else if(set == D_ADDED) return addedSet.GetItsDim(COLUMN);
	else if(set == D_ELIM) return elimSet.GetItsDim(COLUMN);
	else return dataSet.GetItsDim(COLUMN);
}

void CTSET::PrintSeqs(CCTRL &tC, CCFS &tCFs)
{
	if(!tC.prTrain) return;
	int nn;
	CNVEC_INT tmpV;
	char buff[MAXLINE];
	cout << "Training set" << endl;
	cout << "Num: sequence number" << endl;
	cout << "RealE, E: energies in the training set" << endl;
	cout << "Ep: CE-derived energy" << endl;
	sprintf(buff,"%6s %12s %9s  %s","Num","RealE","E-Ep","Seq");
	cout << buff << endl;
	double r_en;
	for(nn=0;nn<GetNSeq(D_ITS);nn++) {
		GetSeq(D_ITS,nn,tmpV,r_en);
		sprintf(buff,"%6d %12.3f %9.3f  ",nn+1,r_en,tCFs.dfInEn[nn]);
		cout << buff << GetSeqText(tC,tmpV) << endl;
	}

	if(GetNSeq(D_ELIM)==0) return;
	cout << "Excluded Training set" << endl;
	sprintf(buff,"%6s %12s  %s","Num","RealE","Seq");
	cout << buff << endl;
	for(nn=0;nn<GetNSeq(D_ELIM);nn++) {
		GetSeq(D_ELIM,nn,tmpV,r_en);
		sprintf(buff,"%6d %12.3f",nn+1,r_en);
		cout << buff << GetSeqText(tC,tmpV) << endl;
	}
}

void CTSET::PrintAddSeqs(CCTRL &tC)
{
	int nn;
	CNVEC_INT tmpV;
	char buff[MAXLINE];
	double en, r_en;

	for(nn=0;nn<GetNSeq(D_ITS);nn++) {
		GetSeq(D_ITS,nn,tmpV,r_en);
		sprintf(buff,"%12.5f %s",r_en,GetSeqText(tC,tmpV).GetString());
		cout << buff << endl;
	}
	for(nn=0;nn<GetNSeq(D_ADDED);nn++) {
		en = 0.0;
		tmpV = GetSeq(D_ADDED,nn);
		sprintf(buff,"%12.5f %s",en,GetSeqText(tC,tmpV).GetString());
		cout << buff << endl;
	}
}

void CTSET::CalcEnergy(CCTRL &tC, CCFS & tCFs)
{
	PRINTLOG(<RUN> Calculate energy of sequences);
	PRINTMIDREMVAR();
	char exec[MAXLINE];
	double val, en;
	ostrstream *strout;
	streambuf* sbuf;
	if(!tC.prLogEn){
		strout = new ostrstream;
		sbuf = cout.rdbuf();
		cout.rdbuf(strout->rdbuf());
	}

	int totnn = GetNSeq(D_ADDED);
	int toten = addedEn.GetItsDim();
	if(totnn == toten) return;
	CNVEC_INT tmpV;
	CSTRING seq;
	for(int n =toten; n<totnn; n++) {
		// current calculated/total to be calculated/eliminated/total sequence calculated except eliminated sequences
		cout << "\t Running : " << n +1 - toten<< "/" << totnn - toten<< "/" << GetNSeq(D_ELIM);
		cout << "/" << GetNSeq(D_ITS) + addedEn.GetItsDim() +1 <<" ...";
		// Prepare sequence text
		tmpV = GetSeq(D_ADDED,n);
		seq = GetSeqText(tC, tmpV);
		// Calculate energy
		sprintf(exec,"%s %s",tC.fnExec.GetString(), seq.GetString()); 
		val = Execution(exec);
		addedEn.Add(val);

		// Print calculated energy
		if(tCFs.lsSeCf.GetItsDim() !=0 ){
			en = tCFs.CalcCEEnergy(tmpV);
			cout << " done : " << val << " (" << en << ") ";
			APPENDTIME 
			cout << endl;
		}else {
			cout << " done : " << val << " ";
			APPENDTIME 
			cout << endl;
		}
	}
	if(!tC.prLogEn){
		cout.rdbuf(sbuf);
		delete strout;
	}
}

void CTSET::ReCalcEnergy(CCTRL &tC)
{
	PRINTLOG(<RUN> Recalculate energy of sequences);
	PRINTMIDREMVAR();
	char exec[MAXLINE];
	double val, r_en;


	int totnn = GetNSeq(D_ITS);
	CNVEC_INT tmpV;
	CSTRING seq;
	for(int n =0; n<totnn; n++) {
		// current calculated/total to be calculated/eliminated/total sequence calculated except eliminated sequences
		cout << "\t Recalc Running : " << n +1 << "/" << totnn <<" ...";
		// Prepare sequence text
		GetSeq(D_ITS,n,tmpV,r_en);
		seq = GetSeqText(tC, tmpV);
		// Calculate energy
		sprintf(exec,"%s %f %s",tC.fnExec.GetString(), r_en, seq.GetString()); 
		val = Execution(exec);
		itsEn[n] = val;

		// Print calculated energy
		cout << " done : " << val << " ";
		APPENDTIME 
		cout << endl;
	}
}

double CTSET::Execution( const CSTRING & command)
{
	FILE*           fp;
	char            buffer [MAXLINE];
	CSTRING  str;

	fp = NULL;
	if ((fp = POPEN(command.GetString(), "r")) == NULL) {__ERRORTHROW(ERR_POPEN);}
	fgets(buffer, sizeof (buffer), fp);
	str = buffer;
	PCLOSE(fp);

	return str.CvtDbl();
} 

void CTSET::Union() 
{
	if( GetNSeq(D_ADDED) != 0 ) {
		itsSet.Add(COLUMN, addedSet);
		itsEn.Add(addedEn);
		addedEn.SetDim(0);
		addedSet.SetDim(0,0);
	}
}

void CTSET::ReadTrainingSet(FILE *fp, CCTRL &tC)
{
	int k;
	CSTRING str;

	char *tok;
	char cseps[] = " \r\n\t";
	char buffer[MAXLINE];
	CSTRING str2;
	int st, tot;

	// skip file header and find sequence starting column number
	while(!feof(fp)) {
		str.ReadLine(fp);
		tot = str.SplitLine(cseps);
		if(tot < 2) continue;
		if(!str(0).IsNumeric()) continue;
		for(k=0;k<tot;k++) {
			if(!str(k).IsNumeric()) break;
		}
		st = k;
		break;
	}
	int nrow=tot-st;
	if(nrow==0) {__ERRORTHROW(ERR_SEQFORM); }
	int ncol=0;
	do {
		++ncol;
		dataSet.SetDim(nrow, ncol);
		dataEn.SetDim(ncol);
		strcpy(buffer,str.GetString());
		if((tok = strtok(buffer, cseps)) == NULL) {
			dataSet.SetDim(nrow, ncol-1); 
			dataEn.SetDim(ncol-1);
			break;
		}
		str2 = tok;
		dataEn[ncol-1] = str2.CvtDbl();
		// skip columns
		for(k=1;k<st;k++) {
			if((tok = strtok(NULL, cseps)) == NULL) break;
		}
		// read data
		for(k=st;k<tot;k++) {
			if((tok = strtok(NULL, cseps)) == NULL) break;
			int nn = tok[0] - 'A'+1;
			if (nn >= 27 || nn < 1) {
				cout << "Error occured when reading a following line" << endl;
				cout << str << endl;
				__ERRORTHROW(ERR_SEQFORM);
			}
			dataSet(k-st,ncol-1) = tC.rvDsSt(k-st,nn);
		}
		if(nrow != k-st) {
			dataSet.SetDim(nrow, ncol-1); 
			dataEn.SetDim(ncol-1);
			break;
		}
	} while (str.ReadLine(fp)); 
}

