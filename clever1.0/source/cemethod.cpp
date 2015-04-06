/*
	cemethod.cpp
*/
#include "cemethod.hpp"
EXTERN

CAPPCEM::CAPPCEM() {
	nECIStep = 0;
	sStep = 0;
	cStep = 0;
	nStep = 0;
	nStep_cs = 0;
	itsNum = 0;
	itsCount = 0;
	strVer = VERSION;
}

void CAPPCEM::PrintHeader()
{
	cout << "###############################################################" << endl;
	cout << "  CLEVER can be used without restriction for academic purposes." << endl;
	cout << "  http://web.mit.edu/biology/keating                           " << endl;
	cout << "  version: 1.0 (Feb 5, 2010)                                   " << endl;
	cout << "  correspondence:                                              " << endl;
	cout << "     Amy Keating (keating@mit.edu)                             " << endl;
	cout << "  Major Author:                             " << endl;
	cout << "     Seungsoo Hahn (hahnss@kaist.ac.kr, coding source and writing docs)" << endl;
	cout << "  Minor Author:                             " << endl;
	cout << "     Orr Ashenberg (writing docs and comment)                  " << endl;
	cout << "  Other contributor:                             " << endl;
	cout << "     Gevorg Grigoryan (comment)                                " << endl;
	cout << "###############################################################\n" << endl;
}

void CAPPCEM::CalcCE_File()
{
	FILE *fp;
	switch(nStep) {
		case  1: cStep = nStep_cs; goto START_POS_1;
		case  2: cStep = nStep_cs; goto START_POS_2;
		case  3: return;
		case 10: sStep = cStep = nStep_cs; goto START_POS_CEI;
	}
	// Read sequence file (when mode = 1)
	if (sStep == 0) {
		fp = tUtil.SafeOpen(tC.fnSeqFile.GetString(),"r");
		tSet.ReadTrainingSet(fp,tC);
		fclose(fp);
		if(tC.numSeqs > tSet.GetNSeq(D_DATA)) {
			cout << "WARNING: Number of sequences in the data file is less than input number" << endl;
			tC.numSeqs = tSet.GetNSeq(D_DATA);
		}
		tSet.addedEn = tSet.dataEn.GetSelect(0,tC.numSeqs-1);
		tSet.addedSet = tSet.dataSet.GetSelect(COLUMN,0,tC.numSeqs-1);
		tSet.itsEn.SetDim(0);
		tSet.itsSet.SetDim(0,0);
	}
	// Initialize Cluster function information
	RUN_THIS(1, tCFs.GenInitCFs(tC);)
nStep = 1; 	nStep_cs = cStep;
Save();
START_POS_1:
	RUN_THIS(0, tSet.Union();)     // append additional sequences into training set.
	do {
		RUN_THIS(2, tCFs.UpdateCFs(tC, tSet);)
		RUN_THIS(1, tCFs.CalcRemCFsECIs(tC,tSet);)
//	RUN_THIS(1, tCFs.AnalSeqErr(tC,tSet);) 
START_POS_CEI:
		if(!tCFs.ChooseErrInteract(tC,tSet)) break; 
nStep = 2; 	nStep_cs = cStep;
Save();
START_POS_2:
		RUN_THIS(1, tCFs.ElimSeqs(tC,tSet);) 
	} while(1);
nStep = 3; 	nStep_cs = cStep;
Save();
}

void CAPPCEM::CalcCE_Rand()
{
	switch(nStep) {
		case  1: cStep = nStep_cs; goto START_POS_1;
		case  2: cStep = nStep_cs; goto START_POS_2;
		case  3: return;
		case 10: sStep = cStep = nStep_cs; goto START_POS_CEI;
	}
	// Initialize Cluster function information
	RUN_THIS(1, tCFs.GenInitCFs(tC);)
	// Generate random sequences and calculate energy
	RUN_THIS(1, tMC.GenTSet(tC, tCFs, tSet);) 
	RUN_THIS(1, tSet.CalcEnergy(tC, tCFs);)  //calculate energy of additional sequences.
nStep = 1; 	nStep_cs = cStep;
Save();
START_POS_1:
	RUN_THIS(0, tSet.Union();)     // append additional sequences into training set.
	do {
		RUN_THIS(2, tCFs.UpdateCFs(tC, tSet);)
		RUN_THIS(1, tCFs.CalcRemCFsECIs(tC,tSet);)
START_POS_CEI:
		if(!tCFs.ChooseErrInteract(tC,tSet)) break; 
nStep = 2; 	nStep_cs = cStep;
Save();
START_POS_2:
		RUN_THIS(1, tCFs.ElimSeqs(tC,tSet);) 
	} while(1);
nStep = 3; 	nStep_cs = cStep;
Save();
}

void CAPPCEM::GenRandSeqs()
{
	// Initialize Cluster function information
	tCFs.GenInitCFs(tC);
	// Generate random sequences and calculate energy
	tMC.GenTSet(tC, tCFs, tSet); 
}

void CAPPCEM::Save()
{
	Save(1);
}

void CAPPCEM::Save(int np)
{
	if(np == 0) return;
	ofstream fout;
	if(np == 1) {
		fout.open(tC.fnOutPutDump.GetString(),std::ios::binary);
	} else if(np == 2 && tC.saveOption) {
		char filename[MAXLINE];
		sprintf(filename,"%s_%04d",tC.fnOutPutDump.GetString(),tCFs.nUpCFs);
		fout.open(filename,std::ios::binary);
	} else {
		Save(--np);
		return;
	}
	fout << strVer;
	fout << tC;
	WRITEDATA(fout, sStep);
	WRITEDATA(fout, nUpCFsOk);
	WRITEDATA(fout, rmsEci);
	WRITEDATA(fout, countN);
	WRITEDATA(fout, nStep);
	WRITEDATA(fout, nStep_cs);
	WRITEDATA(fout, itsNum);
	WRITEDATA(fout, itsCount);
	WRITEDATA(fout, nECIStep);
	fout << tCFs;
	fout << tSet;
	fout.close();
	Save(--np);
}

void CAPPCEM::SaveTrainResults()
{
	Save();
	// Eliminate some informations which are not used in energy prediction
	tCFs.sePm.SetDim(0,0);
	tCFs.pmFoAlCf.SetDim(0,0);
	tCFs.svRs.SetDim(0,0);
	tCFs.crVl.SetDim(0);
	tSet.addedEn.SetDim(0);
	tSet.addedSet.SetDim(0,0);
	tSet.elimEn.SetDim(0);
	tSet.elimSet.SetDim(0,0);
	tSet.dataEn.SetDim(0);
	tSet.dataSet.SetDim(0,0);
	tSet.itsEn.SetDim(0);
	tSet.itsSet.SetDim(0,0);
	CCFLIST cflst;
	int ncf, n;
	for(n=0;n<tCFs.lsSeCf.GetItsDim();n++) {
		ncf = tCFs.lsSeCf[n];
		cflst.Add(tCFs.cfDefs[ncf]);
	}
	for(n=0;n<tCFs.elimCFsList.GetItsDim();n++) {
		ncf = tCFs.elimCFsList[n];
		cflst.Add(tCFs.cfDefs[ncf]);
	}
	tCFs.cfDefs = cflst;

	ofstream fout;
	fout.open(tC.fnEciFile.GetString(),std::ios::binary);

	fout << strVer;
	fout << tC;
	WRITEDATA(fout, sStep);
	WRITEDATA(fout, nUpCFsOk);
	WRITEDATA(fout, rmsEci);
	WRITEDATA(fout, countN);
	WRITEDATA(fout, nStep);
	WRITEDATA(fout, nStep_cs);
	WRITEDATA(fout, itsNum);
	WRITEDATA(fout, itsCount);
	WRITEDATA(fout, nECIStep);
	fout << tCFs;
	fout << tSet;
	fout.close();
	Read();
}

bool CAPPCEM::Read()
{
	if(!tC.reDo) return false;
	// temperary variables for compatibility of different versions
	// Check whether a saved file exists
	ifstream fin(tC.fnOutPutDump.GetString(),std::ios::in|std::ios::binary);
	if(!fin.is_open()) return false;

	fin >> rFileVer;
	if(!rFileVer.IsStrIn("ver")) {
		return false;
//		fin.seekg(0,std::ios::beg);
	}
	tC.SetVer(rFileVer);
	tCFs.SetVer(rFileVer);
	tSet.SetVer(rFileVer);
	CCTRL         tc_old;          // Old tC data
	tc_old = tC;
	fin >> tC;
	READDATA(fin, sStep);
	READDATA(fin, nUpCFsOk);
	READDATA(fin, rmsEci);
	READDATA(fin, countN);
	READDATA(fin, nStep);
	READDATA(fin, nStep_cs);
	READDATA(fin, itsNum);
	READDATA(fin, itsCount);
	READDATA(fin, nECIStep);
	fin >> tCFs;
	fin >> tSet;

	int change = 0;
	if (tc_old.nMode != tC.nMode ) change = 1;
	tC.nMode = tc_old.nMode;

	fin.close();
	cStep = 0;

	tC.SetVer(strVer);
	tCFs.SetVer(strVer);
	tSet.SetVer(strVer);
	// save control data to proper varables.
	CNVEC_INT tmpList;
	for(int i=0;i < MAXINPUTITEM; i++) {
		if(!tC.lstInput[i].mFlag) continue;
		if (tC.lstInput[i].Compare(tc_old.lstInput[i])) {
			if(tC.lstInput[i].nMode==0) change = 1;
			else change = 2;
		}
	}
	if(change == 1) {
		nStep_cs = 0;
		cStep = 0;
		nStep = 0;
		sStep = 0;
		nECIStep = 0;
		itsNum = 0;
		itsCount = 0;
		nUpCFsOk = 1;
		tCFs.pnVl.SetDim(0);
	} else if (change ==2) {
		if (tC.nMode == 3 || tC.nMode == 2) {
			if (nStep >=6 ) nStep = 10;
		} else if (tC.nMode == 1 || tC.nMode == 4) {
			if (sStep >= 3) sStep = 3;
			if (nStep >=3 ) nStep = 10;
		}
	}
	tMC.UpdateParams(tC);
	return true;
}

void CAPPCEM::PrintRandSeqs() 
{
	tSet.PrintAddSeqs(tC);
}

void CAPPCEM::PrintResults()
{
	int nn, mm, mm1;
	double cvrms, dblpr;
	char buff[MAXLINE];
	int nsave;

	cout << "\n### Results ###" << endl;
	// print  overal information
	sprintf(buff,"%4s %9s %9s %9s %9s %9s %9s","step","nSeq","nCFs","CVRMS",
		"Max(NHR)","DOF(f)","RMSD");
	cout << buff << endl;
	for(nn=0;nn<tCFs.svRsSmr.GetItsDim(ROW);nn++) {
		sprintf(buff,"%4.0f %9.0f %9.0f %9.3f %9.5f %9.3f %9.3f",tCFs.svRsSmr(nn,0),
			tCFs.svRsSmr(nn,1),tCFs.svRsSmr(nn,2),tCFs.svRsSmr(nn,3),
			tCFs.svRsSmr(nn,6),tCFs.svRsSmr(nn,9),tCFs.svRsSmr(nn,13));
		cout << buff << endl;
	}
	cout << endl;
	// print not fitted sequence to the given backbone
	for(nn=0;nn<tCFs.elimCFsList.GetItsDim();nn++) {
		if(nn==0) {cout << "Bad sequences" << endl;}
		mm = tCFs.elimCFsList[nn];
		cout << tCFs.cfDefs[mm].GetCFText(tC) << "  ";
		cout << tSet.GetSeqText(tC,tSet.sel_BadSeqs.GetVec(COLUMN,nn)) << endl;
	}
	if(tCFs.elimCFsList.GetItsDim() !=0) { cout << endl; }
	// Some additional calculations
	CNVEC_INT tmpV;
	CSTRING tmpStr;
	tmpV = tCFs.lsSeCf.GetSortIdx(ASCEND);
	// print selected cluster functions
	cvrms = tCFs.crVlSmr.GetEndVal();
//	sprintf(buff,"\nCVRMS value : %12.5f",cvrms);
//	cout << buff << endl;
//	sprintf(buff,"Mag Sq. Err : %12.5f",tCFs.magErr);
//	cout << buff << endl;
	cout << "Selected cluster functions" << endl;
	cout << "1/IPR: participation number (representing the number of sequences independently" << endl;
	cout << "       involved in the CF in a training set)" << endl;
	sprintf(buff,"%5s %12s %9s %9s %s","Num","ECI","NH-ratio","1/IPR","CFInfo");
	cout << buff << endl;

	CNVEC_DBL tmpPR;
	CNVEC_INT tmplCFs;
	if(tCFs.svRs.GetItsDim(ROW) != 0) {
		nsave   = tCFs.svRs.GetItsDim(COLUMN)-2;
		tmpPR   = tCFs.svRs.GetVec(COLUMN,nsave);
		nsave   = tCFs.svRsIdx.GetItsDim(COLUMN) -1;
		tmplCFs = tCFs.svRsIdx.GetVec(COLUMN,nsave);
		for(nn=0;nn<tmpV.GetItsDim();nn++) {
			int nn1 = tmpV[nn];
			mm = tCFs.lsSeCf[nn1];
			tmpStr = tCFs.cfDefs[mm].GetCFText(tC);
			mm1 = tmplCFs.FindElem(mm);
			dblpr = tmpPR[mm1];
			double zval = tCFs.htEciCo(nn1,nn1)/(cvrms*cvrms*tCFs.eciCo(nn1,nn1));
			double cosval = (zval-1.0)*sqrt(dblpr)/tCFs.magErr;
			sprintf(buff,"%5d %12.6f %9.5f %9.2f%s",nn+1,tCFs.lsEci[nn1], cosval, dblpr, tmpStr.GetString());
			cout << buff << endl;
		}
	}
	cout << "\n### Definition and other information ###" << endl;
	// print ECI values of other CFs
/*	if(tC.prRemCFs && tCFs.svRs.GetItsDim(ROW)!=0) {
		cout << "Remaining CFs" << endl;
		sprintf(buff,"%5s %12s %9s %9s %9s %s","Num","ECI value","DelAlpha","1/IPR","Cos(X)","CFInfo");
		cout << buff << endl;
		double teci, tda, tpr, tcos;
		int ns = tCFs.svRs.GetItsDim(COLUMN)-4;
		int ns1 = tCFs.svRsIdx.GetItsDim(COLUMN)-1;
		mm=0;
		for(nn=0;nn<tCFs.svRs.GetItsDim(ROW);nn++) {
			nn1 = tCFs.svRsIdx(nn,ns1);
			if(tCFs.lsSeCf.FindElem(nn1) != -1) continue; // skip existent CFs. 
			mm++;
			teci = tCFs.svRs(nn,ns);
			tda  = tCFs.svRs(nn,ns+1);
			tpr  = tCFs.svRs(nn,ns+2);
			tcos = tCFs.svRs(nn,ns+3);
			tmpStr = tCFs.cfDefs[nn1].GetCFText(tC);
			sprintf(buff,"%5d %12.6f %9.3f %9.3f %9.3f %s",nn1,teci, tda, tpr, tcos, tmpStr.GetString());
			cout << buff << endl;
		}
	}
*/
	// print definition of cluster function 
	tCFs.cfDefs.PrintDefCFs(tC);

	// Print cv score
	if(tC.prCVRMS) {
		cout << "\nCVRMS change" << endl;
		cout << "INum: iteration number" << endl;
		sprintf(buff,"%6s %9s","INum", "CVRMS");
		cout << buff << endl;
		for(nn=0;nn<tCFs.crVl.GetItsDim();nn++) {
			sprintf(buff,"%6d %9.3f",nn+1, tCFs.crVl[nn]);
			cout << buff << endl;
		}
	}
	// print training sequences 
	cout << endl;
	tSet.PrintSeqs(tC, tCFs);
	cout << endl;
	

/*	// print ECI value fluctuation
	if(tC.prCEDetail && tCFs.svRs.GetItsDim(ROW)!=0) {
		cout << "ECI variance" << endl;
		char name0[20], name1[20], name2[20], name3[20], name4[20];
		double val1, val2, val3, val4;
		int val0;
		int nentry = 3;
		if(tCFs.svRs.GetItsDim(COLUMN) !=0) {
			int num = 0;
			int numst = (tCFs.svRs.GetItsDim(COLUMN)-1)/(nentry*4);
			for(nn=0; nn <=numst; nn++) {
				strcpy(buff,"");
				for(mm = 1; mm<=nentry; mm++){
					if(mm+nn*nentry > tCFs.svRs.GetItsDim(COLUMN)/4) continue; 
					sprintf(name0,"nCF%d",mm+nn*nentry);
					sprintf(name1,"ECI%d",mm+nn*nentry);
					sprintf(name2,"DAlp%d",mm+nn*nentry);
					sprintf(name3,"PR%d",mm+nn*nentry);
					sprintf(name4,"Cos%d",mm+nn*nentry);
					sprintf(buff,"%s %6s %9s %9s %9s %9s",buff,name0,name1,name2,name3,name4);
				}
				cout << buff << endl;
				for(mm=0;mm<tCFs.svRs.GetItsDim(ROW);mm++) {
					strcpy(buff,"");
					for(int mm1 = 1; mm1 <=nentry; mm1++) {
						if(mm1+nn*nentry > tCFs.svRs.GetItsDim(COLUMN)/4) continue;
						int mm2 = nn*nentry*4 + (mm1-1)*4;
						val0 = tCFs.svRsIdx(mm,mm1+nn*nentry-1);
						val1 = tCFs.svRs(mm,mm2);
						val2 = tCFs.svRs(mm,mm2+1);
						val3 = tCFs.svRs(mm,mm2+2);
						val4 = tCFs.svRs(mm,mm2+3);
						sprintf(buff,"%s %6d %9.3f %9.3f %9.1f %9.5f",buff,val0,val1,val2,val3,val4);
					}
					cout << buff << endl;
				}
				cout << endl;
			}
		}
	} //if(tC.prCEDetail)*/
}

void CAPPCEM::PrintTestResults()
{
	// print log
	int nn;
	CNVEC_INT tmpV;
	char buff[MAXLINE];
	cout << "Test set" << endl;
	double r_en;
	for(nn=0;nn<tSet.GetNSeq(D_ADDED);nn++) {
		tSet.GetSeq(D_ADDED,nn, tmpV, r_en);
		sprintf(buff,"%6d %12.3f  ",nn+1,r_en);
		cout << buff << tSet.GetSeqText(tC,tmpV) << endl;
	}
	tCFs.cfDefs.PrintDefCFs(tC);
}

void CAPPCEM::ReadInput(CSTRING fname)
{
	FILE *fp;
	fp = tUtil.SafeOpen(fname.GetString(),"r");
	tC.ReadCtrl(fp);
	fclose(fp);

	fp = tUtil.SafeOpen(tC.fnDSites.GetString(),"r");
	tC.ReadDesign(fp);
	fclose(fp);

	// Update MC parameters
	tMC.UpdateParams(tC);
}

void CAPPCEM::ReadDesign(CSTRING fname)
{
	FILE *fp;
	fp = tUtil.SafeOpen(fname.GetString(),"r");
		tC.ReadDesign(fp);
	fclose(fp);
	// Update MC parameters
	tMC.UpdateParams(tC);
}

void CAPPCEM::TrainResults()
{
	tCFs.CalcEnDiff(tSet.itsEn,tSet.itsSet);
}


void CAPPCEM::PrintCEFromSeq()
{
	FILE *fp;
	int nn;
	// Read sequence file (when mode = 1)
	fp = tUtil.SafeOpen(tC.fnSeqFile.GetString(),"r");
	tSet.ReadTrainingSet(fp,tC);
	fclose(fp);

	char buff[MAXLINE];
	double cvrms = tCFs.crVlSmr.GetEndVal();

	sprintf(buff,"Degrees of Free : %12d",tSet.GetNSeq(D_ITS)-tCFs.lsSeCf.GetItsDim());
	cout << buff << endl;
	sprintf(buff,"CVRMS value     : %12.5f",cvrms);
	cout << buff << endl;

	sprintf(buff,"%12s %12s %12s %12s %12s %12s","Ener","CEEner","E-Ep", 
		"HeteroSig","ScaledEnDiff","Seq");
	cout << buff << endl;
	tSet.itsEn = tSet.dataEn;
	tSet.itsSet = tSet.dataSet;

	nn = tSet.GetNSeq(D_ITS);
	CNVEC_INT tmSeq;
	double diff, gammahetero, en, heterosig, r_en;
	CNVEC_INT tmNu;
	for(int i=0;i<nn;i++) {
		tSet.GetSeq(D_ITS,i,tmSeq, r_en);
		en = tCFs.CalcCEEnergy(tmSeq);
		tCFs.CalcSCD(tmSeq,tmNu);
		gammahetero = tmNu*tCFs.htEciCo*tmNu;
		heterosig = sqrt(cvrms*cvrms + gammahetero);
		diff = tSet.itsEn[i] - en;
		sprintf(buff,"%12.4f %12.4f %12.4f %12.4f %12.4f  ",r_en,en,diff,
			heterosig,diff/heterosig);
		cout << buff << tSet.GetSeqText(tC,tmSeq) << endl; 
	}
}

void CAPPCEM::PrintECIs(CCFS &tCFs)
{
	ofstream fout("tmpvECIvalues.dat");
	fout << "dblvECIs\n" << endl;
	int i, j;
	int nrow, ncol;
	nrow = tCFs.svRs.GetItsDim(ROW);
	ncol = tCFs.svRs.GetItsDim(COLUMN);
	for(i=0;i<ncol;i++) {
		for(j=0;j<nrow;j++) {
			fout << "     " << tCFs.svRs(j,i);
		}
		fout << endl;
	}
	fout.close();
}
