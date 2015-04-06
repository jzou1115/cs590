// MATRIX definition
TEMPLATE_CNMAT
int INST_CNMAT::nMat = 0;

TEMPLATE_CNMAT 
ofstream & operator<< (ofstream& fout, INST_CNMAT &rhs)
{
	CSTRING str = "CNMAT";
	fout << str;
	int row, col;
	row = rhs.GetItsDim(ROW);
	col = rhs.GetItsDim(COLUMN);
	WRITEDATA(fout, row);
	WRITEDATA(fout, col);
	int nc;
	bool ok = rhs.IsTmp();
	if(ok) rhs.SetTmp(false);
	for(nc=0;nc<col;nc++) fout << rhs.GetVec(COLUMN,nc);
	if(ok) delete &rhs;
	return fout;
}

TEMPLATE_CNMAT 
ifstream & operator>> (ifstream& fin,  INST_CNMAT &rhs)
{
	CSTRING str;
	fin >> str;
	if(str != "CNMAT") {__ERRORTHROW(ERR_DATAFORM);}
	int row, col;
	READDATA(fin, row);
	READDATA(fin, col);
	rhs.SetDim(row, col);
	int nc;
	for(nc=0;nc<col;nc++) {
		fin >> rhs.GetVec(COLUMN,nc);
	}
	return fin;
}

TEMPLATE_CNMAT 
ostream & operator<< (ostream& cout, INST_CNMAT &rhs)
{
	int ncol = 1000;
	cout << "ItsNRow: " << rhs.itsNRow << "\n";
	cout << "ItsNCol: " << rhs.itsNCol << "\n";
	int count;
	int row, col;
	row = rhs.GetItsDim(ROW);
	col = rhs.GetItsDim(COLUMN);
	if(col > ncol) {
		for(int i=0; i< row; i++) {
			cout << "Row " << i << "\n";
			count = 0;
			for(int j=0; j< col; j++){ 
				cout << "  " << rhs(i,j);
				count++;
				if(count%ncol == 0) cout << "\n";
			}
			cout << "\n";
		}
		cout << endl;
	} else {
		for(int i=0; i< row; i++) {
			for(int j=0; j< col; j++) 
				cout << "  " << rhs(i,j);
			cout << "\n";
		}
		cout << endl;
	}
	if(rhs.IsTmp()) delete &rhs;
	return cout;
}


TEMPLATE_CNMAT
INST_CNMAT::CNMAT()
{
	++nMat;
	if(nMat == 1) PRINTLOWONE("Matrix created",sizeof(T));
	PRINTHIGHTWO("Matrix",sizeof(T),nMat);
	itsNCol= 0;
	itsNRow= 0;
	Tmp = false;
}

TEMPLATE_CNMAT
INST_CNMAT::CNMAT(int m, int n)
{
	++nMat;
	if(nMat == 1) PRINTLOWONE("Matrix created",sizeof(T));
	PRINTHIGHTWO("Matrix",sizeof(T),nMat);
	itsNCol= 0;
	itsNRow= 0;
	Tmp = false;
	SetDim(m, n);
}

TEMPLATE_CNMAT
INST_CNMAT::CNMAT(INST_CNMAT &rhs)
{
	++nMat;
	if(nMat == 1) PRINTLOWONE("Matrix created",sizeof(T));
	PRINTHIGHTWO("Matrix",sizeof(T),nMat);
	itsNCol= 0;
	itsNRow= 0;
	Tmp = false;
	(*this) = rhs;
	if(rhs.IsTmp()) delete &rhs;
}

TEMPLATE_CNMAT
INST_CNMAT::~CNMAT() 
{
	--nMat;
	PRINTHIGHTWO("Matrix destroyed",sizeof(T),nMat);
	if(nMat == 0) PRINTLOWONE("All matrix destroyed",sizeof(T));
	if(itsNCol != 0) {Clear();}
}

TEMPLATE_CNMAT
int INST_CNMAT::Len(DIRECTION drc) { if(drc) return itsNCol; else return itsNRow; }

TEMPLATE_CNMAT
void INST_CNMAT::Clear()
{
	for(int i=0;i<itsNCol;i++) delete itsMat[i];
	itsMat.SetDim(0);
	itsNRow = itsNCol = 0;
}

TEMPLATE_CNMAT
void INST_CNMAT::SetDim( int m, int n)
{
	int i;

	// if no change in matrix
    if (itsNRow == m && itsNCol == n) { return; }
	// no elements
	if (m == 0 || n == 0) { 
		if(itsNCol != 0) {Clear();}
		return;
	}
	// Adjust number of columns
	if(itsNCol < n) {
		itsMat.SetDim(n);
		for(i=itsNCol;i<n;i++) { 
			itsMat[i] = new INST_CNVEC;
			itsMat[i]->SetDim(itsNRow);
		}
		itsNCol = n;
	} else if (itsNCol > n) {
		for(i=n;i<itsNCol;i++) {
			delete itsMat[i];
			itsMat[i] = NULL;
		}
		itsNCol = n;
	}

	// Adjust number of rows
	if(itsNRow != m) {
		for(i=0;i<itsNCol;i++) { 
			itsMat[i]->SetDim(m);
		}
		itsNRow = m;
	} 
}

TEMPLATE_CNMAT // not complete
INST_CNVEC& INST_CNMAT::GetVec(DIRECTION drc, const int offset)
{
	int nn, n;
	INST_CNVEC* vVec;
	INST_CNVEC* cp;
	T *vec;
	switch(drc) {
	case ROW:
		if(offset>=itsNRow) {__ERRORTHROW(ERR_MEMBDR);}
		nn = GetItsDim(COLUMN);
		if(nn==0) {__ERRORTHROW(ERR_NODIM);}
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsNCol);
		vec = vVec->GetItsVec();
		for(n=0;n<itsNCol;n++) vec[n] = (*itsMat[n])[offset];
		vVec->SetTmp(true);
		cp = vVec;
		break;
	case COLUMN:
		if(offset>=itsNCol) {__ERRORTHROW(ERR_MEMBDR);}
		nn = GetItsDim(ROW);
		if(nn==0) {__ERRORTHROW(ERR_NODIM);}
		if(Tmp) {
			SAFEALLOC(INST_CNVEC,vVec);
			*vVec = *itsMat[offset];
			vVec->SetTmp(true);
			cp = vVec;
		} else {
			cp = itsMat[offset];
		}
		break;
	}
	if(Tmp) delete this;
	return *cp;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::GetSelect(DIRECTION drc, CNVEC_INT & vec)
{
	INST_CNMAT* vMat;
	SAFEALLOC(INST_CNMAT,vMat);
	vMat->SetDim(0,0);
	int i, nn;
	bool ok = Tmp;
	if(ok) Tmp=false;
	for(i=0;i<vec.GetItsDim();i++) {
		nn = vec[i];
		if(nn >= GetItsDim(drc)) {__ERRORTHROW(ERR_MEMBDR);}
		vMat->Add(drc,GetVec(drc,nn));
	}
	if(vec.IsTmp()) delete &vec;
	if(ok) delete this;
	vMat->SetTmp(true);
	return *vMat;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::GetSelect(DIRECTION drc, int st, int end)
{
	CNVEC_INT vs;
	vs.SetDim(0);
	for(int n=st;n<=end;n++) vs.Add(n);
	return GetSelect(drc,vs);
}

// m.Fill(3.0) // Set every elements in m to 3.0
TEMPLATE_CNMAT
void INST_CNMAT::Fill(T val) 
{ 
	for(int i=0;i<itsNCol;i++) itsMat[i]->Fill(val);
}


TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::Tp()
{
	INST_CNMAT* vMat;
	SAFEALLOC(INST_CNMAT,vMat);

	vMat->SetDim(itsNCol,itsNRow);
	int i, j;
	T * pi;
	for(i=0;i<itsNRow;i++) {
		pi = (vMat->GetVec(COLUMN,i)).GetItsVec();
		for(j=0;j<itsNCol;j++) {
			pi[j] = (*itsMat[j])[i];
		}
	}
	if(Tmp) delete this;
	vMat->SetTmp(true);
	return *vMat;
}

// return summation results of each vector
// drc = ROW, return summation for each row vector
TEMPLATE_CNMAT
INST_CNVEC & INST_CNMAT::Sum(DIRECTION drc)
{
	INST_CNVEC* vVec;
	SAFEALLOC(INST_CNVEC,vVec);
	if(itsNCol == 0) {vVec->SetDim(0); return *vVec;}
	vVec->SetZero();
	int nc, nr;
	T *pi, *vv;
	T val;
	switch(drc) {
	case ROW:
		vVec->SetDim(itsNRow);
		vv = vVec->GetItsVec();
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			for(nr=0;nr<itsNRow; nr++) vv[nr] += pi[nr];
		}
		break;
	case COLUMN:
		vVec->SetDim(itsNCol);
		vv = vVec->GetItsVec();
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			val = (T)0.0;
			for(nr=0;nr<itsNRow; nr++) val += pi[nr];
			vv[nc] = val;
		}
		break;
	}
	if(Tmp) delete this;
	vVec->SetTmp(true);
	return *vVec;
}

// return MT*diag*M
TEMPLATE_CNMAT
CNMAT_DBL & INST_CNMAT::SMTransMat(CNVEC_DBL &diag)
{
	if(itsNRow != diag.GetItsDim()) {__ERRORTHROW(ERR_MATMAT); }
	CNMAT_DBL* vDblMat;
	SAFEALLOC(CNMAT_DBL,vDblMat);

	// following function set vector dimension and initialize every elements to be zero when created 
	vDblMat->SetDim(itsNCol,itsNCol);
	double val1;
	double *lamda;
	T *ti, *tk;
	lamda = diag.GetItsVec();
	int i, j ,k;
	for(i=0;i<itsNCol; i++) {
		ti = itsMat[i]->GetItsVec();
		val1 = 0.0;
		for(j=0;j<itsNRow;j++) val1 += lamda[j]*ti[j]*ti[j];
		(*vDblMat)(i,i) = val1;
		for(k=i+1; k<itsNCol; k++) {
			tk = itsMat[k]->GetItsVec();
			val1 = 0.0;
			for(j=0;j<itsNRow;j++) val1 += lamda[j]*ti[j]*tk[j];
			(*vDblMat)(k,i) = (*vDblMat)(i,k) = val1;
		}
	}
	if(diag.IsTmp()) delete &diag;
	vDblMat->SetTmp(true);
	return *vDblMat;
}

// return reference of m, n element
TEMPLATE_CNMAT
T& INST_CNMAT::operator() (const int m, const int n) 
{ 
	if(m >= itsNRow || n >=itsNCol) {__ERRORTHROW(ERR_MEMBDR);}
	return (*itsMat[n])[m]; 
}


TEMPLATE_CNMAT
void INST_CNMAT::Add(DIRECTION drc, CNMAT & mat)
{
	if(mat.GetItsDim(ROW) == 0) {
		TAKECARE(Possible error: No addition);
		return;
	}
	int n;
	if(mat.GetItsDim(COLUMN) == 0) return;
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	switch(drc) {
	case ROW:
		if(itsNCol !=0 && itsNCol != mat.GetItsDim(COLUMN)) {	__ERRORTHROW(ERR_MATDIM);}
		if(itsNCol == 0) {
			*this = mat;
		} else if (itsNCol == mat.GetItsDim(COLUMN)) {
			for(n=0;n<mat.GetItsDim(COLUMN);n++) itsMat[n]->Add(mat.GetVec(COLUMN,n));
			itsNRow = itsNRow + mat.GetItsDim(ROW);
		}
		break;
	case COLUMN:
		if(itsNRow != mat.GetItsDim(ROW) && itsNRow != 0) { __ERRORTHROW(ERR_MATDIM); }
		if(itsNRow == 0) {
			*this = mat;
		} else if (itsNRow == mat.GetItsDim(ROW)) {
			int ncol = itsNCol;
			SetDim(itsNRow,itsNCol+mat.GetItsDim(COLUMN));
			for(n=ncol;n<itsNCol;n++) *itsMat[n] = mat.GetVec(COLUMN,n-ncol);
		}
		break;
	}
	if(ok) delete &mat;
}
TEMPLATE_CNMAT
void INST_CNMAT::Add(DIRECTION drc, INST_CNVEC & vec)
{
	if(vec.GetItsDim() == 0) {
		TAKECARE(Dimension of added vector is zero);
		return;
	}
	int n, mm;
	switch(drc) {
	case ROW:
		if(itsNCol !=0 && itsNCol != vec.GetItsDim()) {	__ERRORTHROW(ERR_MATVEC);}
		if( itsNCol == 0) {
			SetDim(1,vec.GetItsDim());
			for(n=0;n<vec.GetItsDim();n++) (*itsMat[n])[0]=vec[n];
		} else if(itsNCol == vec.GetItsDim()){
			for(n=0;n<vec.GetItsDim();n++) itsMat[n]->Add(vec[n]);
			itsNRow = itsNRow + 1;
		}
		break;
	case COLUMN:
		mm = vec.GetItsDim();
		if(itsNRow != 0 && itsNRow != mm) {__ERRORTHROW(ERR_MATVEC); }
		if(itsNRow == 0) SetDim(vec.GetItsDim(),1);
		else SetDim(itsNRow,itsNCol+1);

		*itsMat[itsNCol-1] = vec;
		break;
	}
	if(vec.IsTmp()) delete &vec;
}

TEMPLATE_CNMAT
void INST_CNMAT::Remove(DIRECTION drc, int offset)
{
	// if there are no data to delete then throw error
	if(offset >= GetItsDim(drc)) {__ERRORTHROW(ERR_MEMBDR);}
	int nc;
	switch(drc) {
	case ROW:
		// if there are only one line to delete then delete the line and return
		if(itsNRow == 1) { this->SetDim(0,0); return; }
		for(nc=0;nc<itsNCol;nc++) itsMat[nc]->Remove(offset);
		--itsNRow;
		break;
	case COLUMN:
		// if there are only one line to delete then delete the line and return
		if(itsNCol == 1) { this->SetDim(0,0); return; }
		delete itsMat[offset];
		itsMat.Remove(offset);
		--itsNCol;
		break;
	}
}

TEMPLATE_CNMAT
void INST_CNMAT::Remove(DIRECTION drc, int offset1, int offset2)
{
	if(offset1 > offset2) {	int nt = offset1; offset1 = offset2; offset2= nt; }
	// if there are no data to delete then throw error
	if(offset2 >= GetItsDim(drc)) {__ERRORTHROW(ERR_MEMBDR);}
	int nc, nr, nn;
	switch(drc) {
	case ROW:
		// if there are only one line to delete then delete the line and return
		nn = itsNRow - (offset2-offset1+1);
		if(nn == 0) { this->SetDim(0,0); return; }
		for(nc=0;nc<itsNCol;nc++) itsMat[nc]->Remove(offset1,offset2);
		itsNRow = nn;
		break;
	case COLUMN:
		// if there are only one line to delete then delete the line and return
		nn = itsNCol - (offset2-offset1+1);
		if(nn == 0) { this->SetDim(0,0); return; }
		for(nr=offset1;nr<=offset2;nr++) delete itsMat[nr];
		itsMat.Remove(offset1, offset2);
		itsNCol = nn;
		break;
	}
}

TEMPLATE_CNMAT
void INST_CNMAT::Remove(DIRECTION drc, CNVEC_INT &listV)
{
	// Sort the elimination list to descending order
	CNVEC_INT tmList;
	tmList = listV.GetSortIdx(DESCEND);
	// If list contains index above boundary then throw error
	if(listV[tmList[0]] >= GetItsDim(drc)) {__ERRORTHROW(ERR_MEMBDR);}
	// Remove data contains same elements
	int nc, nn;
	switch(drc) {
	case ROW:
		// if there are only one line to delete then delete the line and return
		nn =itsNRow - listV.GetItsDim();
		if(nn == 0) { this->SetDim(0,0); return; }
		for(nc=0;nc<itsNCol;nc++) itsMat[nc]->Remove(listV);
		itsNRow = nn;
		break;
	case COLUMN:
		// if there are only one line to delete then delete the line and return
		nn =itsNCol - listV.GetItsDim();
		if(nn == 0) { this->SetDim(0,0); return; }
		for(nc=0;nc<listV.GetItsDim();nc++) delete itsMat[listV[nc]];
		itsMat.Remove(listV);
		itsNCol = nn;
		break;
	}
}

TEMPLATE_CNMAT
void INST_CNMAT::SetZero()
{
	for(int n=0;n<itsNCol;n++) itsMat[n]->SetZero();
}

TEMPLATE_CNMAT
INST_CNVEC & INST_CNMAT::operator[] (int offset)
{
	INST_CNVEC* cp;
	if(offset>=itsNCol) {__ERRORTHROW(ERR_MEMBDR);}
	if(Tmp) {
		SAFEALLOC(INST_CNVEC,cp);
		*cp = *itsMat[offset];
		cp->SetTmp(true);
	} else {
		cp = itsMat[offset];
	}
	if(Tmp) delete this;
	return *cp;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator= (INST_CNMAT &rhs)
{
	if (this == &rhs) return *this;
	SetDim(rhs.GetItsDim(ROW),rhs.GetItsDim(COLUMN));
	bool ok = rhs.IsTmp();
	if(ok) { rhs.SetTmp(false); }
	for(int n=0;n<rhs.GetItsDim(COLUMN);n++) {
		*itsMat[n] = rhs.GetVec(COLUMN,n);
	}
	itsNRow = rhs.GetItsDim(ROW);
	if(ok) delete &rhs;
	return *this;
}
TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator- ()
{
	if(itsNCol == 0) {return *this;}
	INST_CNMAT* vMat;
	int nc, nr;
	T *pi, *pj;
	if(Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = -pi[nr];
		}
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = -pi[nr];
		}
		vMat->SetTmp(true);
		return *vMat;
	}
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator- (const T & val)
{
	if(itsNCol == 0) {return *this;}
	INST_CNMAT* vMat;
	int nc, nr;
	T *pi, *pj;
	if(Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = pi[nr] - val;
		}
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = pi[nr] - val;
		}
		vMat->SetTmp(true);
		return *vMat;
	}
}
TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator- (INST_CNMAT &mat)
{
	if(itsNCol == 0) {return *this;}
	INST_CNMAT* vMat;
	if(itsNCol != mat.GetItsDim(COLUMN) || itsNRow != mat.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATMAT); }
	int nc, nr;
	T *pi, *pj, *pk;
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	if(Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pk = mat.GetVec(COLUMN,nc).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = pi[nr] - pk[nr];
		}
		if(ok) delete &mat;
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			pk = mat.GetVec(COLUMN,nc).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = pi[nr] - pk[nr];
		}
		if(ok) delete &mat;
		vMat->SetTmp(true);
		return *vMat;
	}
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator+ ()
{
	return *this;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator+ (const T & val)
{
	if(itsNCol == 0) {return *this;}
	INST_CNMAT* vMat;
	int nc, nr;
	T *pi, *pj;
	if (Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = pi[nr] + val;
		}
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = pi[nr] + val;
		}
		vMat->SetTmp(true);
		return *vMat;
	}
}
TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator+ (INST_CNMAT &mat)
{
	if(itsNCol == 0) {return *this;}
	if(itsNCol != mat.GetItsDim(COLUMN) || itsNRow != mat.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATMAT); }
	INST_CNMAT* vMat;
	int nc, nr;
	T *pi, *pj, *pk;
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	if(Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pk = mat.GetVec(COLUMN,nc).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = pi[nr] + pk[nr];
		}
		if(ok) delete &mat;
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			pk = mat.GetVec(COLUMN,nc).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = pi[nr] + pk[nr];
		}
		if(ok) delete &mat;
		vMat->SetTmp(true);
		return *vMat;
	}
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator* (INST_CNMAT &mat)
{
	if(itsNCol == 0) {return *this;}
	if(itsNCol != mat.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATMAT); }
	int nrow = itsNRow;
	int ncol = mat.GetItsDim(COLUMN);

	INST_CNMAT* vMat;
	SAFEALLOC(INST_CNMAT,vMat);
	vMat->SetDim(nrow,ncol);
	int ni, nj, nk;
	T *av, *bv, *cv, bval;
	bool ok = mat.IsTmp();
	if(ok)mat.SetTmp(false);
	for(nk=0;nk<ncol; nk++) {
		cv = (vMat->GetVec(COLUMN,nk)).GetItsVec();
		bv = mat.GetVec(COLUMN,nk).GetItsVec();
		for(ni=0;ni<itsNCol; ni++) {
			av = itsMat[ni]->GetItsVec();
			bval = bv[ni];
			for(nj=0;nj<nrow;nj++) {
				cv[nj] += av[nj]*bval;
			}
		}
	}
	if(ok) delete &mat;
	if(Tmp) delete this;
	vMat->SetTmp(true);
	return *vMat;
}

TEMPLATE_CNMAT
INST_CNVEC & INST_CNMAT::operator* (INST_CNVEC &rhs)
{
	if(itsNCol != rhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	INST_CNVEC* vVec;
	SAFEALLOC(INST_CNVEC,vVec);
	if(itsNCol == 0 || rhs.GetItsDim() == 0) {return *vVec;}
	vVec->SetDim(itsNRow);
	vVec->SetZero();
	int ni, nj;
	T *av, *bv, *cv, bval;
	cv = vVec->GetItsVec();
	bv = rhs.GetItsVec();
	for(ni=0;ni<itsNCol; ni++) {
		av = itsMat[ni]->GetItsVec();
		bval = bv[ni];
		for(nj=0;nj<itsNRow;nj++) {
			cv[nj] += av[nj]*bval;
		}
	}
	if(rhs.IsTmp()) delete &rhs;
	if(Tmp) delete this;
	vVec->SetTmp(true);
	return *vVec;
}

// Return double type vector
TEMPLATE_CNMAT
CNVEC_DBL & INST_CNMAT::DotDbl (CNVEC_DBL &rhs)
{
	if(itsNCol != rhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	CNVEC_DBL* vDblVec;
	SAFEALLOC(CNVEC_DBL,vDblVec);
	if(itsNCol == 0 || rhs.GetItsDim() == 0) {return *vDblVec;}
	vDblVec->SetDim(itsNRow);
	int ni, nj;
	T *av;
	double *cv, *bv, bval;
	cv = vDblVec->GetItsVec();
	bv = rhs.GetItsVec();
	for(ni=0;ni<itsNCol; ni++) {
		av = itsMat[ni]->GetItsVec();
		bval = bv[ni];
		for(nj=0;nj<itsNRow;nj++) {
			cv[nj] += ((double)av[nj])*bval;
		}
	}
	if(rhs.IsTmp()) delete &rhs;
	if(Tmp) delete this;
	vDblVec->SetTmp(true);
	return *vDblVec;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator/ (const T rhs)
{
	if(itsNCol == 0) {return *this;}
	if(rhs == (T)0.0) {__ERRORTHROW(ERR_DIVZERO);}
	INST_CNMAT* vMat;
	int nc, nr;
	T *pi, *pj;
	if(Tmp) {
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pi[nr] = pi[nr]/rhs;
		}
		return *this;
	} else {
		SAFEALLOC(INST_CNMAT,vMat);
		vMat->SetDim(itsNRow,itsNCol);
		for(nc=0;nc<itsNCol; nc++) {
			pi = itsMat[nc]->GetItsVec();
			pj = (vMat->GetVec(COLUMN,nc)).GetItsVec();
			for(nr=0;nr<itsNRow; nr++) pj[nr] = pi[nr]/rhs;
		}
		vMat->SetTmp(true);
		return *vMat;
	}
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator+= (T rhs)
{
	if(itsNCol == 0) {SetDim(1,1);}
	int nc, nr;
	T *pi;
	for(nc=0;nc<itsNCol; nc++) {
		pi = itsMat[nc]->GetItsVec();
		for(nr=0;nr<itsNRow; nr++) pi[nr] += rhs;
	}
	return *this;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator+= (INST_CNMAT &rhs)
{
	if(itsNCol != rhs.GetItsDim(COLUMN) || itsNRow != rhs.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATDIM);}
	if(itsNCol == 0 ) {__ERRORTHROW(ERR_NODIM);}
	int nc, nr;
	T *pi, *pj;
	bool ok = rhs.IsTmp();
	if(ok) rhs.SetTmp(false);
	for(nc=0;nc<itsNCol; nc++) {
		pi = itsMat[nc]->GetItsVec();
		pj = (rhs.GetVec(COLUMN,nc)).GetItsVec();
		for(nr=0;nr<itsNRow; nr++) pi[nr] += pj[nr];
	}
	if(ok) delete &rhs;
	return *this;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator-= (T rhs)
{
	if(itsNCol == 0) {SetDim(1,1);}
	int nc, nr;
	T *pi;
	for(nc=0;nc<itsNCol; nc++) {
		pi = itsMat[nc]->GetItsVec();
		for(nr=0;nr<itsNRow; nr++) pi[nr] -= rhs;
	}
	return *this;
}

TEMPLATE_CNMAT
INST_CNMAT & INST_CNMAT::operator-= (INST_CNMAT &rhs)
{
	if(itsNCol != rhs.GetItsDim(COLUMN) || itsNRow != rhs.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATDIM);}
	if(itsNCol == 0 ) {__ERRORTHROW(ERR_NODIM);}
	int nc, nr;
	T *pi, *pj;
	bool ok = rhs.IsTmp();
	if(ok) rhs.SetTmp(false);
	for(nc=0;nc<itsNCol; nc++) {
		pi = itsMat[nc]->GetItsVec();
		pj = (rhs.GetVec(COLUMN,nc)).GetItsVec();
		for(nr=0;nr<itsNRow; nr++) pi[nr] -= pj[nr];
	}
	if(ok) delete &rhs;
	return *this;
}
