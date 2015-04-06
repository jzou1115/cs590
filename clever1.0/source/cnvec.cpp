TEMPLATE_CNVEC
int INST_CNVEC::nVec = 0;

/******************** Global Functions *******************/
TEMPLATE_CNVEC 
ofstream & operator<< (ofstream& fout, INST_CNVEC &rhs)
{
	CSTRING str;
	str = "CNVEC";
	fout << str;
	int dim = rhs.GetItsDim();
	WRITEDATA(fout, dim);
	if(dim != 0){
		T *rp = rhs.GetItsVec();
		WRITEARRDATA(fout, rp, dim);
	}
	if(rhs.IsTmp()) delete &rhs;
	return fout;
}

TEMPLATE_CNVEC 
ifstream & operator>> (ifstream& fin,  INST_CNVEC &rhs)
{
	CSTRING str;
	fin >> str;
	if(str != "CNVEC") {__ERRORTHROW(ERR_DATAFORM);}
	int dim;
	READDATA(fin, dim);
	rhs.SetDim(dim);
	if(dim != 0){
		T *rp = rhs.GetItsVec();
		READARRDATA(fin, rp, dim);
	}
	return fin;
}

TEMPLATE_CNVEC 
ostream & operator<< (ostream& cout, INST_CNVEC &rhs)
{
	int dim = rhs.GetItsDim();
	cout << "Dim " << dim << endl;
	for(int i=0; i< dim; i++) 
		cout << i << "  " << rhs[i] << endl;
	if (rhs.IsTmp()) delete &rhs;
	return cout;
}

TEMPLATE_CNVEC 
istream & operator>> (istream& cin,  INST_CNVEC &rhs)
{
	int j, dim;
	char tmpstr[30];
	cin >> tmpstr >> dim;
	rhs.SetDim(dim);
	for(int i=0; i< dim; i++) 
		cin >> j >> rhs[i];
	return cin;
}
TEMPLATE_CNVEC 
INST_CNVEC & operator+ (T lhs,  INST_CNVEC &rhs)
{
	int n;
	int dim = rhs.GetItsDim();
	if(dim == 0) return rhs;
	T * rv = rhs.GetItsVec();
	if (rhs.IsTmp()){
		for(n=0;n<dim;n++) rv[n] = lhs + rv[n];
		return rhs;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(dim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<dim;n++) rp[n] = lhs + rv[n];
		vVec->SetTmp(true);
		return *vVec;
	}
}
TEMPLATE_CNVEC 
INST_CNVEC & operator- (T lhs,  INST_CNVEC &rhs)
{
	int n;
	int dim = rhs.GetItsDim();
	if(dim == 0) return rhs;
	T * rv = rhs.GetItsVec();
	if (rhs.IsTmp()){
		for(n=0;n<dim;n++) rv[n] = lhs - rv[n];
		return rhs;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(dim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<dim;n++) rp[n] = lhs - rv[n];
		vVec->SetTmp(true);
		return *vVec;
	}
}
TEMPLATE_CNVEC 
INST_CNVEC & operator* (T lhs,  INST_CNVEC &rhs)
{
	int n;
	int dim = rhs.GetItsDim();
	if(dim == 0) return rhs;
	T * rv = rhs.GetItsVec();
	if (rhs.IsTmp()){
		for(n=0;n<dim;n++) rv[n] *= lhs;
		return rhs;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec)
		vVec->SetDim(dim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<dim;n++) rp[n] = rv[n]*lhs;
		vVec->SetTmp(true);
		return *vVec;
	}
}


/******************* Member Functions  *******************/
TEMPLATE_CNVEC
INST_CNVEC::CNVEC()
{
	++nVec;
	if(nVec == 1) PRINTLOWONE("Vector created",sizeof(T));
	PRINTHIGHTWO("Vector",sizeof(T),nVec);
	itsDim = 0;
	maxDim = 0;
	itsVec= NULL;
	Tmp = false;
}

TEMPLATE_CNVEC
INST_CNVEC::CNVEC(int n)
{
	++nVec;
	if(nVec == 1) PRINTLOWONE("Vector created",sizeof(T));
	PRINTHIGHTWO("Vector",sizeof(T),nVec);
	maxDim = 0;
	itsVec= NULL;
	itsDim = 0;
	Tmp = false;
	SetDim(n);
}

TEMPLATE_CNVEC
INST_CNVEC::CNVEC(INST_CNVEC &rhs)
{
	++nVec;
	if(nVec == 1) PRINTLOWONE("Vector created",sizeof(T));
	PRINTHIGHTWO("Vector",sizeof(T),nVec);
	maxDim = 0;
	itsVec= NULL;
	itsDim = 0;
	Tmp = false;
	(*this) = rhs;
}

TEMPLATE_CNVEC
INST_CNVEC::~CNVEC() 
{
	--nVec;
	PRINTHIGHTWO("Vector destroyed",sizeof(T),nVec);
	if(nVec == 0) PRINTLOWONE("All vector destroyed",sizeof(T));
	if(itsVec != NULL) delete[] itsVec;
}

// adjust dimension of vector
TEMPLATE_CNVEC
void INST_CNVEC::SetDim( int n) {
	if(n<=maxDim) {
		if (n < itsDim) memset(itsVec+n,0,sizeof(T)*(itsDim-n));
	} else {
		T* tmpV;
		maxDim = n+MEM_EXPAND;
		tmpV = new T[maxDim];
		memset(tmpV,0,sizeof(T)*maxDim);
		if(itsVec != NULL) {
			memcpy(tmpV,itsVec,sizeof(T)*itsDim);
			delete[] itsVec;
		}
		itsVec = tmpV;
	}
	itsDim = n;
}

// no additional space
TEMPLATE_CNVEC
void INST_CNVEC::SetMaxDim( int n) {
	if (n == itsDim && n == maxDim) return;       // no change
	if (n == 0 && itsDim > 0 ) {   // set dimension to zero
		delete[] itsVec;
		itsVec = NULL;
		itsDim = maxDim = 0;
		return;
	}
	T* tmpV;
	tmpV = new T[n];
	memset(tmpV,0,sizeof(T)*n);
	if(itsDim > 0) {
		if(n > itsDim) memcpy(tmpV,itsVec,sizeof(T)*itsDim);
		else memcpy(tmpV,itsVec,sizeof(T)*n);
		delete[] itsVec;
	}
	itsVec = tmpV;
	itsDim = maxDim = n;
}

TEMPLATE_CNVEC
void INST_CNVEC::Reduce() {
	SetMaxDim(itsDim);
}

TEMPLATE_CNVEC
void INST_CNVEC::SetZero()
{
	memset(itsVec,0,sizeof(T)*itsDim);
}

// v1.Fill(3.0) // Set every elements in v1 to 3.0
TEMPLATE_CNVEC
void INST_CNVEC::Fill(T val) 
{ 
	for(int i=0;i<itsDim;i++) itsVec[i] = val;
}

// v1.Fill(3.0, n, m) // Set elements between "n" and "m" in v1 to 3.0
TEMPLATE_CNVEC
void INST_CNVEC::Fill( T val, int n, int m) { 
	if(n > m) {SWAP_T(int,n,m);}
	for(int i=n;i<=m;i++) {
		if(i>=itsDim){__ERRORTHROW(ERR_MEMBDR);}
		itsVec[i] = val;
	}
}

// v1.Fill(3.0, v) // Set list of elements to 3.0
TEMPLATE_CNVEC
void INST_CNVEC::Fill( T val, CNVEC_INT & vec) { 
	for(int i=0;i<vec.GetItsDim();i++) {
		int n1 = vec[i];
		if(n1>=itsDim){__ERRORTHROW(ERR_MEMBDR);}
		itsVec[n1] = val;
	}
	if(vec.IsTmp()) delete &vec;
}

// v1.FillInc() // Set elements from zero to dimension
TEMPLATE_CNVEC
void INST_CNVEC::FillInc() {
	for(int i=0;i<itsDim;i++) itsVec[i] = (T)i;
}

TEMPLATE_CNVEC
void INST_CNVEC::FillInc(T val) {
	for(int i=0;i<itsDim;i++) itsVec[i] = (T)i + val;
}

// v1.DirectDot(v2) = v1.*v2 in matlab
TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::DirectDot(CNVEC &rhs)
{
	int n;
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM);}
	if (Tmp || rhs.IsTmp()){
		T * rp = rhs.GetItsVec();
		if(Tmp) { 
			for(n=0;n<itsDim;n++) itsVec[n] = itsVec[n] * rp[n];
			if(rhs.IsTmp() && this != &rhs) delete &rhs;
			return *this;
		} else {
			for(n=0;n<itsDim;n++) rp[n] = itsVec[n] * rp[n];
			return rhs;
		}
	} else {
		INST_CNVEC *tmpVec;
		SAFEALLOC(INST_CNVEC,tmpVec);
		tmpVec->SetDim(rhs.GetItsDim());
		T *rp = rhs.GetItsVec();
		T *sp = tmpVec->GetItsVec();
		for(n=0;n<itsDim;n++) { sp[n] = itsVec[n] * rp[n]; }
		tmpVec->SetTmp(true);
		return *tmpVec;
	}
}

TEMPLATE_CNVEC
INST_CNMAT & INST_CNVEC::Dyadic (INST_CNVEC &rhs)
{
	if(itsDim == 0 || rhs.GetItsDim() == 0) {__ERRORTHROW(ERR_NODIM); }
	INST_CNMAT * tmpMat;
	SAFEALLOC(INST_CNMAT,tmpMat);
	// following function set vector dimension and initialize every elements to be zero when created 
	tmpMat->SetDim(itsDim, rhs.GetItsDim());
	T *rp, *rv;
	int i,j;
	rv = rhs.GetItsVec();
	for(i=0;i<rhs.GetItsDim(); i++) {
		rp = (tmpMat->GetVec(COLUMN,i)).GetItsVec();
		for(j=0; j<itsDim; j++) {
			rp[j] = itsVec[j]*rv[i];
		}
	}
	if(rhs.IsTmp() && this != &rhs) delete &rhs;
	if(Tmp) delete this;
	tmpMat->SetTmp(true);
	return *tmpMat;
}

TEMPLATE_CNVEC
double INST_CNVEC::Average()
{
	int n;
	T val;
	val = (T)0.0;
	if(itsDim == 0) return val;
	for(n=0;n<itsDim;n++) { val += itsVec[n]; }
	if(Tmp) delete this;
	return (double)val/(double)itsDim;
}

TEMPLATE_CNVEC
T INST_CNVEC::Sum ()
{
	int n;
	T val=(T)0.0;
	for(n=0;n<itsDim;n++) val += itsVec[n];
	return val;
}

TEMPLATE_CNVEC
int INST_CNVEC::CountElemExcept(T val)
{
	int n, count;
	if(itsDim == 0) return 0;
	count=0;
	for(n=0;n<itsDim;n++) { if(itsVec[n] != val) ++count; }
	return count;
}

// Copy vector rhs to itsVec, two vectors must have same dimension
TEMPLATE_CNVEC
bool INST_CNVEC::Copy(INST_CNVEC &rhs)
{
	if (this == &rhs) return true;
	int dim = rhs.GetItsDim();
	if( dim == 0) {
		SetDim(0);
		if(rhs.IsTmp()) delete &rhs;
		return true;
	}
	T *rp = rhs.GetItsVec();
	if(rhs.IsTmp()) {
		if(itsVec != NULL) delete itsVec;
		itsVec = rp;
		rhs.SetItsVec(NULL);
		itsDim = dim;
		maxDim = rhs.GetItsMaxDim();
		delete &rhs;
	} else {
		SetDim(dim);
		memcpy(itsVec,rp,sizeof(T)*itsDim);
	}
	return true;
}

// Copy part of vector from offset1 to offset2
TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::GetSelect(int offset1, int offset2)
{
	// set offset2 have larger number
	if(offset1 > offset2) { int tmp = offset1; offset2 = offset1; offset1=tmp; }
	if(offset2 >= itsDim) {__ERRORTHROW(ERR_MEMBDR); }
	INST_CNVEC *tmpV;
	SAFEALLOC(INST_CNVEC,tmpV)
	if(offset1 == offset2) {
		tmpV->Add(itsVec[offset1]);
	} else {
		int dim = offset2-offset1+1;
		tmpV->SetDim(dim);
		T *rp = tmpV->GetItsVec();
		memcpy(rp,itsVec+offset1,sizeof(T)*dim);
	}
	if(Tmp) delete this;
	tmpV->SetTmp(true);
	return *tmpV;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::GetSelect(CNVEC_INT & vec)
{
	INST_CNVEC* tmpV;
	SAFEALLOC(INST_CNVEC, tmpV);
	int dim = vec.GetItsDim();
	if(dim != 0) {
		tmpV->SetDim(dim);
		int *rv = vec.GetItsVec();
		T *rp = tmpV->GetItsVec();
		int i, nn;
		for(i=0;i<dim;i++) {
			nn = rv[i];
			if(nn >= GetItsDim()) {__ERRORTHROW(ERR_MEMBDR);}
			rp[i] = itsVec[nn];
		}
	}
	if(vec.IsTmp()) delete &vec;
	if(Tmp) delete this;
	tmpV->SetTmp(true);
	return *tmpV;
}

// Find an element that is matched with given value "val"
// if matched then return offset of the value in the vector, if not then return -1;
TEMPLATE_CNVEC
int INST_CNVEC::FindElem(T val)
{
	int offset = -1;
	for(int i=0;i<itsDim;i++) {
		if(itsVec[i] == val) { offset = i; break; }
	}
	if(Tmp) delete this;
	return offset;
}

//TEMPLATE_CNVEC
//INST_CNVEC & INST_CNVEC::DotUINT(CNMAT_INT &mat)
//{
//	if(mat.GetItsDim(ROW) != itsDim) {__ERRORTHROW(ERR_MATVEC); }
//	INST_CNVEC * tmpVec;
//	SAFEALLOC(INST_CNVEC,tmpVec);
//	int dim = mat.GetItsDim(COLUMN);
//	tmpVec->SetDim(dim);
//	T val1;
//	int i, j;
//	int *rv;
//	T * rp = tmpVec->GetItsVec();
//	bool ok = mat.IsTmp();
//	if(ok) mat.SetTmp(false);
//	for(i=0;i<dim; i++) {
//		rv = (mat.GetVec(COLUMN,i)).GetItsVec();
//		val1 = (T)0.0;
//		for(j=0; j<itsDim; j++) {
//			val1 += itsVec[j]*((T)rv[j]);
//		}
//		rp[i] = val1;
//	}
//	if( ok) delete &mat;
//	if(Tmp) delete this;
//	tmpVec->SetTmp(true);
//	return *tmpVec;
//}

TEMPLATE_CNVEC
CNVEC_DBL & INST_CNVEC::UnitVec()
{
	if(itsDim == 0) {__ERRORTHROW(ERR_NOELEM); }
	CNVEC_DBL * tmpVec;
	SAFEALLOC(CNVEC_DBL,tmpVec);
	tmpVec->SetDim(itsDim);
	int i;
	double val=0.0;
	for(i=0;i<itsDim; i++) {
		val += (double)(itsVec[i]*itsVec[i]);
	}
	val = sqrt(val);
	if(val == 0.0) {__ERRORTHROW(ERR_DIVZERO); }
	double *rp = tmpVec->GetItsVec();
	for(i=0;i<itsDim; i++) {
		rp[i] = (double)(itsVec[i])/val;
	}
	if(Tmp) delete this;
	tmpVec->SetTmp(true);
	return *tmpVec;
}

TEMPLATE_CNVEC
double INST_CNVEC::Magnitude()
{
	double val=0.0;
	for(int i=0;i<itsDim; i++) {
		val += (double)(itsVec[i]*itsVec[i]);
	}
	if(Tmp) delete this;
	return sqrt(val);
}

TEMPLATE_CNVEC
T INST_CNVEC::Max()
{
	if(itsDim == 0) return -999.0;
	T val=itsVec[0];
	for(int i=1;i<itsDim; i++) {
		if(val < itsVec[i]) val = itsVec[i];
	}
	if(Tmp) delete this;
	return val;
}

TEMPLATE_CNVEC
double INST_CNVEC::Skewness()
{
	if(itsDim < 3) return 0.0;
	double val, val1, val2, val3, skew, denom;
	val1 = val2 = val3 = 0.0;
	for(int i=0;i<itsDim; i++) {
		val = (double)(itsVec[i]);
		val1 += val;
		val2 += val*val;
		val3 += val*val*val;
	}
	val1 /= (double)itsDim;
	val2 /= (double)itsDim;
	val3 /= (double)itsDim;
	denom = (val2 - val1*val1)*sqrt(val2 - val1*val1);
	if(denom == 0.0) return 0.0;
	skew = (val3 + 2.0*val1*val1*val1 - 3.0*val1*val2)/denom;
	if(Tmp) delete this;
	return skew;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator+ (CNVEC &rhs)
{
	int n;
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0 ) return *this;
	INST_CNVEC *vVec;
	T *rv = rhs.GetItsVec();
	if (Tmp || rhs.IsTmp()){
		if(Tmp) {
			for(n=0;n<itsDim;n++) { itsVec[n] += rv[n]; }
			if(rhs.IsTmp() && this != &rhs) delete &rhs;
			return *this;
		} else {
			for(n=0;n<itsDim;n++) { rv[n] += itsVec[n]; }
			return rhs;
		}
	}else {
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) rp[n] = itsVec[n]+rv[n];
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator+ (const T &rhs)
{
	int n;
	if (Tmp){
		for(n=0;n<itsDim;n++) 
			itsVec[n] += rhs;
		return *this;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) {
			rp[n] = itsVec[n] + rhs;
		}
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator- (CNVEC &rhs)
{
	int n;
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0 ) return *this;
	INST_CNVEC *vVec;
	T *rv = rhs.GetItsVec();
	if (Tmp || rhs.IsTmp()){
		if(Tmp) {
			for(n=0;n<itsDim;n++) { itsVec[n] -= rv[n]; }
			if(rhs.IsTmp() && this != &rhs) delete &rhs;
			return *this;
		} else {
			for(n=0;n<itsDim;n++) { rv[n] = itsVec[n] - rv[n]; }
			return rhs;
		}
	}else {
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) rp[n] = itsVec[n]-rv[n];
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator- (const T rhs)
{
	int n;
	if (Tmp){
		for(n=0;n<itsDim;n++) 
			itsVec[n] -= rhs;
		return *this;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) {
			rp[n] = itsVec[n] - rhs;
		}
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator- ()
{
	int n;
	if (Tmp){
		for(n=0;n<itsDim;n++) itsVec[n] = -itsVec[n];
		return *this;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) {
			rp[n] = -itsVec[n];
		}
		vVec->SetTmp(true);
		return *vVec;
	}
}


TEMPLATE_CNVEC
T INST_CNVEC::operator* (CNVEC &rhs)
{
	int n;
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	T val = (T)0.0;
	if(itsDim == 0) return val;
	T *rv = rhs.GetItsVec();
	for(n=0;n<itsDim;n++) val += itsVec[n]*rv[n];
	if(rhs.IsTmp() && this != &rhs) delete &rhs;
	if(Tmp) delete this;
	return val;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator* (const T &rhs)
{
	int n;
	if (Tmp){
		for(n=0;n<itsDim;n++) itsVec[n] *= rhs;
		return *this;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) rp[n] = itsVec[n]*rhs;
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator* (INST_CNMAT &rhs)
{
	int row = rhs.GetItsDim(ROW);
	int col = rhs.GetItsDim(COLUMN);
	if(itsDim != row) {__ERRORTHROW(ERR_MATVEC); }
	INST_CNVEC * vVec;
	SAFEALLOC(INST_CNVEC,vVec);
	if(itsDim == 0) return *vVec;
	int i, j;
	T val, *rv;
	vVec->SetDim(col);
	T *rp = vVec->GetItsVec();
	for(j=0;j<col; j++) {
		rv = (rhs.GetVec(COLUMN,j)).GetItsVec();
		val = (T)0.0;
		for(i=0; i<row; i++) val += rv[i]*itsVec[i];
		rp[j] = val;
	}
	if(rhs.IsTmp()) delete &rhs;
	if(Tmp) delete this;
	vVec->SetTmp(true);
	return *vVec;
}

// Direct divide
TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator/ (INST_CNVEC &rhs)
{
	int n;
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0 ) return *this;
	INST_CNVEC *vVec;
	T *rv = rhs.GetItsVec();
	if (Tmp || rhs.IsTmp()){
		if(Tmp) {
			for(n=0;n<itsDim;n++) { 
				if(rv[n]==(T)0.0) {__ERRORTHROW(ERR_DIVZERO);}
				itsVec[n] /= rv[n]; 
			}
			if(rhs.IsTmp() && this != &rhs) delete &rhs;
			return *this;
		} else {
			for(n=0;n<itsDim;n++) { 
				if(rv[n]==(T)0.0) {__ERRORTHROW(ERR_DIVZERO);}
				rv[n] = itsVec[n]/rv[n]; 
			}
			return rhs;
		}
	}else {
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) { 
			if(rv[n]==(T)0.0) {__ERRORTHROW(ERR_DIVZERO);}
			rp[n] = itsVec[n]/rv[n]; 
		}
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator/ (const T &rhs)
{
	int n;
	if(rhs == (T)0.0) {__ERRORTHROW(ERR_DIVZERO);}
	if (Tmp){
		for(n=0;n<itsDim;n++) itsVec[n] /= rhs;
		return *this;
	}else {
		INST_CNVEC *vVec;
		SAFEALLOC(INST_CNVEC,vVec);
		vVec->SetDim(itsDim);
		T *rp = vVec->GetItsVec();
		for(n=0;n<itsDim;n++) rp[n] = itsVec[n]/rhs;
		vVec->SetTmp(true);
		return *vVec;
	}
}

TEMPLATE_CNVEC
CNVEC_DBL & INST_CNVEC::CvtDbl()
{
	int n;
	CNVEC_DBL *vDblVec;
	SAFEALLOC(CNVEC_DBL,vDblVec);
	vDblVec->SetDim(itsDim);
	double *rp = vDblVec->GetItsVec();
	for(n=0;n<itsDim;n++) rp[n] = (double) itsVec[n];
	if(Tmp) delete this; 
	vDblVec->SetTmp(true);
	return *vDblVec;
}

TEMPLATE_CNVEC
CNVEC_INT & INST_CNVEC::CvtUINT()
{
	int n;
	CNVEC_INT *vVec;
	SAFEALLOC(CNVEC_INT,vVec);
	vVec->SetDim(itsDim);
	int *rp = vVec->GetItsVec();
	for(n=0;n<itsDim;n++) rp[n] = (int) itsVec[n];
	if(Tmp) delete this; 
	vVec->SetTmp(true);
	return *vVec;
}

TEMPLATE_CNVEC
int INST_CNVEC::CompVec(INST_CNVEC &vec)
{
	if(this == &vec) {
		if(Tmp) delete this;
		return 0;
	}
	if(itsDim != vec.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0) {
		if(vec.IsTmp()) delete &vec;
		if(Tmp) delete this;
		return 0;
	}
	int i;
	T * rv = vec.GetItsVec();
	int val=0;
	for(i=0;i<itsDim;i++) {
		if(itsVec[i] != rv[i]) {
			if(itsVec[i] < rv[i]) val = -1;
			else val = 1;
			break;
		}
	}
	if(vec.IsTmp()) delete &vec;
	if(Tmp) delete this;
	return val;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator+= (INST_CNVEC &rhs)
{
	if (itsDim != rhs.GetItsDim()){__ERRORTHROW(ERR_VECDIM); }
	T * rv = rhs.GetItsVec();
	for(int i=0; i<itsDim;i++) itsVec[i] += rv[i];
	if(rhs.IsTmp() && this != &rhs)	delete &rhs;
	return *this;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator+= (const T rhs)
{
	if (itsDim == 0){SetDim(1);}
	for(int i=0; i<itsDim;i++) itsVec[i] += rhs;
	return *this;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator-= (INST_CNVEC &rhs)
{
	if (itsDim != rhs.GetItsDim()){__ERRORTHROW(ERR_VECDIM); }
	T * rv = rhs.GetItsVec();
	for(int i=0; i<itsDim;i++) itsVec[i] -= rv[i];
	if(rhs.IsTmp() && this != &rhs)	delete &rhs;
	return *this;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator-= (const T rhs)
{
	if (itsDim == 0){SetDim(1);}
	for(int i=0; i<itsDim;i++) itsVec[i] -= rhs;
	return *this;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator= (INST_CNVEC &rhs)
{
	if (this == &rhs) return *this;
	if (Tmp) {
		delete this;
		return rhs;
	}
	Copy(rhs);
	return *this;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::operator= (T rhs)
{
	if(itsDim == 0)	SetDim(1);
	for(int n=0; n<itsDim;n++) itsVec[n] = rhs;
	return *this;
}

TEMPLATE_CNVEC
CNVEC_INT & INST_CNVEC::operator|| (CNVEC_INT &rhs)
{
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0) return *this;
	int i;
	T * rv = rhs.GetItsVec();
	if (rhs.IsTmp() || Tmp) {
		if (Tmp) {
			for(i=0;i<itsDim;i++) itsVec[i] = itsVec[i] || rv[i];
			if(rhs.IsTmp() && this != &rhs)	delete &rhs;
			return *this;
		}else {
			for(i=0;i<itsDim;i++) rv[i] = itsVec[i] || rv[i];
			return rhs;
		}
	} else {
		CNVEC_INT *tmpV;
		SAFEALLOC(CNVEC_INT, tmpV);
		tmpV->SetDim(itsDim);
		int *rp = tmpV->GetItsVec();
		for(i=0;i<itsDim;i++) rp[i] = itsVec[i] || rv[i];
		tmpV->SetTmp(true);
		return *tmpV;
	}
}

TEMPLATE_CNVEC
CNVEC_INT & INST_CNVEC::operator&& (CNVEC_INT &rhs)
{
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0) return *this;
	int i;
	T * rv = rhs.GetItsVec();
	if (rhs.IsTmp() || Tmp) {
		if (Tmp) {
			for(i=0;i<itsDim;i++) itsVec[i] = itsVec[i] && rv[i];
			if(rhs.IsTmp() && this != &rhs)	delete &rhs;
			return *this;
		}else {
			for(i=0;i<itsDim;i++) rv[i] = itsVec[i] && rv[i];
			return rhs;
		}
	} else {
		CNVEC_INT *tmpV;
		SAFEALLOC(CNVEC_INT, tmpV);
		tmpV->SetDim(itsDim);
		int *rp = tmpV->GetItsVec();
		for(i=0;i<itsDim;i++) rp[i] = itsVec[i] && rv[i];
		tmpV->SetTmp(true);
		return *tmpV;
	}
}

TEMPLATE_CNVEC
bool INST_CNVEC::operator== (INST_CNVEC &rhs)
{
	if(this == &rhs) {
		if(Tmp) delete this;
		return true;
	}
	if(itsDim != rhs.GetItsDim()) {__ERRORTHROW(ERR_VECDIM); }
	if(itsDim == 0) return true;
	bool ok = true;
	int i;
	T * rv = rhs.GetItsVec();
	for(i=0;i<itsDim;i++) {
		if(itsVec[i] != rv[i]) {
			ok = false;
			break;
		}
	}
	if(rhs.IsTmp() && this != &rhs) delete &rhs;
	if(Tmp) delete this;
	return ok;
}

TEMPLATE_CNVEC
T& INST_CNVEC::operator[] (const int n)
{
	if(n>=itsDim) {__ERRORTHROW(ERR_MEMBDR); }
	return *(itsVec + n);
}

TEMPLATE_CNVEC
void INST_CNVEC::Print() 
{
	fprintf(stdout,"\nitsDim: %4d\n",itsDim);
	for(int i=0; i< itsDim; i++) 
		cout << i << "  " << itsVec[i] << "\n";
	fprintf(stdout,"\n");
}

TEMPLATE_CNVEC
void INST_CNVEC::Write(ofstream &fout) 
{
	fout << itsVec[0];
}

TEMPLATE_CNVEC
void INST_CNVEC::Read(ifstream &fin)
{
	fin >> itsVec[0];
}

TEMPLATE_CNVEC
void INST_CNVEC::Add(INST_CNVEC &vec)
{
	if (Tmp) {
		if (vec.IsTmp() && this != &vec) delete &vec;
		delete this; return;
	}
	int dim = vec.GetItsDim();
	if(dim == 0) { 
		if (vec.IsTmp() && this != &vec) delete &vec;
		return;
	}
	int tlen = itsDim+dim;
	int cdim = itsDim;
	if( tlen > maxDim) SetDim(tlen);
	else itsDim = tlen;
	memcpy(itsVec+cdim,vec.GetItsVec(),sizeof(T)*dim);

	if(vec.IsTmp()) delete &vec;
}

TEMPLATE_CNVEC
void INST_CNVEC::Add(T val)
{
	if (Tmp) {delete this; return;}
	int tlen = itsDim+1;
	int cdim = itsDim;
	if( tlen > maxDim) SetDim(tlen);
	else itsDim = tlen;
	itsVec[cdim] = val;
}

TEMPLATE_CNVEC
void INST_CNVEC::Insert(int offset, T val)
{
	if (Tmp) {delete this; return;}
	if(offset >= itsDim) {
		if(offset < maxDim) {itsDim = offset+1; itsVec[offset] = val;}
		else {
			SetDim(offset+1);
			itsVec[offset] = val;
		}
	} else {
		if(itsDim == maxDim) SetDim(itsDim+1);
		else ++itsDim;
		memmove(itsVec+offset+1,itsVec+offset,sizeof(T)*(itsDim-1-offset));
		itsVec[offset] = val;
	}
}

TEMPLATE_CNVEC
void INST_CNVEC::Remove(int offset)
{
	if (Tmp) {delete this; return;}
	// if there are no data to delete then throw error
	if(offset >= itsDim) {__ERRORTHROW(ERR_MEMBDR);}
	if(itsDim == 1) { this->SetDim(0); return; }

	// copy block below the deleted row vector
	if(offset != itsDim-1) {
		memmove(itsVec+offset,itsVec+(offset+1),sizeof(T)*(itsDim-offset-1));
	}
	SetDim(itsDim-1);
}

TEMPLATE_CNVEC
void INST_CNVEC::Remove(int offset1, int offset2)
{
	if (Tmp) {delete this; return;}
	if(offset1 > offset2) {	SWAP_T(int,offset1,offset2); } // nt = offset1; offset1 = offset2; offset2= nt; }
	// if there are no data to delete then throw error
	if(offset2 >= itsDim) {__ERRORTHROW(ERR_MEMBDR);}
	// if there are only one point to delete then delete the point and return

	int nn = (itsDim-(offset2-offset1+1));
	if(nn == 0) { this->SetDim(0); return; }
	// copy block below the deleted row vector
	if(offset2 != itsDim-1) {
		memmove(itsVec+offset1,itsVec+(offset2+1),sizeof(T)*(itsDim-offset2-1));
	}
	SetDim(nn);
}

TEMPLATE_CNVEC
void INST_CNVEC::Remove(CNVEC_INT &listV)
{
	if (Tmp) {
		if(listV.IsTmp()) delete &listV;
		delete this; return;
	}
	// GetSortIdx the elimination list to descending order
	int dim = listV.GetItsDim();
	if(dim == 0) { if(listV.IsTmp()) delete &listV; return;}
	CNVEC_INT tmList(listV);
	tmList.Sort(DESCEND);
	// If list contains index above boundary then throw error
	if(tmList[0] >= GetItsDim()) {__ERRORTHROW(ERR_MEMBDR);}
	// Remove data contains same elements
	int nn;
	int o1, n1;
	int *rv = tmList.GetItsVec();
	o1 = rv[0];
	Remove(o1);
	for(nn=1; nn < dim ; nn++){
		n1 = rv[nn];
		if(o1 != n1) {
			Remove(n1);
			o1 = n1;
		}
	}
	if (listV.IsTmp()) delete &listV;
}

// Remove elements in vector which are appeared in listV
TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::RemElem(INST_CNVEC &listV)
{
	int dim = listV.GetItsDim();
	// GetSortIdx the elimination list to descending order
	INST_CNVEC* tmpV;
	SAFEALLOC(INST_CNVEC, tmpV);
	int nn, mm;
	bool ok = Tmp;
	if(ok) Tmp = false;
	CNVEC_INT elimlist;
	if( this != &listV ) {
		if( dim != 0 && itsDim != 0) {
			// Remove data which are the same elements in listV
			T *rv = listV.GetItsVec();
			for(nn = 0; nn < itsDim; ++nn) {
				T val = itsVec[nn];
				for(mm = 0; mm < dim; ++mm) {
					if(val == rv[mm]) { elimlist.Add(nn); tmpV->Add(val); break;}
				}
			}
			Remove(elimlist);
		}
		if(listV.IsTmp()) delete &listV;
	} else {
		(*tmpV) = (*this);
		SetDim(0);
	}
	if(ok) delete this;
	tmpV->SetTmp(true);
	return *tmpV;
}

TEMPLATE_CNVEC
INST_CNVEC & INST_CNVEC::RemElem(T val)
{
	// Remove data contains same elements
	int nn;
	INST_CNVEC* tmpV;
	SAFEALLOC(INST_CNVEC, tmpV);
	tmpV->SetDim(0);
	CNVEC_INT elimlist;
	for(nn = 0; nn < itsDim; ++nn) {
		if(val== itsVec[nn]) { 
			elimlist.Add(nn);
			tmpV->Add(val);
		}
	}
	Remove(elimlist);
	if(Tmp) delete this;
	tmpV->SetTmp(true);
	return *tmpV;
}

TEMPLATE_CNVEC
void INST_CNVEC::Sort(SORT_METHOD srt)
{
	if(itsDim < 2) {return ;}
	int i, j, ord, ord2;
	INST_CNVEC tmpV;
	tmpV.SetDim(itsDim);
	T *rv = tmpV.GetItsVec();
	ord=1;
	while(ord < itsDim) {
		ord2 = 2*ord;
		for(i=0;i<itsDim;i+=ord2){
			int j1=i;
			int j2 = i+ ord;
			if(j2>=itsDim) {
				while(j1<itsDim) {
					rv[j1] = itsVec[j1];
					j1++;
				}
				continue;
			}
			j=i;
			while(j<ord2+i && j!=itsDim) {
				if(srt) {
					if(itsVec[j1]<itsVec[j2]) {
						rv[j++]=itsVec[j1];
						j1++;
						if(j1 == i+ord) {
							while(j2<ord2+i && j != itsDim) rv[j++] = itsVec[j2++];
						}
					}else {
						rv[j++]=itsVec[j2];
						j2++;
						if (j2 == i+ord2 || j2 == itsDim) {
							while(j1<ord+i) rv[j++] = itsVec[j1++];
						}
					}
				}else {
					if(itsVec[j1]>itsVec[j2]) {
						rv[j++]=itsVec[j1];
						j1++;
						if(j1 == i+ord) {
							while(j2<ord2+i && j != itsDim) rv[j++] = itsVec[j2++];
						}
					}else {
						rv[j++]=itsVec[j2];
						j2++;
						if (j2 == i+ord2 || j2 == itsDim) {
							while(j1<ord+i) rv[j++] = itsVec[j1++];
						}
					}
				}
			}
		}
		ord = ord2;
		*this = tmpV;
	}
}


TEMPLATE_CNVEC
CNVEC_INT & INST_CNVEC::GetSortIdx(SORT_METHOD srt)
{
	CNVEC_INT *idxvec;
	SAFEALLOC(CNVEC_INT,idxvec);
	idxvec->SetDim(itsDim);
	idxvec->FillInc();
	if(itsDim < 2) {
		idxvec->SetTmp(true);
		return *idxvec;
	}
	int i, j, ord, ord2;
	CNVEC_INT tmpV;
	tmpV.SetDim(itsDim);
	int *rv = tmpV.GetItsVec();
	ord=1;
	while(ord < itsDim) {
		int *idx = idxvec->GetItsVec();
		ord2 = 2*ord;
		for(i=0;i<itsDim;i+=ord2){
			int j1=i;
			int j2 = i+ ord;
			if(j2>=itsDim) {
				while(j1<itsDim) {
					rv[j1] = idx[j1];
					j1++;
				}
				continue;
			}
			j=i;
			while(j<ord2+i && j!=itsDim) {
				if(srt) {
					if(itsVec[idx[j1]]<itsVec[idx[j2]]) {
						rv[j++]=idx[j1];
						j1++;
						if(j1 == i+ord) {
							while(j2<ord2+i && j != itsDim) rv[j++] = idx[j2++];
						}
					}else {
						rv[j++]=idx[j2];
						j2++;
						if (j2 == i+ord2 || j2 == itsDim) {
							while(j1<ord+i) rv[j++] = idx[j1++];
						}
					}
				}else {
					if(itsVec[idx[j1]]>itsVec[idx[j2]]) {
						rv[j++]=idx[j1];
						j1++;
						if(j1 == i+ord) {
							while(j2<ord2+i && j != itsDim) rv[j++] = idx[j2++];
						}
					}else {
						rv[j++]=idx[j2];
						j2++;
						if (j2 == i+ord2 || j2 == itsDim) {
							while(j1<ord+i) rv[j++] = idx[j1++];
						}
					}
				}
			}
		}
		ord = ord2;
		*idxvec = tmpV;
	}
	idxvec->SetTmp(true);
	return *idxvec;
}
