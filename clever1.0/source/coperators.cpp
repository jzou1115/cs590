// Global operators
#include "cnvec.hpp"
#include "cnmat.hpp"
__LOGOUT
// Golbal matrix operation
CNMAT_DBL & operator* (CNMAT_DBL &mat1, CNMAT_INT &mat2)
{
	if(mat1.GetItsDim(COLUMN) != mat2.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATMAT); }
	int dim = mat1.GetItsDim(COLUMN);
	if( dim == 0) {
		if(mat2.IsTmp()) delete &mat2;
		return mat1;
	}
	CNMAT_DBL* vMat;
	SAFEALLOC(CNMAT_DBL,vMat);
	int nrow = mat1.GetItsDim(ROW);
	int ncol = mat2.GetItsDim(COLUMN);
	vMat->SetDim(nrow,ncol);

	int ni, nj, nk;
	double *av, *cv, bval;
	int *bv;
	bool ok1 = mat1.IsTmp();
	bool ok2 = mat2.IsTmp();
	if(ok1)mat1.SetTmp(false);
	if(ok2)mat2.SetTmp(false);
	for(nk=0;nk<ncol; nk++) {
		cv = (vMat->GetVec(COLUMN,nk)).GetItsVec();
		bv = mat2.GetVec(COLUMN,nk).GetItsVec();
		for(ni=0;ni<dim; ni++) {
			av = mat1.GetVec(COLUMN,ni).GetItsVec();
			bval = (double)bv[ni];
			for(nj=0;nj<nrow;nj++) {
				cv[nj] += av[nj]*bval;
			}
		}
	}
	if(ok1) delete &mat1;
	if(ok2) delete &mat2;
	vMat->SetTmp(true);
	return *vMat;
}

CNMAT_DBL & operator* (CNMAT_INT &mat1, CNMAT_DBL &mat2)
{
	if(mat1.GetItsDim(COLUMN) != mat2.GetItsDim(ROW)) {__ERRORTHROW(ERR_MATMAT); }
	if(mat2.GetItsDim(COLUMN) == 0) {
		if(mat1.IsTmp()) delete &mat1;
		return mat2;
	}
	CNMAT_DBL* vMat;
	SAFEALLOC(CNMAT_DBL,vMat);
	int dim = mat1.GetItsDim(COLUMN);
	int nrow = mat1.GetItsDim(ROW);
	int ncol = mat2.GetItsDim(COLUMN);
	vMat->SetDim(nrow,ncol);

	int ni, nj, nk;
	double *cv,*bv, bval;;
	int *av;
	bool ok1 = mat1.IsTmp();
	bool ok2 = mat2.IsTmp();
	if(ok1)mat1.SetTmp(false);
	if(ok2)mat2.SetTmp(false);
	for(nk=0;nk<ncol; nk++) {
		cv = (vMat->GetVec(COLUMN,nk)).GetItsVec();
		bv = mat2.GetVec(COLUMN,nk).GetItsVec();
		for(ni=0;ni<dim; ni++) {
			av = mat1.GetVec(COLUMN,ni).GetItsVec();
			bval = bv[ni];
			for(nj=0;nj<nrow;nj++) {
				cv[nj] += (double)av[nj]*bval;
			}
		}
	}
	if(ok1) delete &mat1;
	if(ok2) delete &mat2;
	vMat->SetTmp(true);
	return *vMat;
}

// Matrix*vector definition
CNVEC_DBL & operator* (CNMAT_DBL &mat, CNVEC_INT &rhs)
{
	int col = mat.GetItsDim(COLUMN);
	if(col != rhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	if(col == 0) {
		if(mat.IsTmp()) delete &mat;
		return *vVec;
	}
	int dim = mat.GetItsDim(ROW);
	vVec->SetDim(dim);
	int ni, nj;
	double *av, *cv, bval;
	int *bv;
	cv = vVec->GetItsVec();
	bv = rhs.GetItsVec();
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	for(ni=0;ni<col; ni++) {
		av = mat[ni].GetItsVec();
		bval = bv[ni];
		for(nj=0;nj<dim;nj++) {
			cv[nj] += av[nj]*bval;
		}
	}
	if(rhs.IsTmp()) delete &rhs;
	if(ok) delete &mat;
	vVec->SetTmp(true);
	return *vVec;
}

CNVEC_DBL & operator* (CNMAT_INT &mat, CNVEC_DBL &rhs)
{
	int col = mat.GetItsDim(COLUMN);
	if(col != rhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	if(col == 0) {
		if(mat.IsTmp()) delete &mat;
		return rhs;
	}
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	int dim = mat.GetItsDim(ROW);
	vVec->SetDim(dim);
	int ni, nj;
	int *av;
	double *cv, bval;
	double *bv;
	cv = vVec->GetItsVec();
	bv = rhs.GetItsVec();
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	for(ni=0;ni<col; ni++) {
		av = mat[ni].GetItsVec();
		bval = bv[ni];
		for(nj=0;nj<dim;nj++) {
			cv[nj] += av[nj]*bval;
		}
	}
	if(rhs.IsTmp()) delete &rhs;
	if(ok) delete &mat;
	vVec->SetTmp(true);
	return *vVec;
}

// Global vector operation
CNVEC_DBL & operator* (CNVEC_INT &lhs, CNMAT_DBL &mat)
{
	int row = mat.GetItsDim(ROW);
	if(row != lhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	if(row == 0) {
		if(mat.IsTmp()) delete &mat;
		return *vVec;
	}
	int dim = mat.GetItsDim(COLUMN);
	vVec->SetDim(dim);
	int ni, nj;
	int *av;
	double *cv, *bv, bval;
	cv = vVec->GetItsVec();
	av = lhs.GetItsVec();
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	for(ni=0;ni<dim; ni++) {
		bv = mat[ni].GetItsVec();
		bval = 0.0;
		for(nj=0;nj<row;nj++) {
			bval += av[nj]*bv[nj];
		}
		cv[ni] = bval;
	}
	if(lhs.IsTmp()) delete &lhs;
	if(ok) delete &mat;
	vVec->SetTmp(true);
	return *vVec;
}


CNVEC_DBL & operator* (CNVEC_DBL &lhs, CNMAT_INT &mat)
{
	int row = mat.GetItsDim(ROW);
	if(row != lhs.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	if(row == 0) {
		if(mat.IsTmp()) delete &mat;
		return lhs;
	}
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	int dim = mat.GetItsDim(COLUMN);
	vVec->SetDim(dim);
	int ni, nj;
	double *av, *cv, bval;
	int *bv;
	cv = vVec->GetItsVec();
	av = lhs.GetItsVec();
	bool ok = mat.IsTmp();
	if(ok) mat.SetTmp(false);
	for(ni=0;ni<dim; ni++) {
		bv = mat[ni].GetItsVec();
		bval = 0.0;
		for(nj=0;nj<row;nj++) {
			bval += av[nj]*bv[nj];
		}
		cv[ni] = bval;
	}
	if(lhs.IsTmp()) delete &lhs;
	if(ok) delete &mat;
	vVec->SetTmp(true);
	return *vVec;
}




double operator* (CNVEC_DBL & vec1, CNVEC_INT &vec2)
{
	int dim = vec1.GetItsDim();
	if(dim != vec2.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	double *v1 = vec1.GetItsVec();
	int    *v2 = vec2.GetItsVec();
	double value=0.0;
	for(int i=0;i<dim; i++) {
		value += v1[i]*v2[i];
	}
	if(vec1.IsTmp()) delete &vec1;
	if(vec2.IsTmp()) delete &vec2;
	return value;
}

double operator* (CNVEC_INT &vec1, CNVEC_DBL & vec2)
{
	int dim = vec1.GetItsDim();
	if(dim != vec2.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	int    *v1 = vec1.GetItsVec();
	double *v2 = vec2.GetItsVec();
	double value=0.0;
	for(int i=0;i<dim; i++) {
		value += v1[i]*v2[i];
	}
	if(vec1.IsTmp()) delete &vec1;
	if(vec2.IsTmp()) delete &vec2;
	return value;
}

CNVEC_DBL & operator- (CNVEC_DBL & vec1, CNVEC_INT &vec2)
{
	int dim = vec1.GetItsDim();
	if(dim != vec2.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	if(dim == 0 ) {
		if(vec2.IsTmp()) delete &vec2;
		return vec1;
	}
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	*vVec = vec1;
	double *cv = vVec->GetItsVec();
	int    *v2 = vec2.GetItsVec();
	for(int i=0;i<dim; i++) {
		cv[i] -= v2[i];
	}
	if(vec1.IsTmp()) delete &vec1;
	if(vec2.IsTmp()) delete &vec2;
	vVec->SetTmp(true);
	return *vVec;
}

CNVEC_DBL & operator- (CNVEC_INT & vec1, CNVEC_DBL &vec2)
{
	int dim = vec1.GetItsDim();
	if(dim != vec2.GetItsDim()) {__ERRORTHROW(ERR_MATVEC); }
	if(dim == 0 ) {
		if(vec1.IsTmp()) delete &vec1;
		return vec2;
	}
	CNVEC_DBL* vVec;
	SAFEALLOC(CNVEC_DBL,vVec);
	vVec->SetDim(dim);
	double *cv = vVec->GetItsVec();
	int    *v1 = vec1.GetItsVec();
	double *v2 = vec2.GetItsVec();
	for(int i=0;i<dim; i++) {
		cv[i] = v1[i] - v2[i];
	}
	if(vec1.IsTmp()) delete &vec1;
	if(vec2.IsTmp()) delete &vec2;
	vVec->SetTmp(true);
	return *vVec;
}

