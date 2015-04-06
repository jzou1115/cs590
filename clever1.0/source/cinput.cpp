#include "cinput.hpp"

CINPUT::CINPUT(void)
{
	 idName = " "; ptrVar = NULL; nMode = 0; mFlag = false; mType=0;
}

CINPUT::~CINPUT(void)
{
}

void CINPUT::operator() (CSTRING tname, void *tptr, int tmode, int tvalue, bool tflag)
{
	idName = tname; 	ptrVar = tptr; 	mFlag = tflag;	nMode = tmode;	mType = 0;
	*((int *)tptr) = tvalue;
}

void CINPUT::operator() (CSTRING tname, void *tptr, int tmode, double tvalue, bool tflag)
{
	idName = tname; 	ptrVar = tptr; 	mFlag = tflag;	nMode = tmode;	mType = 1;
	*((double *)tptr) = tvalue;
}
void CINPUT::operator() (CSTRING tname, void *tptr, int tmode, CSTRING tvalue, bool tflag)
{
	idName = tname; 	ptrVar = tptr; 	mFlag = tflag;	nMode = tmode;	mType = 2;
	*((CSTRING *)tptr) = tvalue;
}

void CINPUT::SetValue(CSTRING tvalue)
{
	switch(mType) {
		case 0:
			*((int *)ptrVar) = tvalue.CvtInt(); 
			break;
		case 1:
			*((double *)ptrVar) = tvalue.CvtDbl(); 
			break;
		case 2:
			*((CSTRING *)ptrVar) = tvalue; 
			break;
	}
}

// If values are equal then return 0 or else return 1
int CINPUT::Compare(CINPUT & input)
{
	switch(mType) {
		case 0:
			if(*((int *)(ptrVar)) == *((int *)(input.ptrVar))) return 0;
			else return 1;
		case 1:
			if(*((double *)(ptrVar)) == *((double *)(input.ptrVar))) return 0;
			else return 1;
		case 2:
			if(*((CSTRING *)(ptrVar)) == *((CSTRING *)(input.ptrVar))) return 0;
			else return 1;
	}
	return 0;
}

CINPUT & CINPUT::operator= (CINPUT &rhs)
{
	idName = rhs.idName; 	mFlag = rhs.mFlag;	nMode = rhs.nMode;	mType = rhs.mType;
	if(rhs.ptrVar == NULL) return *this;
	switch(mType) {
		case 0:
			*((int *)ptrVar) = *((int *)rhs.ptrVar); 
			break;
		case 1:
			*((double *)ptrVar) = *((double *)rhs.ptrVar); 
			break;
		case 2:
			*((CSTRING *)ptrVar) = *((CSTRING *)rhs.ptrVar); 
			break;
	}
	return *this;
}

