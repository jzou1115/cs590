#include <iostream>
#include <fstream>
#include "exception.hpp"
using std::ofstream;
using std::endl;
using std::cout;
//extern ostream cout;

/*
 Definition of exception code
	MEMORYBOUNDARY  = offset exceeds memory boundary
*/
#define PRINTERR(x) {cout << "\n" << #x << "\nat line: " << itsLine << "\n in file:" << itsFile << endl;}

CEXCEPTION::CEXCEPTION(ERRORCODE code, char *file , const int line)
{
	xCode = code;
	itsFile = file;
	itsLine = line;
}

void CEXCEPTION::Print()
{
	switch(xCode) {
	case ERR_MEMBDR: 
		PRINTERR(Error: offset exceeds memory boundary);
		break;
	case ERR_XSORT:
		PRINTERR(Error: sorted vector should have elements);
		break;
	case ERR_VECDIM:
		PRINTERR(Error: dimension of two vectors should be same);
		break;
	case ERR_MATDIM:
		PRINTERR(Error: dimension of two matrixes should be same);
		break;
	case ERR_MATVEC:
		PRINTERR(Error: dimension mismatch between vector and matrix);
		break;
	case ERR_MATMAT:
		PRINTERR(Error: dimension mismatch between matrices);
		break;
	case ERR_MEMALOC:
		PRINTERR(Error: memory allocation error);
		break;
	case ERR_NODIM:
		PRINTERR(Error: check dimension of vector or matrix);
		break;
	case ERR_DIVZERO:
		PRINTERR(Error: divided by zero. Increase mnNSeq value in control file);
		break;
	case ERR_CONVFAIL:
		PRINTERR(Error: fail in converting type);
		break;
	case ERR_SEQFORM:
		PRINTERR(Error: wrong sequence file format);
		break;
	case ERR_DSNFORM:
		PRINTERR(Error: wrong design file format);
		break;
	case ERR_DATAFORM:
		PRINTERR(Error: wrong data format);
		break;
	case ERR_CCFS_BND:
		PRINTERR(Error: error in CCFS class);
		break;
	case ERR_POPEN:
		PRINTERR(Error: Process is not created);
		break;
	case ERR_PCLOSE:
		PRINTERR(Error: Cannot close pipe by using pclose function);
		break;
	case ERR_NOELEM:
		PRINTERR(Error: No vector or matrix elements);
		break;
	case ERR_GENSEQ:
		PRINTERR(Error: Error in Generating sequences);
		break;
	case ERR_ELIMSEQS:
		PRINTERR(Error: Error in ElimSeqs routine in CCFS class);
		break;
	case ERR_NOCFDEFS:
		PRINTERR(Error: There is no defintion of a cluster function);
		break;
	}
}
