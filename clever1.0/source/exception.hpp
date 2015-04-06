#ifndef __EXCEPTION_HPP
#define __EXCEPTION_HPP

#define __STARTEXCEPTION__  try {
#define __ENDEXCEPTION__  } catch(CEXCEPTION & theException) { theException.Print(); }
#define __ERRORTHROW(x) throw CEXCEPTION(x, __FILE__, __LINE__)

enum ERRORCODE { 
	ERR_MEMBDR,         // Offset exceeds memory boundary 
	ERR_XSORT,          // Sorted vector should have elements. 
	ERR_VECDIM,         // Dimension of two vectors should be same
	ERR_MATDIM,         // Dimension of two matrix should be same
	ERR_MATVEC,         // Dimension mismatch between vector and matrix
	ERR_MATMAT,         // Dimension mismatch between matrix and matrix
	ERR_MEMALOC,        // Memory allocation error
	ERR_NODIM,          // Check dimension of vector or matrix
	ERR_DIVZERO,        // value divided by zero
	ERR_CONVFAIL,       // type conversion failed
	ERR_SEQFORM,        // wrong sequence file format
	ERR_DSNFORM,        // wrong design file format
	ERR_DATAFORM,       // error occured in reading vector and matrix 
	ERR_CCFS_BND,       // error in CCFS class
	ERR_POPEN,          // POpen error
	ERR_PCLOSE,         // pclose error
	ERR_NOELEM,         // No vector and matrix elements
	ERR_GENSEQ,         // Error in generating sequences
	ERR_ELIMSEQS,       // Error in ElimSeqs routine in CCFS class
	ERR_NOCFDEFS        // There is no definition of cluster function
};

class CEXCEPTION
{
public:
	CEXCEPTION(ERRORCODE code) {xCode = code; }
	CEXCEPTION(ERRORCODE code, char *, const int);
	void Print();
	void SetCode(ERRORCODE code) { xCode = code; }

private:
	ERRORCODE xCode;
	char *itsFile;
	int itsLine;
};

#endif // __EXCEPTION_HPP

