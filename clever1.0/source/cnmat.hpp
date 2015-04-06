#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "general.hpp"
#include "cstring.hpp"

__LOGOUT
#ifndef MATH_TEMPLATE
  #define MATH_TEMPLATE
  #define TEMPLATE_CNVEC template<class T, int _i>
  #define INST_CNVEC CNVEC<T,_i>
  #define TEMPLATE_CNMAT template<class T, int _i>
  #define INST_CNMAT CNMAT<T,_i>
#endif

TEMPLATE_CNVEC class CNVEC;

// Define Matrix class
// CNMAT definition
TEMPLATE_CNMAT
class CNMAT 
{
private:
	CNVEC<INST_CNVEC*, _i> itsMat;  // Save column pointer
	int itsNRow;                   // dimension of Row
	int itsNCol;                   // dimension of Column
	bool    Tmp;                    // denote temperary matrix or not
public:
	static int nMat;
	CNMAT();
	CNMAT(int, int);
	CNMAT(CNMAT &);
	~CNMAT();
	void Clear();
	void SetTmp(bool tmp) { Tmp = tmp; }
	bool IsTmp() { return Tmp; }
	void SetDim(int i, int j);          // set matrix dimension
	int Len(DIRECTION);
	int const GetItsDim(DIRECTION drc) const // return dimension of matrix
		{  if (drc) return itsNCol; else return itsNRow;}
	void Fill(T n);                           // Set values with "n"
	CNMAT &Tp();                              // return a transposed matrix
	void SetZero();                           // set elements to zero 
	void Add(DIRECTION, CNMAT &);             // append a matrix to this matrix
	void Add(DIRECTION, INST_CNVEC &);        // append a vector to this matrix
	void Remove(DIRECTION, int offset);      // remove a vector from this matrix
	void Remove(DIRECTION, int offset1, int offset2);  // remove vectors from "offset1" to "offset2"
	void Remove(DIRECTION drc, CNVEC_INT &listV);       // remove vectors with offsets in "listV"
	INST_CNVEC & Sum(DIRECTION drc);		  // return summation results of each vector
	void Read(ifstream &);                    // read a matrix from a file
	void Write(ofstream &);                   // write a matrix to a file
	void Print();                             // print this matrix to standard output
	INST_CNVEC & GetVec(DIRECTION drc, const int offset);      // Return vector with the direction and offset
	CNMAT & GetSelect(DIRECTION drc, CNVEC_INT & vec);         // return part of matrix
	CNMAT & GetSelect(DIRECTION drc, int st, int end);        // return part of matrix
	CNMAT_DBL & SMTransMat(CNVEC_DBL &diag);                    // return matrix, MT*diag*M


	// Define operators
	INST_CNVEC & operator[] (int offset);
	CNMAT & operator+ (CNMAT &);
	CNMAT & operator+ (const T &);    // add a given number to a matrix
	CNMAT & operator+ ();
	CNMAT & operator+= (const T);     // add a given number to a matrix
	CNMAT & operator+= (CNMAT &rhs);  // add a given number to a matrix
	CNMAT & operator-= (const T);  
	CNMAT & operator-= (CNMAT &rhs);  
	CNMAT & operator- (CNMAT &);
	CNMAT & operator- (const T &);
	CNMAT & operator- ();
	CNMAT & operator* (CNMAT & rhs);
	INST_CNVEC & operator* (INST_CNVEC &);      // It returns (M*v) in matlab, v is regarded as a column vector
	CNMAT & operator* (const T &);              // Matrix*constant
	CNVEC_DBL & DotDbl (CNVEC_DBL &);           // Return double type vector
	CNMAT & operator/ (CNMAT &);                // It returns(M./M) in matlab(?)
	CNMAT & operator/ (const T);                // Matrix/constant
	CNMAT & operator= (CNMAT &);
	T& operator() (int, int);                 // return an element of matrix
	CNMAT& operator() (int, int, int, int); // vector dimension setting

};

TEMPLATE_CNMAT 
ofstream & operator<< (ofstream& fout, INST_CNMAT &rhs);
TEMPLATE_CNMAT 
ifstream & operator>> (ifstream& fin,  INST_CNMAT &rhs);
TEMPLATE_CNMAT 
ostream & operator<< (ostream& cout, INST_CNMAT &rhs);
TEMPLATE_CNMAT 
istream & operator>> (istream& cin,  INST_CNMAT &rhs);

#include "cnmat.cpp"

// global operators
CNVEC_DBL & operator* (CNMAT_DBL &mat, CNVEC_INT &vec);
CNVEC_DBL & operator* (CNVEC_INT &vec, CNMAT_DBL &mat);
CNVEC_DBL & operator* (CNMAT_INT &mat, CNVEC_DBL &vec);
CNVEC_DBL & operator* (CNVEC_DBL &vec, CNMAT_INT &mat);

CNMAT_DBL & operator* (CNMAT_DBL & mat1, CNMAT_INT &mat2);
CNMAT_DBL & operator* (CNMAT_INT & mat1, CNMAT_DBL &mat2);

