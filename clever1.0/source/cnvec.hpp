#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "general.hpp"
#include "cstring.hpp"

// Define Vector class
// CNVEC definition
__LOGOUT
#ifndef MATH_TEMPLATE
  #define MATH_TEMPLATE
  #define TEMPLATE_CNVEC template<class T, int _i>
  #define INST_CNVEC CNVEC<T,_i>
  #define TEMPLATE_CNMAT template<class T, int _i>
  #define INST_CNMAT CNMAT<T,_i>
#endif

TEMPLATE_CNMAT class CNMAT;

TEMPLATE_CNVEC
class CNVEC
{

private:
	  T*   itsVec;                // first point position (starting point of vector)
	int   itsDim;                // vector dimension
	int   maxDim;                // variable maximum vector dimension (real memory allocation)
	bool      Tmp;                // true: temporary vector after sending values will be destroid
		                          // This parameter is introduced to reduce memory loss
public:
	static int nVec;
	CNVEC();
	CNVEC(int);
	CNVEC(CNVEC &);
	~CNVEC();
	void SetTmp(bool tmp) { Tmp = tmp; }
	bool IsTmp() { return Tmp; }
	int const GetItsDim() const { return itsDim; };  // return dimension of vector
	int const GetItsMaxDim() const { return maxDim; };  // return dimension of vector
	int Len() {return itsDim;}             // return dimension of vector
	  T* GetItsVec() const { return itsVec;}          // return a pointer of the array of the vector
	void SetItsVec (T* vec) { itsVec = vec; }         // set the pointer to the vector
	void SetDim(int i);                              // change dimension of vector
	void SetMaxDim(int i);                           // adjust real memory
	void Reduce();                                    // remove additional memory space
	CNVEC_DBL &CvtDbl();            // convert vector elements to double
	CNVEC_INT &CvtUINT();          // convert vector elements to int
	T Sum();                        // return summation value of vector elements
	CNVEC & DirectDivide(CNVEC &);          // return v./m in matlab
	CNVEC & GetSelect(CNVEC_INT & vec);    // return a vector which contains the selected elements of the vector
	CNVEC & GetSelect(int st, int end);   // elements from "st" to "end"
	int CompVec(CNVEC &vec);                // return it < vec : -1, it == vec : 0, it > vec : 1

	// set vector dimension
	double Average();                            // return average value;
	CNVEC & DirectDot(CNVEC &);             // returns (v.*m) in matlab
	bool Copy(CNVEC &);               // copy a vector to this vector
	CNVEC & Copy(int offset1, int offset2); // copy part of vector from offset1 to offset2 and return the vector
	void Fill(T val);                      // Set all vector elements to val
	void Fill( T val, int n, int m);     // Set elements between n and m to val
	void Fill(T val, CNVEC_INT &);        // Set elements between n and m to val
	void FillInc();                        // Set elements from 0 to n-1 with increment by 1
	void FillInc(T val);                   // Set elements from val to n+val-1 with increments by 1
	void Read(ifstream &);                 // Read a vector object from file stream
	void Write(ofstream &);                // Write a vector object to file
	void Print();                          // Print vector elements to stdard output
	void Add(CNVEC &);                     // Append a vector to the end of this vector
	void Add(T val);                       // Append an element to the end of this vector
	void Insert(int offset, T val);        // Append an element to the end of this vector
	void Remove(int offset);              // Remove an element located at "offset"
	void Remove(int offset1, int offset2);  // Remove the elements located between "offset1" and "offset2"
	void Remove(CNVEC_INT &listV);        // Remove elements according to the offsets in the list vector
	CNVEC & RemElem(CNVEC &listV);         // remove elements which have the same value in listV
	CNVEC & RemElem(T val);                // remove elements which have the value "val"
	void SetZero();                        // initialize all elements to zero
	int CountElemExcept(T val);           // Count number of positions which contain element except "val"
	T GetEndVal() {return itsVec[itsDim-1]; }; // return the value of the last element
	CNVEC_INT & GetSortIdx(SORT_METHOD srt);    // Sort elements in a vector and return a vector containing ordered "offset" values
	                                       // this operation does not change the original vector.
	void Sort(SORT_METHOD srt);			   // Sort elements in this vector
	CNVEC_DBL & UnitVec();                 // return a unit vector of this vector
	double Magnitude();                    // return magnitude(norm) of this vector
	T Max();                               // return the maximum value
	double Skewness();                     // return skewness of distribution of vector elements 
	INST_CNMAT & Dyadic (CNVEC &rhs);      // return dyadic matrix, Mij = vi*vj
	int  FindElem(T val);                  // find an element "val" in the vector and return offset value, not finding return -1;

	// Define operators
	CNVEC & operator+ (CNVEC &);
	CNVEC & operator+ (const T &);
	CNVEC & operator+ (){return *this;}
	CNVEC & operator- (CNVEC &);
	CNVEC & operator- (const T);
	CNVEC & operator- ();
	T operator* (CNVEC &);               // It returns (v*m) in matlab
	CNVEC & operator* (const T &);       // vector*constant
	CNVEC & operator* (INST_CNMAT &rhs); // vector*Matrix
	CNVEC & operator/ (CNVEC &);         // It returns(v./m) in matlab
	CNVEC & operator/ (const T &);       // vector/constant
	CNVEC & operator= (CNVEC &rhs);
	CNVEC & operator= (T );
	bool operator== (CNVEC &);
	CNVEC & operator+= (CNVEC &rhs);
	CNVEC & operator+= (const T );
	CNVEC & operator-= (CNVEC &rhs);
	CNVEC & operator-= (const T );
	CNVEC_INT & operator|| (CNVEC_INT &rhs);
	CNVEC_INT & operator&& (CNVEC_INT &rhs);

	T &operator[] (int offset);
};

/************************ Related global functions ********************************/
TEMPLATE_CNVEC
ofstream & operator<< (ofstream& fout, INST_CNVEC &rhs);
TEMPLATE_CNVEC
ifstream & operator>> (ifstream& fin,  INST_CNVEC &rhs);
TEMPLATE_CNVEC
ostream & operator<< (ostream& cout, INST_CNVEC &rhs);
TEMPLATE_CNVEC 
istream & operator>> (istream& cin,  INST_CNVEC &rhs);
TEMPLATE_CNVEC 
INST_CNVEC & operator* (T lhs,  INST_CNVEC &rhs);
TEMPLATE_CNVEC 
INST_CNVEC & operator+ (T lhs,  INST_CNVEC &rhs);
TEMPLATE_CNVEC 
INST_CNVEC & operator- (T lhs,  INST_CNVEC &rhs);

#include "cnvec.cpp"
// global operators

double operator* (CNVEC_INT &vec2, CNVEC_DBL &vec1);
double operator* (CNVEC_DBL &vec1, CNVEC_INT &vec2);

CNVEC_DBL & operator- (CNVEC_DBL & vec1, CNVEC_INT &vec2);
CNVEC_DBL & operator- (CNVEC_INT & vec1, CNVEC_DBL &vec2);

