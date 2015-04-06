#include<fstream>
#include<iostream>
#include<strstream>
#include<time.h>
#include <sys/timeb.h>
#include "exception.hpp"

#ifndef __GENERAL_HPP
#define __GENERAL_HPP
using std::ofstream;
using std::ifstream;
using std::ostream;
using std::istream;
using std::endl;
using std::cout;
using std::streambuf;
using std::ostrstream;

#define VERSION "ver 1.0"

#ifdef __LINUX
	#define POPEN popen
	#define PCLOSE pclose
#else 
	#define POPEN _popen
	#define PCLOSE _pclose
#endif

//#define _CRT_SECURE_NO_WARNINGS 
#define NONE   0
#define LOW    1
#define MEDIUM 2
#define HIGH   3
#define DEBUGLVL NONE                // debugging mode

#define MAXINPUTITEM 100          // Available number of input parameters used in reading a control file in CCTRL class
#define MAXLINE     1000          // Maximum length of a line 
#define MAXCVRMS   999.0          // Maximum CVRMS value
#define INVCUTOFF    1.0          // discard cluster funtion with d value less than INVCUTOFF, d = nu'*nu - mu'*t;
#define MEM_EXPAND    20          // memory allocation unit for efficiency
#define MAXCFORDER     6          // maximum order of CFs

#ifndef __EXCEPTION_HPP
	#define __STARTEXCEPTION__
	#define __ENDEXCEPTION__
	#define __ERRORTHROW(x)
#endif

// define global variables :start
#define EXTERN  extern struct timeb  __time_start; extern CNVEC_INT __tmpVecUINT;\
				extern CNVEC_DBL __tmpVecDBL; extern CNMAT_INT __tmpMatUINT;\
				extern CNMAT_DBL __tmpMatDBL; 
#define __TIMESTART__  ftime(&__time_start);
#define __TIMEEND__  {	struct timeb __time_end;\
						ftime(&__time_end);\
						int etime = (int)(__time_end.time-__time_start.time)*1000 +\
						            (__time_end.millitm - __time_start.millitm);\
					    int ms = etime%1000; etime = etime/1000;\
						int ss = etime%60; etime = etime/60;\
						int mm = etime%60; int hh = etime/60;\
						char timebuff[300];\
						sprintf(timebuff,"%02d:%02d:%02d.%03d",hh,mm,ss,ms);\
						cout << "Normal termination [running time: "<<timebuff<<"]" << endl;}
#define __PRINTTIME__  {	time_t __time_end;\
						time(&__time_end);\
						cout << ctime(&__time_end);}
#define EXTERNDEF struct timeb __time_start; CNVEC_INT __tmpVecUINT; CNVEC_DBL __tmpVecDBL;\
					CNMAT_INT __tmpMatUINT; CNMAT_DBL __tmpMatDBL;
// define global variables :end

#define UINT unsigned int
#define USRT unsigned short

#define DIRECTION bool
#define ROW    false
#define COLUMN true

#define SORT_METHOD bool
#define DESCEND false
#define ASCEND  true

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define SWAP_T(type,a,b) { type temp=(a);(a)=(b);(b)=temp; }

#define WRITEDATA(fp,x) fp.write((char *) &x, sizeof(x))
#define READDATA(fp,x) fp.read((char *) &x, sizeof(x))
#define	WRITEARRDATA(fp,x,y) fp.write((char *) x, sizeof(*x)*(y)) 
#define	READARRDATA(fp,x,y) fp.read((char *) x, sizeof(*x)*(y)) 
#define APPENDTIME  {	struct timeb __time_end;\
						ftime(&__time_end);\
						int etime = (int)(__time_end.time-__time_start.time)*1000 +\
						            (__time_end.millitm - __time_start.millitm);\
					    int ms = etime%1000; etime = etime/1000;\
						int ss = etime%60; etime = etime/60;\
						int mm = etime%60; int hh = etime/60;\
						char timebuff[300];\
						sprintf(timebuff,"%02d:%02d:%02d.%03d",hh,mm,ss,ms);\
						cout << "["<<timebuff<<"]";}
#define LOG(x) { cout << #x <<endl; }
#define ERRMSG(x)\
    {cout<< "\nERROR: " << #x << "\n"<< endl;}
#define PRINTLOG(x) {	cout << #x; \
						struct timeb __time_End;\
						ftime(&__time_End);\
						cout << " [" << __time_End.time - __time_start.time<< "s]" << endl;}
#define EVALLOG(x)\
    {cout<< "\nEVAL: (" << #x << ")\n  on line : "<< __LINE__<< \
                   "\n  in file : "<< __FILE__ << "\n" << x << endl;}
//#define APPENDTIME  {	struct timeb __time_End;\
//						ftime(&__time_End);\
//						cout << "[" << __time_End.time - __time_start.time<< "s]";}

#define RUN_THIS(np,xcode)	if(sStep == cStep) { xcode; ++sStep; ++cStep; Save(np);} else{ ++cStep;}

#define __LOGOUT 

#ifdef ASSERT
	#undef ASSERT
#endif
#ifdef EVAL
	#undef EVAL
#endif
#ifdef PRINT
	#undef PRINT
#endif


#if DEBUGLVL<LOW
	#define ASSERT(x)  { }
	#define TAKECARE(x) { }
	#define BECARE(x)  { }
	#define PRINTLOW(x)  { }
	#define EVALLOW(x) { }
	#define PRINTLOWONE(x,y)  { }
	#define PRINTLOWTWO(x,y,z) { }
	#define __DEBUGTIMEEND__ { }
#else
	#define ASSERT(x)\
		if(!(x))\
		{\
			cout << "\n ASSERT ERROR:" << #x << " is failed\n"; \
			cout << " on line " << __LINE__ << endl; \
			cout << " in file " << __FILE__ << endl; \
		}
	#define TAKECARE(x)\
		{\
			cout << "\n Possible error :" << #x << "\n"; \
			cout << " on line " << __LINE__ << endl; \
			cout << " in file " << __FILE__ << endl; \
		}
	#define PRINTLOW(x) {	cout << #x; \
							struct timeb __time_End;\
							ftime(&__time_End);\
							cout << " : (" << __time_End.time - __time_start.time<< " s)" << endl;}
	#define EVALLOW(x)\
        {cout<< "\nEVAL: (" << #x << ")\n  on line : "<< __LINE__<< \
                       "\n  in file : "<< __FILE__ << "\n" << x << endl;}
	#define PRINTLOWONE(x,y) {	cout << x << "(" << y << ")" << endl;}
	#define PRINTLOWTWO(x,y,z) {cout << x << "(" << y << "," << z << ")" << endl;}
	#define __DEBUGTIMEEND__  {	struct timeb __time_end;\
						ftime(&__time_end);\
						int etime = (int)(__time_end.time-__time_start.time)*1000 +\
						           -(__time_end.millitm - __time_start.millitm);\
					    int ms = etime%1000; etime = etime/1000;\
						int ss = etime%60; etime = etime/60;\
						int mm = etime%60; int hh = etime/60;\
						char timebuff[300];\
						sprintf(timebuff,"%02d:%02d:%02d.%03d",hh,mm,ss,ms);\
						cout << "Normal termination [running time: "<<timebuff<<"]" << endl;}
#endif

#if DEBUGLVL<MEDIUM
	#define PRINTMED(x)  { }
	#define EVALMED(x) { }
	#define PRINTMEDONE(x,y) { }
	#define PRINTMEDTWO(x,y,z) { } 
	#define PRINTMIDREMVAR() { }
#else
	#define PRINTMED(x) {	cout << #x; \
							time_t __time_End;\
							time(&__time_End);\
							cout << " : (" << __time_End - __time_start<< " s)" << endl;}
	#define EVALMED(x)\
        {cout<< "\nEVAL: (" << #x << ")\n  on line : "<< __LINE__<< \
                       "\n  in file : "<< __FILE__ << "\n" << x << endl;}
	#define PRINTMEDONE(x,y) {	cout << x << "(" << y << ")" << endl;}
	#define PRINTMEDTWO(x,y,z) {cout << x << "(" << y << "," << z << ")" << endl;}
	#define PRINTMIDREMVAR() { cout << "(Str,Vec4,Vec8,Mat4,Mat8) : ( "; \
							cout << ", " << __tmpVecUINT.nVec;\
							cout << ", " << __tmpVecDBL.nVec;\
							cout << ", " << __tmpMatUINT.nMat;\
							cout << ", " << __tmpMatDBL.nMat << ")" << endl;}
#endif

#if DEBUGLVL<HIGH
	#define EVALHIGH(x) { }
	#define PRINTHIGH(x)  { }
	#define PRINTHIGHONE(x,y)  { }
	#define PRINTHIGHTWO(x,y,z)  { }
#else
	#define EVALHIGH(x)\
        {cout<< "\nEVAL: (" << #x << ")\n  on line : "<< __LINE__<< \
                       "\n  in file : "<< __FILE__ << "\n" << x << endl;}
	#define PRINTHIGH(x) { cout << #x <<endl; }
	#define PRINTHIGHONE(x,y) {	cout << x << "(" << y << ")" << endl;}
	#define PRINTHIGHTWO(x,y,z) {cout << x << "(" << y << "," << z << ")" << endl;}
#endif

// Some algorithm
#define SAFEALLOC(type,object) object = new type;  ASSERT(object !=NULL);if(object == NULL) {__ERRORTHROW(ERR_MEMALOC); }

// Template types
#define CNVEC_INT  CNVEC<  int,  0>
#define CNVEC_DBL  CNVEC<double, 0>

#define CNMAT_INT  CNMAT<   int, 0>
#define CNMAT_DBL  CNMAT<double, 0>

#endif // __GENERAL_HPP

