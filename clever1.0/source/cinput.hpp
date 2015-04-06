#pragma once
#include "general.hpp"
#include "cstring.hpp"

class CINPUT { 
public:
	CSTRING idName; 
	void *ptrVar; 
	int  nMode; 
	int  mType;
	bool mFlag;
	~CINPUT(void);
	CINPUT();
	void SetValue(CSTRING tvalue);
	int Compare(CINPUT &);
	void operator() (CSTRING idName, void *ptrVar, int tmode, int     tvalue, bool tflag);
	void operator() (CSTRING idName, void *ptrVar, int tmode, double  tvalue, bool tflag);
	void operator() (CSTRING idName, void *ptrVar, int tmode, CSTRING tvalue, bool tflag);
	CINPUT & operator= (CINPUT &rhs);
};
