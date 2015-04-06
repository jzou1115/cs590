#pragma once
#include "general.hpp"
#include <string>
#include <vector>

class CSTRING : public std::string
{
public:
	CSTRING();
	CSTRING(char);
	CSTRING(char *);
	CSTRING(std::string);
	~CSTRING();
	char & operator[] (int offset);
	CSTRING & operator= (CSTRING &);
	CSTRING & operator= (char);
	CSTRING & operator= (char *);
	CSTRING & operator= (const char *);
	CSTRING & operator= (std::string);
	CSTRING   operator+ (CSTRING &);    // concatenation of two string
	CSTRING & operator+=(char);
	CSTRING & operator+=(std::string);
	CSTRING & operator+=(CSTRING &);

	CSTRING  operator() (int offset);  // return substring value;
	bool ReadLine(FILE *fp);            // Read a line in the file pointer and save in the string
	bool RemoveComment();               // Remove comments in the line: <ex>"// #"
	int SplitLine(const char* seps);   // split a string with the criterion "seps" and save in the subString
	                                    // return value : number of elements in the subString.
	bool IsStrIn(CSTRING str);          // find "str" in string; return value: if "str" in the string (true) else (false)
	CSTRING  Lower();                   // convert capital characters to lower characters; return converted string
	int GetLen();                      // return length of string 
	int CvtInt();                       // convert the string into integer
	double CvtDbl();                    // convert the string into double
	bool IsNumeric();                   // Check whether the string is numeric
	const char * GetString() const;     // return string pointer with "c" string type
	void SetString(const char *rhs);    // set string with "rhs"

public:
	std::vector<std::string> subString;  // it contains divided strings
};

ofstream & operator<< (ofstream& fout, CSTRING &rhs);
ifstream & operator>> (ifstream& fin, CSTRING &rhs);
ostream & operator<< (ostream& cout, CSTRING &rhs);
