#include "cstring.hpp"
__LOGOUT

ofstream & operator<< (ofstream& fout, CSTRING &rhs)
{
	int len;
	len = rhs.GetLen();
	WRITEDATA(fout, len);
	if(rhs.GetLen() != 0)
		WRITEARRDATA(fout, rhs.c_str(), sizeof(char)*(rhs.GetLen()+1));
	return fout;
}

ifstream & operator>> (ifstream& fin, CSTRING &rhs)
{
	char tmp[MAXLINE];
	int len;

	READDATA(fin, len);
	if(len > MAXLINE) {__ERRORTHROW(ERR_DATAFORM);}
	if(len != 0) 
		READARRDATA(fin, tmp,sizeof(char)*(len+1));
	rhs = tmp;
	return fin;
}

ostream & operator<< (ostream& cout, CSTRING &rhs)
{
	cout << rhs.GetString();
	return cout;
}

CSTRING::~CSTRING() {}

CSTRING::CSTRING() : std::string(){}

CSTRING::CSTRING(char ch) : std::string()
{
	*this = ch;
}

CSTRING::CSTRING(char * str) : std::string(str){}

CSTRING::CSTRING(std::string rhs) : std::string(rhs){}

int CSTRING::GetLen() 
{
	return (int)(length());
}

void CSTRING::SetString(const char * rhs)
{
	assign(rhs);
}

const char * CSTRING::GetString() const
{
	return c_str(); 
}

/* Reads one line into a pre-allocated string str (ignores the trailing \n). Returns 1 upon success, 0 upon EOF or failure */
// Return string without new line and carriage return characters
bool CSTRING::ReadLine(FILE* fp) {
	char *tok;
	static char buff[MAXLINE];
	char* res = fgets(buff, MAXLINE, fp);
	SetString("");
	subString.clear();
	if(res == NULL) return false;
	tok = strtok(buff,"\r\n");
	if(tok != NULL) SetString(tok);
	return true;
}

// Eliminate comment in the string
bool CSTRING::RemoveComment() {
  bool value;
  std::string::size_type s1;
  value = false;
  // Eliminate comment lines
  s1 = find("#");
  if(s1 != std::string::npos) {
	  resize(s1);
	  value = true;
  }
  const char *cstr = "//";
  s1 = find(cstr);
  if(s1 != std::string::npos) {
	  resize(s1);
	  value = true;
  }
  return value;
}

/* Splits the given line into substrings based on the given tokens. Returns the number of tokens.
Note, this changes the passed string. */
int CSTRING::SplitLine(const char* seps) {
  std::string::size_type s1, s2;
  subString.clear();
  // Eliminate comment lines
  s1 = 0;
  while (1){
	  s1 = find_first_not_of(seps,s1);
	  if(s1 == std::string::npos) break;  // no more element found
	  s2 = find_first_of(seps,s1);
	  if(s2 == std::string::npos) {
		  subString.push_back(substr(s1,GetLen()-s1)); // get the last element 
		  break;
	  }
	  subString.push_back(substr(s1,s2-s1)); // get elements
	  s1 = s2;
  }
  return (int) subString.size();
}


/* Convert a string to an integer */
int CSTRING::CvtInt() {
  int ans, par;
  ans = 0;
  par = 1;
  int i;
  for(i=0;i<(int)length();i++) {
	  char c = (*this)[i];
	  if( c >= '0' && c <='9') {
		  ans = ans*10 + (int)(c - '0');
	  }else if (c == '-') {
		  par = -1;
	  } else if ( c == ' ' || c == '+') {
		  continue;
	  } else {
		  break;
	  }
  }
  return ans;
}

bool CSTRING::IsStrIn(CSTRING str)
{
	std::string::size_type s1;
	s1 = find(str.GetString());
	if(s1 != std::string::npos) {
	  return true;
	}
	return false;
}

/* Converts a string to an double */
double CSTRING::CvtDbl() {
  double ans, dec, par;
  bool pt = false;
  ans = 0.0;
  dec = 1.0;
  par = 1.0;
  int i;
  for(i=0;i<GetLen();i++) {
	  char c = (*this)[i];
	  if( c >= '0' && c <='9') {
		  ans = ans*10.0 + (int)(c - '0');
		  if(pt) dec *= 10.0;
	  } else if (c == '.') {
		  pt = true;
	  }else if (c == '-') {
		  par = -1.0;
	  } else if ( c == ' ' || c == '+') {
		  continue;
	  } else {
		  break;
	  }
  }
  ans /= (dec*par);
  return ans;
}

/* Check whether the string is numeric */
bool CSTRING::IsNumeric() {
  bool pt = false;  // cursor does not pass point position
  bool nst = false; // cursor does not bass number starting position
  bool nospace = true; // check spaces between number;
  int i;
  for(i=0;i<GetLen();i++) {
	  char c = (*this)[i];
	  if (c == '-' || c == '+') {
		  if(nst == true) return false;
		  nst = true;
	  }else if ( c >= '0' && c <='9') {
		  if(nospace == false) return false;
		  nst = true;
	  }else if (c == '.') {
		  if (pt == true || nospace == false) return false;
		  pt = true;
	  }else if ( c == ' ') {
		  if(nst == true) nospace = false;
		  continue;
	  } else {
		  return false;
	  }
  }
  return true;
}

CSTRING CSTRING::Lower()
{
	int i;
	CSTRING str;
	str = *this;
	for(i=0; i<GetLen(); i++){
		if(str[i] >= 'A' && str[i] <= 'Z') 
			str[i] = str[i] - 'A' + 'a';
	}
	return str;
}

CSTRING & CSTRING::operator=(char rhs)
{
	std::string *p1;
    subString.clear();
	p1 = this;
	*p1 = rhs;
	return *this;
}

CSTRING & CSTRING::operator=(char * rhs)
{
	assign(rhs);
    subString.clear();
	return *this;
}

CSTRING & CSTRING::operator=(const char * rhs)
{
	assign(rhs);
    subString.clear();
	return *this;
}

CSTRING & CSTRING::operator=(CSTRING & rhs)
{
	assign(rhs.GetString());
	subString = rhs.subString;
	return *this;
}

CSTRING & CSTRING::operator=(std::string rhs)
{
	assign(rhs.c_str());
    subString.clear();
	return *this;
}

CSTRING CSTRING::operator+ (CSTRING & rhs)
{
	CSTRING value;
	std::string *p1, *p2;
	p1 = this;
	p2 = &rhs;
	value = *p1 + *p2;
	return value;
}

CSTRING & CSTRING::operator+= (char ch)
{
	std::string *p1;
	p1 = this;
	*p1 += ch;
	return *this;
}

CSTRING & CSTRING::operator+= (std::string rhs)
{
	std::string *p1;
	p1 = this;
	*p1 += rhs;
	return *this;
}

CSTRING & CSTRING::operator+= (CSTRING & rhs)
{
	std::string *p1, *p2;
	p1 = this;
	p2 = &rhs;
	*p1 += *p2;
	return *this;
}

char & CSTRING::operator[] (int offset)
{
	if (offset > GetLen()) {__ERRORTHROW(ERR_MEMBDR); }
	else return at(offset);
}

CSTRING CSTRING::operator() (int offset) 
{
	if (offset >= (int)subString.size()) {__ERRORTHROW(ERR_MEMBDR); }
	else {
		CSTRING value;
		value = subString[offset];
		return value;
	}
}

