#pragma once
#include "cstring.hpp"
#include <vector>
#include <map>

class CARGPARSER
{
/*******************************************************************************************
   1. Initialize arguments
	   Possible arguments, {a...z, A...Z}
	   InitArgs("adefg",0);            posArgs[0...]=[a,d,e,f,g,q,h,t,r,c,s]
	   InitArgs("qhtr",1);     ===>    posNArgs({a,d,q,c})={0,0,1,2}
	   InitArgs("cs",2);
   2. Parse a command line into variables
	   <command> -a -d -q ad.txt -c 1 d de.1 de.2
	   cmdOpts[{a,d,q,f...}]={true,true,true,false,...}
	   cmdMapOpts({a,d,q,c})={"","","ad.txt","1 d"}
	   cmdElems[:]={"de.1","de.2"}
*******************************************************************************************/
private:
	std::map<char,bool>       cmdOpts;  // save arguments which are appeared in a command line
	std::vector<CSTRING>     cmdElems;  // save parameters except arguments from a command line
	std::map<char,CSTRING> cmdMapOpts;  // save a map containing (argument => its value)
	std::map<char,int>       posNArgs;  // required number of parameters for each argument
	std::map<char,bool>       posArgs;  // a list of possible arguments.

public:
	CARGPARSER(void);
	~CARGPARSER(void);
	CARGPARSER(int argc, char* argv[]); 
	bool Parse(int argc, char* argv[]);  // analyze a command line and allocate each value to proper variables
	                                     // fail in analysis return false and else true
	bool CheckArgs(CSTRING);             // check whether some arguments which are appeared in a command line
	                                     // return true when all arguments are appeared and else false
	bool CheckArgs(char);             
	CSTRING & GetParam(CSTRING);         // get argument value
	CSTRING & GetParam(char);          
	void SetParam(CSTRING, CSTRING);         // get argument value
	void SetParam(char, CSTRING);          
	int GetNCommands();
	CSTRING & GetCommand(int nn);
	bool InitArgs(CSTRING args, int num); // setting posArgs and posNArgs;
	bool CheckAvArgs(CSTRING args);      // Check available argments
	bool CheckAvArgs(char arg); 
};
