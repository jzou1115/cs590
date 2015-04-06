#include "argparser.hpp"

CARGPARSER::CARGPARSER(void){}

CARGPARSER::~CARGPARSER(void){}

bool CARGPARSER::Parse(int argc, char* argv[])
{
	int i, j, num;
	CSTRING flag;
	char ch;
	for(i=0;i<argc;i++) {
		flag = argv[i];
		if(flag[0] == '-') {
			if(flag.length() != 2 ) return false;
			ch = flag[1];
			if(!CheckAvArgs(ch)) return false;
			num = posNArgs[ch];
			cmdOpts[ch]=true;
			cmdMapOpts[ch] = "";
			for(j=0;j<num;j++) {
				if((i+1 > argc -1) || (argv[i+1][0] == '-')) return false;
				if(j==0) cmdMapOpts[ch] = argv[i+1];
				else cmdMapOpts[ch] = cmdMapOpts[ch] + " " + argv[i+1];
				i++;
			}
		} else {
			cmdElems.push_back(flag);
		}
	}
	return true;
}

void CARGPARSER::SetParam(CSTRING flag, CSTRING value)
{
	cmdMapOpts[flag[0]] = value;
}

void CARGPARSER::SetParam(char flag, CSTRING value)
{
	cmdMapOpts[flag] = value;
}

CSTRING & CARGPARSER::GetParam(CSTRING flag)
{
	return cmdMapOpts[flag[0]];
}

CSTRING & CARGPARSER::GetParam(char flag)
{
	return cmdMapOpts[flag];
}

int CARGPARSER::GetNCommands()
{
	return (int) cmdElems.size();
}

CSTRING & CARGPARSER::GetCommand(int nn)
{
	if(nn >= (int)cmdElems.size()) {
		cout << "\nERROR : offset exceeds range (CSTRING & CARGPARSER::GetCommand)" << endl;
		exit(-1);
	}
	return cmdElems[nn];
}

bool CARGPARSER::CheckArgs(CSTRING flags)
{
	int i;
	char ch;
	for(i=0;i<(int)flags.length();i++) {
		ch = flags[i];
		if(!cmdOpts[ch]) return false;
	}
	return true;
}

bool CARGPARSER::CheckArgs(char flag)
{
	return cmdOpts[flag];
}

bool CARGPARSER::InitArgs(CSTRING args, int num)
{
	int i;
	char ch;

	for(i=0;i<(int)args.length();i++) {
		ch = args[i];
		if((ch>='a' && ch<='z') || (ch>='A' && ch<='Z') ){
			posArgs[ch]=true;
			posNArgs[ch] = num;
		} else return false;
	}
	return true;
}

bool CARGPARSER::CheckAvArgs(CSTRING args)
{
	int i;
	char ch;
	for(i=0;i<(int)args.length();i++) {
		ch = args[i];
		if(!posArgs[ch]) return false;
	}
	return true;
}
bool CARGPARSER::CheckAvArgs(char arg)
{
	return posArgs[arg];
}
