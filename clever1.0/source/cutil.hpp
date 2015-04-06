#pragma once
#include "general.hpp"

class CUTIL
{
public:
	CUTIL(void);
	~CUTIL(void);
	FILE* SafeOpen(const char *fname, const char *mode);

};
