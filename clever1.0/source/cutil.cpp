#include "cutil.hpp"

CUTIL::CUTIL(void)
{
}

CUTIL::~CUTIL(void)
{
}

/* Opens the specified file using the specified mode. Exits with an error on failure */
FILE* CUTIL::SafeOpen(const char *fname, const char *mode) {
	FILE * fp;
	fp = fopen(fname,mode);
	if (fp == NULL) { 
		printf("Error: failed in opening a file '%s'\n", fname);
		exit(-1);
	}
	return fp;
}

