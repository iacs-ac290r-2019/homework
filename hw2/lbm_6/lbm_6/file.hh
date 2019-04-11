#ifndef FILE_HH
#define FILE_HH

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <limits>

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p);

#endif
