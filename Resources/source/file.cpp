// ****
// *
// * A simple routine to accurately read floating-point values
// * from datafiles
// * 
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 15/04/2005
// *
// ****


#include <stdio.h>
#include <stdlib.h>


void ReadDouble( void *pInputFile, double *pDouble )
{
FILE *pFile;
char buffer[256], *end;

pFile = (FILE*)pInputFile;

fscanf( pFile, "%s", buffer );

*pDouble = strtod( buffer, &end );
}
