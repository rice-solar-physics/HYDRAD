// ****
// *
// * A routine to calculate the emission measure and the differential emission measure 
// * from spatially averaged electron densities and temperatures
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/11/2015
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../../source/file.h"

int main(void)
{
FILE *pCONFIGFile, *pINPUTFile, *pOUTPUTFile;
char szINPUTFilename[256], szOUTPUTFilename[256];
char szBuffer[256];
double *pfEM;
double fdR, fLogTmax, fLogTmin, fdexStep;
double fTimeStamp, fn, fTe, fTH;
double fLogTe, fEM_T, fDEM_T;
int iNumElements, iNumFiles, iNumRows, iBin;
int i, j;

// Open the configuration file
pCONFIGFile = fopen( "config.cfg", "r" );

	// Get the line-of-sight depth
	ReadDouble( pCONFIGFile, &fdR ); fscanf( pCONFIGFile, "%s", szBuffer );
	// Get the temperature range of the emission measure
	ReadDouble( pCONFIGFile, &fLogTmin ); fscanf( pCONFIGFile, "%s", szBuffer );
	ReadDouble( pCONFIGFile, &fLogTmax ); fscanf( pCONFIGFile, "%s", szBuffer );
	ReadDouble( pCONFIGFile, &fdexStep ); fscanf( pCONFIGFile, "%s", szBuffer );
	iNumElements = ( fLogTmax - fLogTmin ) / fdexStep;
	pfEM = (double*)malloc( sizeof(double) * iNumElements );

	// Get the number of input files
	fscanf( pCONFIGFile, "%i", &iNumFiles ); fscanf( pCONFIGFile, "%s", szBuffer );
	for( i=0; i<iNumFiles; i++ )
	{
		// Reset the emission measure
		for( j=0; j<iNumElements; j++ )
			pfEM[j] = 0.0;

		// Get the filename of the input file
		fscanf( pCONFIGFile, "%s", szINPUTFilename ); fscanf( pCONFIGFile, "%s", szBuffer );
		// Open the input file
		pINPUTFile = fopen( szINPUTFilename, "r" );

			// Get the number of rows from the input file
			fscanf( pINPUTFile, "%i", &iNumRows );
			for( j=0; j<iNumRows; j++ )
			{
				// Get the timestamp
				ReadDouble( pINPUTFile, &fTimeStamp );
				// Get the number density
				ReadDouble( pINPUTFile, &fn );
				// Get the electron temperature
				ReadDouble( pINPUTFile, &fTe );
				// Get the hydrogen temperature
				ReadDouble( pINPUTFile, &fTH );

				// Find the temperature bin to fill with the emission measure
				fLogTe = log10( fTe );
				iBin = ( fLogTe - fLogTmin ) / fdexStep;
				if( iBin >= 0 && iBin < iNumElements )
				pfEM[iBin] += fn * fn;
			}

		fclose( pINPUTFile );

		// Get the filename of the output file
		fscanf( pCONFIGFile, "%s", szOUTPUTFilename ); fscanf( pCONFIGFile, "%s", szBuffer );
		// Open the output file
		pOUTPUTFile = fopen( szOUTPUTFilename, "w" );
			fLogTe = fLogTmin + (fdexStep/2.0);
			for( j=0; j<iNumElements; j++ )
			{
				fEM_T = ((fdR*pfEM[j])/(double)iNumRows)/fdexStep;
				fDEM_T = fEM_T / ( log(10) * pow(10.0,fLogTe) );
				fprintf( pOUTPUTFile, "%.4g\t%.4g\t%.4g\n", fLogTe, fEM_T, fDEM_T );
				fLogTe += fdexStep;
			}
		fclose( pOUTPUTFile );
	}

	free( pfEM );

fclose( pCONFIGFile );

return 0;
}