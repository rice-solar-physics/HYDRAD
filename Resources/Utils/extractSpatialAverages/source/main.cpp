// ****
// *
// * A routine to extract spatially averaged quantities 
// * from a series of data files of user-specified format
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/12/2015
// *
// ****

#include <stdio.h>
#include <malloc.h>

#include "../../../../HYDRAD/source/config.h"
#include "../../../source/file.h"


int main(void)
{
FILE *pCONFIGFile, *pAMRFile, *pEXTFile, *pOUTPUTFile;
char szRoot[256], szExtension[32], szAMRFilename[256], szEXTFilename[256], szOUTPUTFilename[256];
char szBuffer[256];
double *pfSum;
double fLLp, fULp, fLL, fUL;
double fTimeStamp, fL, fs, fds;
double fBuffer;
int iRange[2], *piC;
int iNumColumns, iNumAverages, iNumFiles, iNumCells;
int i, j, k, l, m, n;
int iBuffer;

// Open the configuration file
pCONFIGFile = fopen( "config.cfg", "r" );

	// Get the filename structure of the files containing the data to be spatially averaged
	fscanf(pCONFIGFile, "%s", szRoot ); fscanf( pCONFIGFile, "%s", szBuffer );
	fscanf(pCONFIGFile, "%s", szExtension ); fscanf( pCONFIGFile, "%s", szBuffer );
	// Get the number of columns in the files
	fscanf(pCONFIGFile, "%i", &iNumColumns ); fscanf( pCONFIGFile, "%s", szBuffer );
	
	// Get the number of columns to calculate spatial averages for
	fscanf(pCONFIGFile, "%i", &iNumAverages ); fscanf( pCONFIGFile, "%s", szBuffer );
	// Get the particular column numbers to calculate spatial averages for
	piC = (int*)malloc( sizeof(int) * iNumAverages );
	for( i=0; i<iNumAverages; i++ )
		fscanf( pCONFIGFile, "%i", &(piC[i]) );
	fscanf( pCONFIGFile, "%s", szBuffer );
	pfSum = (double*)malloc( sizeof(double) * iNumAverages );
	
	// Get the lower and upper spatial limits as percentages and convert to fractions
	ReadDouble( pCONFIGFile, &fLLp ); fscanf( pCONFIGFile, "%s", szBuffer );
	ReadDouble( pCONFIGFile, &fULp ); fscanf( pCONFIGFile, "%s", szBuffer );
	fLLp /= 100.0; fULp /= 100.0;

	// Get the number of spatial average files to write
	fscanf(pCONFIGFile, "%i", &iNumFiles ); fscanf( pCONFIGFile, "%s", szBuffer );

	for( i=0; i<iNumFiles; i++ )
	{
		// Get the file range
		fscanf( pCONFIGFile, "%i", &(iRange[0]) );
		fscanf( pCONFIGFile, "%i", &(iRange[1]) ); fscanf( pCONFIGFile, "%s", szBuffer );

		// Construct the output filename and open the file
		sprintf( szOUTPUTFilename, "f(t)(%ito%i).txt", iRange[0], iRange[1] );
		pOUTPUTFile = fopen( szOUTPUTFilename, "w" );
			// Write the number of rows to the output file
			fprintf( pOUTPUTFile, "%i\n", (iRange[1]-iRange[0])+1 );
			for( j=iRange[0]; j<=iRange[1]; j++ )
			{
				// Reset the sum totals
				for( k=0; k<iNumAverages; k++ )
					pfSum[k] = 0.0;

				// Construct the .amr filename and open the file
				sprintf( szAMRFilename, "profile%i.amr", j );
				pAMRFile = fopen( szAMRFilename, "r" );
					// Get the timestamp from the .amr file
					ReadDouble( pAMRFile, &fTimeStamp );
					// Get the file number from the .amr file
					fscanf( pAMRFile, "%i", &iBuffer );
					// Get the loop length from the .amr file
					ReadDouble( pAMRFile, &fL );
					fLL = fLLp * fL; fUL = fULp * fL;
					// Get the number of grid cells from the .amr file
					fscanf( pAMRFile, "%i", &iNumCells );

					// Construct the .ext filename and open the file
					sprintf( szEXTFilename, "%s%i.%s", szRoot, j, szExtension );
					pEXTFile = fopen( szEXTFilename, "r" );
						for( k=0; k<iNumCells; k++ )
						{
							// Read the data from the .amr file
							ReadDouble( pAMRFile, &fs );		// Position
							ReadDouble( pAMRFile, &fds );		// Grid cell width
							ReadDouble( pAMRFile, &fBuffer );	// Mass density
							ReadDouble( pAMRFile, &fBuffer );	// Momentum density
							ReadDouble( pAMRFile, &fBuffer );	// Electron energy
							ReadDouble( pAMRFile, &fBuffer );	// Hydrogen energy
							for( l=0; l<=MAX_REFINEMENT_LEVEL; l++ )
								fscanf( pAMRFile, "%i", &iBuffer );


							// Read the data from the .ext file
							for( m=0; m<iNumColumns; m++ )
							{
								ReadDouble( pEXTFile, &fBuffer );
								// Check to see whether the current column is one to be spatially averaged
								for( n=0; n<iNumAverages; n++ )
									if( ( m+1 == piC[n] ) && ( fs >= fLL && fs <= fUL ) )
									{
										// Spatially average this column
										pfSum[n] += fBuffer * fds;
										break;
									}
							}
						}
					fclose( pEXTFile );
				fclose( pAMRFile );

				for( k=0; k<iNumAverages; k++ )
					pfSum[k] /= ( fUL - fLL );

				// Write the data to the output file
				fprintf( pOUTPUTFile, "%.4g", fTimeStamp );
				for( k=0; k<iNumAverages; k++ )
					fprintf( pOUTPUTFile, "\t%.4g", pfSum[k] );
				fprintf( pOUTPUTFile, "\n" );
			}
		fclose( pOUTPUTFile );
	}
fclose( pCONFIGFile );	

free( pfSum );
free( piC );

return 0;
}