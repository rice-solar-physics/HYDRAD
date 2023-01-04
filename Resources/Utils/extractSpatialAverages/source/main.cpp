// ****
// *
// * A routine to extract spatially averaged quantities 
// * from a series of data files of user-specified format
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>

#include "../../../../HYDRAD/source/config.h"
#include "../../../source/file.h"

// #define READ_ELECTRON_MASS_DENSITY		// Required for runs using the optically-thick chromosphere model

int main( int argc, char **argv )
{
FILE *pCONFIGFile, *pAMRFile, *pEXTFile, *pOUTPUTFile;
char szResultsDirectory[256], szRoot[256], szExtension[32], szAMRFilename[256], szEXTFilename[256], szOUTPUTFilename[256];
double *pfSum;
double fLLp, fULp, fLL, fUL;
double fTimeStamp, fL, fs, fds;
double fBuffer;
int iRange[2], *piC;
int iNumColumns, iNumAverages, iNumFiles, iNumCells;
int i, j, k, l, m, n;
int iBuffer;

if( argc == 1 ) {
	printf( "\nA configuration file must be specified. E.g. extractSpatialAverages config.cfg\n");
	exit( EXIT_SUCCESS );
}

// Open the configuration file
pCONFIGFile = fopen( argv[1], "r" );
	// Get the directory containing the numerical results
	fscanf( pCONFIGFile, "%s", szResultsDirectory );
	// Get the filename structure of the files containing the data to be spatially averaged
	fscanf(pCONFIGFile, "%s", szRoot );
	fscanf(pCONFIGFile, "%s", szExtension );
	// Get the number of columns in the files
	fscanf(pCONFIGFile, "%i", &iNumColumns );
	
	// Get the number of columns to calculate spatial averages for
	fscanf(pCONFIGFile, "%i", &iNumAverages );
	// Get the particular column numbers to calculate spatial averages for
	piC = (int*)malloc( sizeof(int) * iNumAverages );
	for( i=0; i<iNumAverages; i++ )
		fscanf( pCONFIGFile, "%i", &(piC[i]) );
	pfSum = (double*)malloc( sizeof(double) * iNumAverages );
	
	// Get the lower and upper spatial limits as percentages and convert to fractions
	ReadDouble( pCONFIGFile, &fLLp );
	ReadDouble( pCONFIGFile, &fULp );
	fLLp /= 100.0; fULp /= 100.0;

	// Get the number of spatial average files to write
	fscanf( pCONFIGFile, "%i", &iNumFiles );

	printf( "\n%s/%s.%s\n", szResultsDirectory, szRoot, szExtension );
	printf( "\nNumber of Columns = %i\n", iNumColumns );
	printf( "\nNumber of Columns to Average = %i\n", iNumAverages );	
	printf( "\tColumns:");
	for( i=0; i<iNumAverages; i++ ) {
		printf( " %i", piC[i] );
	}
	printf( "\n" );
	printf( "\nAverage Between %g and %g of Domain\n", fLLp, fULp );
	printf( "\nNumber of Spatial Average Files = %i\n", iNumFiles );
	printf( "\tFile Number Range:" );

	for( i=0; i<iNumFiles; i++ )
	{
		// Get the file range
		fscanf( pCONFIGFile, "%i", &(iRange[0]) );
		fscanf( pCONFIGFile, "%i", &(iRange[1]) );

		printf( " [%i,%i]", iRange[0], iRange[1] );

		// Construct the output filename and open the file
		sprintf( szOUTPUTFilename, "%s/f(t)(%ito%i).txt", szResultsDirectory, iRange[0], iRange[1] );
		pOUTPUTFile = fopen( szOUTPUTFilename, "w" );
			// Write the number of rows to the output file
			fprintf( pOUTPUTFile, "%i\n", (iRange[1]-iRange[0])+1 );
			for( j=iRange[0]; j<=iRange[1]; j++ )
			{
				// Reset the sum totals
				for( k=0; k<iNumAverages; k++ )
					pfSum[k] = 0.0;

				// Construct the .amr filename and open the file
				sprintf( szAMRFilename, "%s/profile%i.amr", szResultsDirectory, j );
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
					sprintf( szEXTFilename, "%s/%s%i.%s", szResultsDirectory, szRoot, j, szExtension );
					pEXTFile = fopen( szEXTFilename, "r" );
						for( k=0; k<iNumCells; k++ )
						{
							// Read the data from the .amr file
							ReadDouble( pAMRFile, &fs );		// Position
							ReadDouble( pAMRFile, &fds );		// Grid cell width
#ifdef READ_ELECTRON_MASS_DENSITY
							ReadDouble( pAMRFile, &fBuffer );	// Electron mass density
#endif // READ_ELECTRON_MASS_DENSITY
							ReadDouble( pAMRFile, &fBuffer );	// Hydrogen mass density							
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

printf( "\n" );

return 0;
}