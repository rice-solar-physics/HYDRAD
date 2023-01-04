// ****
// *
// * A routine to calculate the emission measure and the differential emission measure 
// * from HYDRAD .amr and .phy files
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../../../../HYDRAD/source/config.h"
#include "../../../source/file.h"

#define LOWER	0
#define UPPER	1

int main( int argc, char **argv )
{
FILE *pCFGFile, *pAMRFile, *pPHYFile, *pDEMFile;
char szResultsDirectory[256], szAMRFilename[256], szPHYFilename[256], szDEMFilename[256];
double *pflog10T, *pfDEM;
double fdexStep, fds, fne, fTe, flog10Te;
double fBuffer;
int *ppiFileSubRange[2];
int iNumBins, iNumOutputFiles, iFileRangeStart, iFileRangeEnd, iSubRangeStep, iNumGridCells, iBin;
int i, j, k, m;
int iBuffer;

if( argc == 1 ) {
	printf( "\nA configuration file must be specified. E.g. calculateDEM config.cfg\n");
	exit( EXIT_SUCCESS );
}

// Open the configuration file
printf( "\nUsing: %s\n", argv[1] );
pCFGFile = fopen( argv[1], "r" );
	// Get the directory containing the numerical results
	fscanf( pCFGFile, "%s", szResultsDirectory );
	// Get the number of temperature bins
	fscanf( pCFGFile, "%i", &iNumBins );
	// Get the bin width (dex)
	ReadDouble( pCFGFile, &fdexStep );
	// Allocate memory for the temperature intervals which define each bin
	pflog10T = (double*)malloc( sizeof(double*) * ( iNumBins + 1 ) );
	// Read the temperature intervals from the configuration file
	for( i=0; i<(iNumBins+1); i++ ) {
		ReadDouble( pCFGFile, &(pflog10T[i]) );
	}
	// Get the number of output files to write 
	// (the specified range is split into this number and summed over to calculate each DEM(T))
	fscanf( pCFGFile, "%i", &iNumOutputFiles );
	// Get the input file range
	fscanf( pCFGFile, "%i", &iFileRangeStart );
	fscanf( pCFGFile, "%i", &iFileRangeEnd );
// Close the configuration file
fclose( pCFGFile );

// Split the input file range to produce the specified number of output files
iSubRangeStep = ( iFileRangeEnd - iFileRangeStart ) / iNumOutputFiles;
// Allocate memory to the input file sub-ranges
ppiFileSubRange[LOWER] = (int*)malloc( sizeof(int*) * iNumOutputFiles );
ppiFileSubRange[UPPER] = (int*)malloc( sizeof(int*) * iNumOutputFiles );
// Define the input file sub-ranges
for( i=0; i<iNumOutputFiles; i++ ) {
	if( !i ) {
		ppiFileSubRange[LOWER][i] = iFileRangeStart + ( iSubRangeStep * i );
	} else {
		ppiFileSubRange[LOWER][i] = ppiFileSubRange[UPPER][i-1] + 1;
	}
	ppiFileSubRange[UPPER][i] = iFileRangeStart + ( iSubRangeStep * ( i + 1 ) );
}
// Ensure that the input file sub-ranges cover the full specified range
if( ppiFileSubRange[UPPER][i-1] < iFileRangeEnd ) {
	ppiFileSubRange[UPPER][i-1] = iFileRangeEnd;
}

printf( "\nTemperature bins:" );
for ( i=0; i<(iNumBins+1); i++ ) {
	printf( " %.2f", pflog10T[i] );
}
printf( "\n" );

printf( "\nNumber of output files: %i\n", iNumOutputFiles );
printf( "Input file sub-ranges:");
for( i=0; i<iNumOutputFiles-1; i++ ) {
	printf( " %i - %i;", ppiFileSubRange[LOWER][i], ppiFileSubRange[UPPER][i] );
}
printf( " %i - %i\n", ppiFileSubRange[LOWER][i], ppiFileSubRange[UPPER][i] );

// Allocate memory for the DEM(T) profile
pfDEM = (double*)malloc( sizeof(double) * iNumBins );

// Calculate the DEM(T) profiles for the input file sub-ranges
for( i=0; i<iNumOutputFiles; i++ ) {
	// Reset (zero) the memory allocated to the DEM(T) profile
	memset( pfDEM, '\0', iNumBins*sizeof(double) );
	for( j=ppiFileSubRange[LOWER][i]; j<=ppiFileSubRange[UPPER][i]; j++ ) {
		// Construct the .amr and .phy filenames
		sprintf( szAMRFilename, "%s/profile%i.amr", szResultsDirectory, j );
		sprintf( szPHYFilename, "%s/profile%i.phy", szResultsDirectory, j );
		// Open the .amr file and read the header information
		pAMRFile = fopen( szAMRFilename, "r" );
			// Read the header
			ReadDouble( pAMRFile, &fBuffer ); 			// Output time-stamp
			fscanf( pAMRFile, "%i", &iBuffer );			// File number (same as j)
			ReadDouble( pAMRFile, &fBuffer );			// Size of domain (loop length)
			fscanf( pAMRFile, "%i", &iNumGridCells );	// Number of grid cells
			// Open the .phy file
			pPHYFile = fopen( szPHYFilename, "r" );
				// Read the grid cell data from the .amr and .phy files and construct the DEM
				for( k=0; k<iNumGridCells; k++ ) {
					// AMR FILE //
					ReadDouble( pAMRFile, &fBuffer );	// Grid cell coordinate
					ReadDouble( pAMRFile, &fds );		// Grid cell size
					for( m=0; m<4; m++ ) {
						// Mass density, momentum density, electron energy, ion energy
						ReadDouble( pAMRFile, &fBuffer );
					}
					for( m=0; m<MAX_REFINEMENT_LEVEL+1; m++ ) {
						// Adaptive grid information
						fscanf( pAMRFile, "%i", &iBuffer );
					}
					// END OF AMR FILE //
					// PHY FILE //
					for( m=0; m<3; m++ ) {
						// Grid cell coordinate, bulk velocity, sound speed
						ReadDouble( pPHYFile, &fBuffer );
					}
					ReadDouble( pPHYFile, &fne );		// Electron density
					for( m=0; m<3; m++ ) {
						// Hydrogen density, electron pressure, hydrogen pressure
						ReadDouble( pPHYFile, &fBuffer );
					}
					ReadDouble( pPHYFile, &fTe );		// Electron temperature
					for( m=0; m<3; m++ ) {
						// Hydrogen temperature, electron heat flux, hydrogen heat flux
						ReadDouble( pPHYFile, &fBuffer );
					}
					// END OF PHY FILE //
					flog10Te = log10(fTe);
					if( flog10Te >= pflog10T[0] && flog10Te <= pflog10T[iNumBins] ) {
						// Find the temperature bin
						iBin = (int)( ( flog10Te - pflog10T[0] ) / fdexStep );
						// log(10) = 2.3025851
						pfDEM[iBin] += ( ( fne * fne ) / ( fTe * 2.3025851 ) ) * ( fds / fdexStep );
					}
				}
			// Close the .phy file
			fclose( pPHYFile );
		// Close the .amr file
		fclose( pAMRFile );
	}
	// Construct the output filename
	sprintf( szDEMFilename, "%s/DEM(%ito%i).txt", szResultsDirectory, ppiFileSubRange[LOWER][i], ppiFileSubRange[UPPER][i] );
	// Write the DEM(T) profile to the output file
	pDEMFile = fopen( szDEMFilename, "w" );
		fprintf( pDEMFile, "%i\n", iNumBins );
		for( j=0; j<iNumBins; j++ ) {
			fprintf( pDEMFile, "%.2f\t%f\n", pflog10T[j] + (fdexStep/2.0), log10(pfDEM[j]) );
		}
	fclose( pDEMFile );
	// Close the DEM file
}

// Free the memory allocated to the DEM(T) profile
free( pfDEM );

// Free the memory allocated to the input file sub-ranges
free( ppiFileSubRange[UPPER] );
free( ppiFileSubRange[LOWER] );

// Free the memory allocated to the temperature intervals
free( pflog10T );

return 0;
}