// ****
// *
// * A utility to write files containing the equilibrium ion populations
// * as a function of position
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../../../Radiation_Model/source/config.h"
#include "../../../../Radiation_Model/source/radiation.h"
#include "../../../source/file.h"

int main( int argc, char **argv )
{
PRADIATION pRadiation;
FILE *pCFGFile, *pAMRFile, *pPHYFile, *pIEFile;
char szResultsDirectory[256], szAMRFilename[256], szPHYFilename[256], szIEFilename[256];
double fs, fn, fT;
double fBuffer;
int iFrom, iTo, iNumCells;
int iBuffer, i, j;

if( argc == 1 ) {
	printf( "\nA configuration file must be specified. E.g. writeIEFile config.cfg\n");
	exit( EXIT_SUCCESS );
}

// Open and read the configuration file
pCFGFile = fopen( argv[1], "r" );
	// Get the directory containing the numerical results
	fscanf( pCFGFile, "%s", szResultsDirectory );
	// Get the range of output files over which to calculate the .ie files
	fscanf( pCFGFile, "%i", &iFrom );
	fscanf( pCFGFile, "%i", &iTo );
fclose( pCFGFile );

pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_neq.cfg" );

for( i=iFrom; i<=iTo; i++ )
{
	// Construct the profile .amr, .phy, and .ie filenames
	sprintf( szAMRFilename, "%s/profile%i.amr", szResultsDirectory, i );
	sprintf( szPHYFilename, "%s/profile%i.phy", szResultsDirectory, i );
	sprintf( szIEFilename, "%s/profile%i.ie", szResultsDirectory, i );	

	// Get the number of grid cells from the .amr file
	pAMRFile = fopen( szAMRFilename, "r" );
		ReadDouble( pAMRFile, &fBuffer );
		fscanf( pAMRFile, "%i", &iBuffer );
		ReadDouble( pAMRFile, &fBuffer );
		fscanf( pAMRFile, "%i", &iNumCells );
	fclose( pAMRFile );

	// Open the .phy and .ie files
	pPHYFile = fopen( szPHYFilename, "r" );
	pIEFile = fopen( szIEFilename, "w" );

	for(j=0; j<iNumCells; j++)
	{
		ReadDouble( pPHYFile, &fs );		// Position
		fprintf( pIEFile, "%.8e", fs );
		ReadDouble( pPHYFile, &fBuffer );	// Velocity
		ReadDouble( pPHYFile, &fBuffer );	// Sound speed
        ReadDouble( pPHYFile, &fn );		// Electron number density
        ReadDouble( pPHYFile, &fBuffer );	// Hydrogen number density
		ReadDouble( pPHYFile, &fBuffer );	// Electron pressure
		ReadDouble( pPHYFile, &fBuffer );	// Hydrogen pressure
		ReadDouble( pPHYFile, &fT );		// Electron temperature
		ReadDouble( pPHYFile, &fBuffer );	// Hydrogen temperature
		ReadDouble( pPHYFile, &fBuffer );	// Electron thermal energy flux
		ReadDouble( pPHYFile, &fBuffer );	// Hydrogen thermal energy flux
		
#ifdef DENSITY_DEPENDENT_RATES
		pRadiation->WriteEquilIonFracToFile( pIEFile, log10(fT), log10(fn) );
#else // DENSITY_DEPENDENT_RATES
		pRadiation->WriteEquilIonFracToFile( pIEFile, log10(fT) );
#endif // DENSITY_DEPENDENT_RATES
	}

	// Close the .phy and .ie files
	fclose( pIEFile );
	fclose( pPHYFile );
}

delete pRadiation;

return 0;
}