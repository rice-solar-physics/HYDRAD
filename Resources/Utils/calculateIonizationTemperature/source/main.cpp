// ****
// *
// * A utility to write files containing the (effective) ionization temperature as a function of position
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../../../Radiation_Model/source/ionfrac.h"
#include "../../../source/file.h"

int main( int argc, char **argv )
{
// Functions
double fionT( PRADIATION pRadiation, PIONFRAC pIonFrac, double *pflogTRange );
// Variables
PRADIATION pRadiation;
PIONFRAC pIonFrac;
FILE *pCFGFile, *pAMRFile, *pPHYFile, *pINEFile, *pionTFile;
char szResultsDirectory[256], szAMRFilename[256], szPHYFilename[256], szINEFilename[256], szionTFilename[256];
double flogTRange[3], fs, fne, fTe, fTH;
double fBuffer;
int iFrom, iTo, iNumCells;
int iBuffer, i, j;

if( argc == 1 ) {
	printf( "\nA configuration file must be specified. E.g. calculateIonizationTemperature config.cfg\n");
	exit( EXIT_SUCCESS );
}

// Open and read the configuration file
pCFGFile = fopen( argv[1], "r" );
	// Get the directory containing the numerical results
	fscanf( pCFGFile, "%s", szResultsDirectory );
	// Get the (log) temperature range and increment (dex) over which to calculate the (effective) ionization temperature
	for( i=0; i<3; i++ )
		ReadDouble( pCFGFile, &(flogTRange[i]) );
	// Get the range of output files over which to calculate the (effective) ionization temperature
	fscanf( pCFGFile, "%i", &iFrom );
	fscanf( pCFGFile, "%i", &iTo );
fclose( pCFGFile );

pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_neq.cfg" );
pIonFrac = new CIonFrac( NULL, (char *)"Radiation_Model/config/elements_neq.cfg", pRadiation );

for( i=iFrom; i<=iTo; i++ )
{
	// Construct the profile .amr, .phy, .ine, and .ionT filenames
	sprintf( szAMRFilename, "%s/profile%i.amr", szResultsDirectory, i );
	sprintf( szPHYFilename, "%s/profile%i.phy", szResultsDirectory, i );
	sprintf( szINEFilename, "%s/profile%i.ine", szResultsDirectory, i );
	sprintf( szionTFilename, "%s/profile%i.ionT", szResultsDirectory, i );

	// Get the number of grid cells from the .amr file
	pAMRFile = fopen( szAMRFilename, "r" );
		ReadDouble( pAMRFile, &fBuffer );
		fscanf( pAMRFile, "%i", &iBuffer );
		ReadDouble( pAMRFile, &fBuffer );
		fscanf( pAMRFile, "%i", &iNumCells );
	fclose( pAMRFile );

	// Open the .phy, .ine, and ionT files
	pPHYFile = fopen( szPHYFilename, "r" );
	pINEFile = fopen( szINEFilename, "r" );
	pionTFile = fopen( szionTFilename, "w" );

	for(j=0; j<iNumCells; j++)
	{
		// Read data from the .phy file
		ReadDouble( pPHYFile, &fs );		// Position
		ReadDouble( pPHYFile, &fBuffer );	// Velocity
		ReadDouble( pPHYFile, &fBuffer );	// Sound speed
        ReadDouble( pPHYFile, &fne );		// Electron number density
        ReadDouble( pPHYFile, &fBuffer );	// Hydrogen number density
		ReadDouble( pPHYFile, &fBuffer );	// Electron pressure
		ReadDouble( pPHYFile, &fBuffer );	// Hydrogen pressure
		ReadDouble( pPHYFile, &fTe );		// Electron temperature
		ReadDouble( pPHYFile, &fTH );		// Hydrogen temperature
		ReadDouble( pPHYFile, &fBuffer );	// Electron thermal energy flux
		ReadDouble( pPHYFile, &fBuffer );	// Hydrogen thermal energy flux

		// Read data from the .ine file
        ReadDouble( pINEFile, &fBuffer );
        pIonFrac->ReadAllIonFracFromFile( pINEFile );

		// Write data to the .ionT file
		fprintf( pionTFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", fs, fne, fTe, fTH, fionT( pRadiation, pIonFrac, flogTRange ) );
	}

	// Close the .phy, .ine, and ionT files
	fclose( pionTFile );
	fclose( pINEFile );
	fclose( pPHYFile );
}

delete pIonFrac;
delete pRadiation;

return 0;
}

double fionT( PRADIATION pRadiation, PIONFRAC pIonFrac, double *pflogTRange )
{
double *pequilState, *pnonequilState, flog10_Teff, fminDiff, flog10_T, fDiff;
int *piA, iNumElements, iElement, iZ;

piA = pRadiation->pGetAtomicNumbers( &iNumElements );

flog10_Teff = pflogTRange[0];
fminDiff = 1e300;
for( flog10_T=pflogTRange[0]; flog10_T<=pflogTRange[1]; flog10_T+=pflogTRange[2] ) {
	fDiff = 0.0;
	for( iElement=0; iElement<iNumElements; iElement++ ) {
		pequilState = (double*)malloc( sizeof(double) * ( piA[iElement] + 1 ) );
		pRadiation->GetEquilIonFrac( piA[iElement], pequilState, flog10_T );
		pnonequilState = pIonFrac->pGetIonFrac( piA[iElement] );

		for( iZ=0; iZ<=piA[iElement]; iZ++ )
			fDiff+= fabs( pequilState[iZ] - pnonequilState[iZ] );

		free( pequilState );
	}

	if( fDiff < fminDiff ) {
		fminDiff = fDiff;
		flog10_Teff = flog10_T;
	}	
}

return pow( 10.0, flog10_Teff );
}