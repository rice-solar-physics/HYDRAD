// ****
// *
// * Strand class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 08/06/2015
// *
// ****


#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "strand.h"
#include "../../HYDRAD/source/config.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/constants.h"


CStrand::CStrand( char *pszWorkingDirectory, int iProfileNumber, double fH )
{
Initialise( pszWorkingDirectory, iProfileNumber, fH );
}

CStrand::~CStrand( void )
{
FreeAll();
}

void CStrand::Initialise( char *pszWorkingDirectory, int iProfileNumber, double fH )
{
FILE *pAMRFile, *pPHYFile, *pINEFile, *pNEQFile;
char szAMRFilename[256], szPHYFilename[256], szINEFilename[256], szBuffer[256];
double fTemp, fBuffer;
int iBuffer, i, j, k;

// Set the strand diameter
fSD = fH;

// Construct the .amr filename
sprintf( szAMRFilename, "%s/profile%i.amr", pszWorkingDirectory, iProfileNumber );
printf( "%s; ", szAMRFilename );

// Construct the .phy filename
sprintf( szPHYFilename, "%s/profile%i.phy", pszWorkingDirectory, iProfileNumber );

// Open the .amr file
pAMRFile = fopen( szAMRFilename, "r" );

// Retrieve the header information
ReadDouble( pAMRFile, &fTimeStamp );
fscanf( pAMRFile, "%i", &iBuffer );
ReadDouble( pAMRFile, &fLength );
fTemp = _PI_ / fLength;
fscanf( pAMRFile, "%i", &iNumCells );

// Allocate sufficient memory to hold the strand data
pPHYData = (PPHYDATA)malloc( sizeof(SPHYData) * iNumCells );

// Open the .phy file
pPHYFile = fopen( szPHYFilename, "r" );

for( i=0; i<iNumCells; i++ )
{
    // Retrieve data from the .amr file
    ReadDouble( pAMRFile, &(pPHYData[i].fs) );
    ReadDouble( pAMRFile, &(pPHYData[i].fds) );
    ReadDouble( pAMRFile, &fBuffer );
    ReadDouble( pAMRFile, &fBuffer );
    ReadDouble( pAMRFile, &fBuffer );
    ReadDouble( pAMRFile, &fBuffer );
    for( j=0; j<=MAX_REFINEMENT_LEVEL; j++ )
        fscanf( pAMRFile, "%i", &iBuffer );

    // Retrieve data from the .phy file
    ReadDouble( pPHYFile, &fBuffer );           // Field-aligned coordinate
    ReadDouble( pPHYFile, &(pPHYData[i].fv) );  // Velocity
    // Project the velocity onto the line-of-sight
    pPHYData[i].fvp = - pPHYData[i].fv * cos( fTemp * pPHYData[i].fs );
    ReadDouble( pPHYFile, &fBuffer );           // Sound speed
    ReadDouble( pPHYFile, &(pPHYData[i].fn) );  // Electron density
    ReadDouble( pPHYFile, &fBuffer );           // Hydrogen density
    ReadDouble( pPHYFile, &fBuffer );           // Electron pressure
    ReadDouble( pPHYFile, &fBuffer );           // Hydrogen pressure
    ReadDouble( pPHYFile, &(pPHYData[i].fTe) ); // Electron temperature
    ReadDouble( pPHYFile, &(pPHYData[i].fTi) ); // Hydrogen temperature
    ReadDouble( pPHYFile, &fBuffer );           // Electron heat flux
    ReadDouble( pPHYFile, &fBuffer );           // Hydrogen heat flux
}

fclose( pPHYFile );

fclose( pAMRFile );

// Check for the existence of .ine files

// Construct the filename
sprintf( szINEFilename, "%s/profile%i.ine", pszWorkingDirectory, iProfileNumber );

pINEFile = NULL;
pINEFile = fopen( szINEFilename, "r" );

if( pINEFile )
{
    printf( "%s\n", szINEFilename );

    // Get the number of elements solved for
    pNEQFile = fopen( "Radiation_Model/config/elements_neq.cfg", "r" );
	
    fscanf( pNEQFile, "%s", szBuffer );
    fscanf( pNEQFile, "%s", szBuffer );
    fscanf( pNEQFile, "%s", szBuffer );
    fscanf( pNEQFile, "%s", szBuffer );
	
    fscanf( pNEQFile, "%i", &iNumNEQElements );
    // Allocate sufficient memory to hold the list of atomic numbers
    piNEQ_Z = (int*)malloc( sizeof(int) * iNumNEQElements );
    for( i=0; i<iNumNEQElements; i++ )
    {
	fscanf( pNEQFile, "%s", szBuffer );
	fscanf( pNEQFile, "%i", &(piNEQ_Z[i]) );
    }
	
    fclose( pNEQFile );

    // Allocate sufficient memory to hold the ion population fractions for each cell
    pppfNEQ = (double***)malloc( sizeof(double) * iNumCells );
    for( i=0; i<iNumCells; i++ )
    {
        pppfNEQ[i] = (double**)malloc( sizeof(double) * iNumNEQElements );
	for( j=0; j<iNumNEQElements; j++ )
            pppfNEQ[i][j] = (double*)malloc( sizeof(double) * (piNEQ_Z[j] + 1) );
    }

    // Retrieve the data from the .ine file
    for( i=0; i<iNumCells; i++ )
    {
	// Get the field-aligned coordinate
	ReadDouble( pINEFile, &fBuffer );
	for( j=0; j<iNumNEQElements; j++ )
	{
            // Get the atomic number
            fscanf( pINEFile, "%i", &iBuffer );
            for( k=0; k<=piNEQ_Z[j]; k++ )
            {
                // Get the ion population fraction
		ReadDouble( pINEFile, &(pppfNEQ[i][j][k]) );
            }
	}
    }

    fclose( pINEFile );
}
else
{
    printf( "%s does not exist\n", szINEFilename );
    iNumNEQElements = 0;
    piNEQ_Z = NULL;
    pppfNEQ = NULL;
}
}

void CStrand::FreeAll( void )
{
int i, j;

if( iNumNEQElements > 0 )
{
    for( i=0; i<iNumCells; i++ )
    {
        for( j=0; j<iNumNEQElements; j++ )
            free( pppfNEQ[i][j] );

        free( pppfNEQ[i] );
    }

    free( pppfNEQ );
    free( piNEQ_Z );
}

free( pPHYData );
}

double CStrand::GetTimeStamp( void )
{
return fTimeStamp;
}

int CStrand::GetNumCells( void )
{
return iNumCells;
}

double CStrand::GetDiameter( void )
{
return fSD;
}

double CStrand::GetLength( void )
{
return fLength;
}

void CStrand::GetPHYData( int iCell, PPHYDATA pStrandPHYData )
{
memcpy( pStrandPHYData, &(pPHYData[iCell]), sizeof(SPHYData) );
}

int CStrand::GetNumNEQElements( void )
{
return iNumNEQElements;
}

int CStrand::GetAtomicNumber( int iElement )
{
return piNEQ_Z[iElement];
}

double CStrand::GetNEQ( int iCell, int iElement, int iIon )
{
return pppfNEQ[iCell][iElement][iIon];
}