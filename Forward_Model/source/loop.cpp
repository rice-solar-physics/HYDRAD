// ****
// *
// * Loop class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/10/2017
// *
// ****


#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "loop.h"
#include "../../Resources/source/constants.h"


CLoop::CLoop( char *pszWorkingDirectory, int iFrom, int iTo, int iStep, double fdH )
{
Initialise( pszWorkingDirectory, iFrom, iTo, iStep, fdH );
}

CLoop::~CLoop( void )
{
FreeAll();
}

void CLoop::Initialise( char *pszWorkingDirectory, int iFrom, int iTo, int iStep, double fdH )
{
int i;
double fSD;

// Set the working directory
sprintf( szWorkingDirectory, "%s", pszWorkingDirectory );

// Calculate the number of strands
iNumStrands = ( ( iTo - iFrom ) / iStep ) + 1;

// Calculate the diameter of each strand
fSD = fdH / sqrt( (double)iNumStrands );

// Set the range of profiles to use as individual strands
iStrandRange[0] = iFrom;
iStrandRange[1] = iTo;
iStrandRange[2] = iStep;

// Allocate sufficient memory for the strands
ppStrand = (PPSTRAND)malloc( sizeof(CStrand*) * iNumStrands );

for( i=0; i<iNumStrands; i++ )
{
    printf( "Adding strand: " );
    ppStrand[i] = new CStrand( szWorkingDirectory, iFrom + ( i * iStep ), fSD );
}

printf( "\n" );
}

void CLoop::FreeAll( void )
{
int i;

for( i=0; i<iNumStrands; i++ )
    delete ppStrand[i];

free( ppStrand );
}

char *CLoop::pGetWorkingDirectory( void )
{
return szWorkingDirectory;
}

int CLoop::GetNumStrands( void )
{
return iNumStrands;
}

void CLoop::GetStrandRange( int *piStrandRange )
{
piStrandRange[0] = iStrandRange[0];
piStrandRange[1] = iStrandRange[1];
piStrandRange[2] = iStrandRange[2];
}

double CLoop::GetTimeStamp( int iStrand )
{
return ppStrand[iStrand]->GetTimeStamp();
}

int CLoop::GetNumCells( int iStrand )
{
return ppStrand[iStrand]->GetNumCells();
}

double CLoop::GetDiameter( int iStrand )
{
return ppStrand[iStrand]->GetDiameter();
}

double CLoop::GetLength( int iStrand )
{
return ppStrand[iStrand]->GetLength();
}

void CLoop::GetPHYData( int iStrand, int iCell, PPHYDATA pStrandPHYData )
{
ppStrand[iStrand]->GetPHYData( iCell, pStrandPHYData );
}

int CLoop::GetNumNEQElements( int iStrand )
{
return ppStrand[iStrand]->GetNumNEQElements();
}

int CLoop::GetAtomicNumber( int iStrand, int iElement )
{
return ppStrand[iStrand]->GetAtomicNumber( iElement );
}

double CLoop::GetNEQ( int iStrand, int iCell, int iElement, int iIon )
{
return ppStrand[iStrand]->GetNEQ( iCell, iElement, iIon );
}