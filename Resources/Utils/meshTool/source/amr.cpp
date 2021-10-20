// ****
// *
// * Function bodies for the HYDRAD .amr file reader
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/20/2021
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "amr.h"
#include "../../../source/file.h"
#include "../../../source/fitpoly.h"

// **** AMR FILE CLASS ****

// Constructor
CAMRFile::CAMRFile( void )
{
}

// Constructor
CAMRFile::CAMRFile( char *pszAMRFilename, int iAMRMaxRL )
{
	iMaxRL = iAMRMaxRL;
	if( !ReadAMRFile( pszAMRFilename ) )
		printf( "\nFailed to ReadAMRFile()\n" );
}

// Destructor
CAMRFile::~CAMRFile( void )
{
	if( !FreeAll() )
		printf( "\nFailed to FreeAll() in ~CAMRFile()\n" );
}

bool CAMRFile::AllocateAMRQuantities( void )
{
	int i;
	// Allocate memory to store the .amr file data
	ppfAMRQuantities = (double**)malloc( sizeof(double*) * AMRFileHeader.iNumberOfCells );
	if( !ppfAMRQuantities ) {
		printf( "\nFailed to allocate memory at ppfAMRQuantities\n" );
		return( false );
	}
	for( i=0; i<AMRFileHeader.iNumberOfCells; i++ ) {
		// Allocate memory to store the .amr file data
		ppfAMRQuantities[i] = (double*)malloc( sizeof(double) * iNumberOfColumns );
		if( !ppfAMRQuantities[i] ) {
			printf( "\nFailed to allocate memory at ppfAMRQuantities[%i]\n", i );
			return( false );
		}
	}
	return( true );
}

bool CAMRFile::AllocateUniqueIDStructure( void )
{
	int i;
	// Allocate memory to store the unique ID structure used by the adaptive mesh algorithm to track cell connectivity
	ppiUniqueIDStructure = (int**)malloc( sizeof(int*) * AMRFileHeader.iNumberOfCells );
	if( !ppiUniqueIDStructure ) {
		printf( "\nFailed to allocate memory at ppiUniqueIDStructure\n" );
		return( false );
	}
	for( i=0; i<AMRFileHeader.iNumberOfCells; i++ ) {
		ppiUniqueIDStructure[i] = (int*)malloc( sizeof(int) * (iMaxRL+1) );
		if( !ppiUniqueIDStructure[i] ) {
			printf( "\nFailed to allocate memory at ppiUniqueIDStructure[%i]\n", i );
			return( false );
		}
	}
	return( true );
}

bool CAMRFile::GenerateUniqueIDStructure( void )
{
	int iRL, iUniqueID, iGroupSize;
	int h, i, j;

	// Set the first column of the unique ID structure to the maximum refinement level
	j = 0;
	for( i=0; i<AMRFileHeader.iNumberOfCells; i++ )
		ppiUniqueIDStructure[i][j] = iMaxRL;

	// Working row-by-row and then column-by-column, set the unique ID numbers for each refinement level and group of connected cells
	iRL = iMaxRL;
	iUniqueID = 0;
	for( j=1; j<=iMaxRL; j++ ) {
		iGroupSize = pow(2,iRL);
		h = 0;
		for( i=0; i<AMRFileHeader.iNumberOfCells; i++ )
		{
			ppiUniqueIDStructure[i][j] = iUniqueID;
			h++;
			if( h == iGroupSize ) {
				iUniqueID++;
				h = 0;
			}
		}
		iRL--;
	}
	return( true );
}

bool CAMRFile::ReadAMRFile( char *pszAMRFilename )
{
	FILE *pAMRFile;
	int i, j, iBuffer;

	// Check to see whether a .amr file has already been opened or defined
	if( ppfAMRQuantities || ppiUniqueIDStructure ) {
		printf( "\nCannot overwrite protected data. Create a new instance of the AMRFile object instead\n" );
		return( false );
	}

	// Save an internal copy of the .amr filename
	strcpy( szAMRFilename, pszAMRFilename );

	// Open the .amr file
	pAMRFile = fopen( szAMRFilename, "r" );
		if( !pAMRFile ) {
			printf ( "\nFailed to open .amr file: %s\n", szAMRFilename );
			return( false );
		}
		// Get the profile time
		ReadDouble( pAMRFile, &(AMRFileHeader.fProfileTime) );
		// Get the profile number
		fscanf( pAMRFile, "%i", &(AMRFileHeader.iProfileNumber) );
		// Get the field line length
		ReadDouble( pAMRFile, &(AMRFileHeader.fFullLength) );
		// Get the number of cells
		fscanf( pAMRFile, "%i", &(AMRFileHeader.iNumberOfCells) );
		// Allocate memory to store the .amr file data
		if( !AllocateAMRQuantities() ) {
			printf( "\nFailed to AllocateAMRQuantities()\n" );
			fclose( pAMRFile );
			return( false );		
		}
		// Allocate memory to store the unique ID structure to track cell connectivity
		if( !AllocateUniqueIDStructure() ) {
			printf( "\nFailed to AllocateUniqueIDStructure()\n" );
			fclose( pAMRFile );
			return( false );
		}
		// Read the .amr file data
		for( i=0; i<AMRFileHeader.iNumberOfCells; i++ ) {
			for( j=0; j<iNumberOfColumns; j++ )
				ReadDouble( pAMRFile, &(ppfAMRQuantities[i][j]) );
			// Skip the cell ID data for the refinement algorithm
			for( j=0; j<=iMaxRL; j++ )
				fscanf( pAMRFile, "%i", &(ppiUniqueIDStructure[i][j]) );
		}
	fclose( pAMRFile );
	return( true );
}

bool CAMRFile::FreeAll( void )
{
	int i;

	if( !ppfAMRQuantities ) {
		printf( "\nNo allocated memory to be freed in CAMRFile::FreeAll() (ppfAMRQuantities)\n" );
		return( false );
	}

	if( !ppiUniqueIDStructure ) {
		printf( "\nNo allocated memory to be freed in CAMRFile::FreeAll() (ppiUniqueIDStructure)\n" );
		return( false );
	}

	for( i=0; i<AMRFileHeader.iNumberOfCells; i++ )
	{
		free( ppiUniqueIDStructure[i] );
		free( ppfAMRQuantities[i] );
	}
	free( ppiUniqueIDStructure );
	free( ppfAMRQuantities );

	return( true );
}

bool CAMRFile::DefineAMRFile( char *pszMeshDefinitionFilename )
{
	FILE *pMeshDefinitionFile;
	double *pfMeshTemplate = NULL, fScaleFactor;
	int iMinNumberOfCells, iGroupSize;
	int i;

	// Check to see whether a .amr file has already been opened or defined
	if( ppfAMRQuantities || ppiUniqueIDStructure ) {
		printf( "\nCannot overwrite protected data. Create a new instance of the AMRFile object instead\n" );
		return( false );
	}

	// Zero the time and number for a new .amr file
	AMRFileHeader.fProfileTime = 0.0;
	AMRFileHeader.iProfileNumber = 0;
	// Reset the length
	AMRFileHeader.fFullLength = 0.0;

	// STEP 1: Open and read the mesh definition file
	pMeshDefinitionFile = fopen( pszMeshDefinitionFilename, "r" );
	if( !pMeshDefinitionFile ) {
		printf ( "\nFailed to open the mesh definition file: %s\n", pszMeshDefinitionFilename );
		return( false );
	}
		// Get the scaling factor to use with the template of relative cell sizes
		ReadDouble( pMeshDefinitionFile, &fScaleFactor );
		// Get the maximum refinement level
		fscanf( pMeshDefinitionFile, "%i", &iMaxRL );
		// Get the minimum number of cells
		fscanf( pMeshDefinitionFile, "%i", &iMinNumberOfCells );
		// Allocate memory to store the template of relative cell sizes
		pfMeshTemplate = (double*)malloc( sizeof(double) * iMinNumberOfCells );
		if( !pfMeshTemplate ) {
			printf( "\nFailed to allocate memory at pfMeshTemplate\n" );
			fclose( pMeshDefinitionFile );
			return( false );
		}
		// Read the template data
		for( i=0; i<iMinNumberOfCells; i++ ) {
			ReadDouble( pMeshDefinitionFile, &(pfMeshTemplate[i]) );
			AMRFileHeader.fFullLength += fScaleFactor * pfMeshTemplate[i];
		}
	fclose( pMeshDefinitionFile );

	// Calculate the number of cells based on the minimum number all being maximally refined
	iGroupSize = pow(2,iMaxRL);
	AMRFileHeader.iNumberOfCells = iMinNumberOfCells * iGroupSize;

	// STEP 2: Allocate memory for the .amr file data and unique ID structure
	// Allocate memory to store the .amr file data
	if( !AllocateAMRQuantities() ) {
		printf( "\nFailed to AllocateAMRQuantities()\n" );
		return( false );		
	}
	// Allocate memory to store the unique ID structure to track cell connectivity
	if( !AllocateUniqueIDStructure() ) {
		printf( "\nFailed to AllocateUniqueIDStructure()\n" );
		return( false );
	}
	// Generate a unique ID structure for the new .amr file to track cell connectivity
	if( !GenerateUniqueIDStructure() ) {
		printf( "\nFailed to GenerateUniqueIDStructure()\n" );
		return( false );
	}

	// STEP 3: Calculate the cell coordinates and sizes
	double *pfNewAMRQuantities;
	double fs, fCellWidth, fEdge;
	int j, k, m;

	pfNewAMRQuantities = (double*)alloca( sizeof(double) * iNumberOfColumns );

	// Keep track of the spatial coordinate
	fs = fEdge = 0.0;
	// Keep track of the cell number to update in the new .amr file
	k = 1;
	for( i=0; i<iMinNumberOfCells; i++ ) {
		// One template cell and its maximally refined group at a time
		fCellWidth = fScaleFactor * pfMeshTemplate[i];
		fEdge += fCellWidth;
		fCellWidth /= (double)iGroupSize;
		fs += fCellWidth / 2.0;
		for( j=0; j<iGroupSize; j++ ) {
			pfNewAMRQuantities[0] = fs;
			pfNewAMRQuantities[1] = fCellWidth;
			for( m=2; m<iNumberOfColumns; m++ ) {
				// Zero the other cell quantities
				pfNewAMRQuantities[m] = 0.0;				
			}
			SetCellData( k, pfNewAMRQuantities );
			k++;
			// Calculate the new spatial coordinate
			fs += fCellWidth;
		}
		// Reset the spatial coordinate to the edge location of the largest cell which contains the group
		fs = fEdge;
	}

	// Free the memory allocated to store the template of relative cell sizes
	free( pfMeshTemplate );

	return( true );
}

bool CAMRFile::InterpolateAMRFile( CAMRFile *pAMRFile )
{
	AMRFILEHEADER AMRHeaderData;
	double *pfNewAMRQuantities, *pfAMRQuantities[2];
	double x[3], y[3];
	int i, j, k;

	// Get the header data from the interpolating .amr file
	pAMRFile->GetHeaderData( &AMRHeaderData );

	// Note: the number of columns must be the same in both .amr files
	pfNewAMRQuantities = (double*)alloca( sizeof(double) * iNumberOfColumns );
	for( i=0; i<2; i++ ) {
		pfAMRQuantities[i] = (double*)alloca( sizeof(double) * iNumberOfColumns );
	}

	j = 1;
	for( i=1; i<=AMRFileHeader.iNumberOfCells; i++ ) {
		// Get the current cell data
		GetCellData( i, pfNewAMRQuantities );
		// Find the cluster of cells in the interpolating .amr file which surround the current cell coordinate
		for( k=j; k<=AMRHeaderData.iNumberOfCells; k++ ) {
			pAMRFile->GetCellData( k, pfAMRQuantities[1] );
			if( pfAMRQuantities[1][0] > pfNewAMRQuantities[0] ) {
				if( k > 1 && k <= AMRHeaderData.iNumberOfCells ) {
					pAMRFile->GetCellData( k-1, pfAMRQuantities[0] );
					// pAMRFile->GetCellData( k, pfAMRQuantities[1] );
					if( j < k - 1 ) j = k - 1;
					break;
				} else {
					// k < 2
					pAMRFile->GetCellData( 1, pfAMRQuantities[0] );
					pAMRFile->GetCellData( 2, pfAMRQuantities[1] );
					j = 1;
					break;
				}
			}
		}
		// Trap the case where no cluster of cells in the interpolating .amr file surround the current cell coordinate
		if( k > AMRHeaderData.iNumberOfCells ) {
			pAMRFile->GetCellData( AMRHeaderData.iNumberOfCells-1, pfAMRQuantities[0] );					
			pAMRFile->GetCellData( AMRHeaderData.iNumberOfCells, pfAMRQuantities[1] );
			j = AMRHeaderData.iNumberOfCells-1;
		}

		// Set the x values for the interpolation
		x[1] = pfAMRQuantities[0][0];
		x[2] = pfAMRQuantities[1][0];

		// Interpolate the quantities from the existing .amr file into the current .amr file
		for( k=2; k<iNumberOfColumns; k++ ) {
			y[1] = pfAMRQuantities[0][k];
			y[2] = pfAMRQuantities[1][k];
			// Linear interpolation is safest
			LinearFit( x, y, pfNewAMRQuantities[0], &(pfNewAMRQuantities[k]) );
		}

		// Update the current cell with the interpolated quantities
		SetCellData( i, pfNewAMRQuantities );
	}

	return( true );
}

void CAMRFile::GetHeaderData( PAMRFILEHEADER pAMRFileHeader )
{
	memcpy( pAMRFileHeader, &AMRFileHeader, sizeof(strAMRFileHeader) );
}

bool CAMRFile::GetCellData( int iCellNumber, double *pfAMRQuantities )
{
	// Validate the requested cell is within the range of cells
	if( iCellNumber < 1 || iCellNumber > AMRFileHeader.iNumberOfCells ) {
		printf( "\nCell number %i out of range (valid range: 1 to %i )\n", iCellNumber, AMRFileHeader.iNumberOfCells );
		return( false );
	}
	// User enters a number between 1 and iNumberOfCells, but cells are indexed from 0
	iCellNumber--;
	if( !memcpy( pfAMRQuantities, &(ppfAMRQuantities[iCellNumber][0]), sizeof(double)*iNumberOfColumns ) ) {
		printf( "\nFailed to copy ppfAMRQuantities[%i] into array passed to GetCellData()\n", iCellNumber );
		return( false );
	}

	return( true );
}

bool CAMRFile::GetCellData( int iCellNumber, int iColumnNumber, double *pfAMRQuantity )
{
	// Validate the requested cell is within the range of cells
	if( iCellNumber < 1 || iCellNumber > AMRFileHeader.iNumberOfCells ) {
		printf( "\nCell number %i out of range (valid range: 1 to %i )\n", iCellNumber, AMRFileHeader.iNumberOfCells );
		return( false );
	}
	// Validate the requested column is within the range of columns
	if( iColumnNumber < 1 || iColumnNumber > iNumberOfColumns ) {
		printf( "\nColumn number %i out of range (valid range: 1 to %i )\n", iColumnNumber, iNumberOfColumns );
		return( false );
	}
	// User enters a number between 1 and iNumberOfCells, and 1 and iNumberOfColumns, but array indices are from 0
	iCellNumber--;
	iColumnNumber--;
	*pfAMRQuantity = ppfAMRQuantities[iCellNumber][iColumnNumber];

	return( true );
}

void CAMRFile::SetHeaderData( AMRFILEHEADER AMRHeader )
{
	memcpy( &AMRFileHeader, &AMRHeader, sizeof(strAMRFileHeader) );
}

bool CAMRFile::SetCellData( int iCellNumber, double *pfAMRQuantities )
{
	// Validate the requested cell is within the range of cells
	if( iCellNumber < 1 || iCellNumber > AMRFileHeader.iNumberOfCells ) {
		printf( "\nCell number %i out of range (valid range: 1 to %i )\n", iCellNumber, AMRFileHeader.iNumberOfCells );
		return( false );
	}
	// User enters a number between 1 and iNumberOfCells, but cells are indexed from 0
	iCellNumber--;
	if( !memcpy( &(ppfAMRQuantities[iCellNumber][0]), pfAMRQuantities, sizeof(double)*iNumberOfColumns ) ) {
		printf( "\nFailed to copy array passed to SetCellData() into ppfAMRQuantities[%i]\n", iCellNumber );
		return( false );
	}

	return( true );
}

bool CAMRFile::SetCellData( int iCellNumber, int iColumnNumber, double fAMRQuantity )
{
	// Validate the requested cell is within the range of cells
	if( iCellNumber < 1 || iCellNumber > AMRFileHeader.iNumberOfCells ) {
		printf( "\nCell number %i out of range (valid range: 1 to %i )\n", iCellNumber, AMRFileHeader.iNumberOfCells );
		return( false );
	}
	// Validate the requested column is within the range of columns
	if( iColumnNumber < 1 || iColumnNumber > iNumberOfColumns ) {
		printf( "\nColumn number %i out of range (valid range: 1 to %i )\n", iColumnNumber, iNumberOfColumns );
		return( false );
	}
	// User enters a number between 1 and iNumberOfCells, and 1 and iNumberOfColumns, but array indices are from 0
	iCellNumber--;
	iColumnNumber--;
	ppfAMRQuantities[iCellNumber][iColumnNumber] = fAMRQuantity;

	return( true );
}

bool CAMRFile::SaveAMRFile( char *pszAMRFilename )
{
	FILE *pAMRFile;
	int i, j;

	// Update the internal copy of the filename
	strcpy( szAMRFilename, pszAMRFilename );

	pAMRFile = fopen( szAMRFilename, "w" );
	if( !pAMRFile ) {
		printf ( "\nFailed to open .amr file: %s\n", szAMRFilename );
		return( false );
	}
		// Write the header information into the .amr file
		fprintf( pAMRFile, "%f\n%i\n%.16e\n%i\n", AMRFileHeader.fProfileTime, AMRFileHeader.iProfileNumber, AMRFileHeader.fFullLength, AMRFileHeader.iNumberOfCells );
		for( i=0; i<AMRFileHeader.iNumberOfCells; i++ ) {
			fprintf( pAMRFile, "%.16e", ppfAMRQuantities[i][0] );
			for( j=1; j<iNumberOfColumns; j++ ) {
				fprintf( pAMRFile, "\t%.16e", ppfAMRQuantities[i][j] );
			}
			for( j=0; j<(iMaxRL+1); j++ ) {
				fprintf( pAMRFile, "\t%i", ppiUniqueIDStructure[i][j] );
			}
			fprintf( pAMRFile, "\n" );
		}
	fclose( pAMRFile );
	return( true );
}