// ****
// *
// * Function bodies for the mesh tool
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/20/2021
// *
// ****

#include <stdio.h>

#include "meshtool.h"
#include "../../../source/file.h"

// **** MESH TOOL CLASS ****

// Constructor
CMeshTool::CMeshTool( char *pszConfigFilename )
{
	if( !ReadConfigFile( pszConfigFilename ) )
		return;

	// Get the existing .amr file
	pAMRFile = new CAMRFile( szAMRFilename, iAMRMaxRL );

	// Create a new .amr file
	pNewAMRFile = new CAMRFile();
	if( !pNewAMRFile->DefineAMRFile( szMeshDefinitionFilename ) )
		printf( "\nFailed to DefineAMRFile()\n" );
	// Interpolate quantities from the existing .amr file into the new .amr file
	pNewAMRFile->InterpolateAMRFile( pAMRFile );
	if( !pNewAMRFile->SaveAMRFile( szNewAMRFilename ) )
		printf( "\nFailed to SaveAMRFile()\n" );
}

// Destructor
CMeshTool::~CMeshTool( void )
{
	if( !FreeAll() )
		printf( "\nFailed to FreeAll() in ~CMeshTool()\n" );
}

bool CMeshTool::ReadConfigFile( char *pszConfigFilename )
{
	FILE *pConfigFile;

	pConfigFile = fopen( pszConfigFilename, "r" );
	if( !pConfigFile ) {
		printf ( "\nFailed to open the meshTool configuration file: %s\n", pszConfigFilename );
		return( false );
	}
		fscanf( pConfigFile, "%s", szMeshDefinitionFilename );
		fscanf( pConfigFile, "%s", szNewAMRFilename );
		fscanf( pConfigFile, "%s", szAMRFilename ); fscanf( pConfigFile, "%i", &iAMRMaxRL );
	fclose( pConfigFile );
	return( true );
}

bool CMeshTool::FreeAll( void )
{
	delete pNewAMRFile;
	delete pAMRFile;

	return( true );	
}