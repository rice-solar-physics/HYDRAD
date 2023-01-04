// ****
// *
// * A utility to generate a new mesh onto which an existing mesh is interpolated
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>

#include "meshtool.h"

int main( int argc, char **argv )
{
	PMESHTOOL pMeshTool;

	if( argc < 2 ) {
		printf( "\nUsage: meshTool \"config filename\"\n" );
		printf( "\nUsing default config file: \"config.cfg\"\n" );
		pMeshTool = new CMeshTool( (char*)"config.cfg" );
	} else {
		printf( "\nUsing config file: \"%s\"\n", argv[1] );
		pMeshTool = new CMeshTool( argv[1] );
	}

	delete pMeshTool;

	return 0;
}