// ****
// *
// * Class definition of the mesh tool
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/13/2021
// *
// ****

#include "amr.h"

// **** MESH TOOL CLASS ****

// Define the mesh tool class
class CMeshTool {

    private:
	
	char szMeshDefinitionFilename[256], szNewAMRFilename[256], szAMRFilename[256];
	PAMRFILE pNewAMRFile, pAMRFile;

	// Open and read the configuration file
	bool ReadConfigFile( char *pszConfigFilename );
	
	// Function to free all allocated memory
	bool FreeAll( void );

	public:
	
    // Constructor
	CMeshTool( char *pszConfigFilename );
	
    // Destructor
    ~CMeshTool( void );
};

// Define a type for the mesh tool class
typedef CMeshTool* PMESHTOOL;