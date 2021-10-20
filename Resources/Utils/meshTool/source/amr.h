// ****
// *
// * Class definition of the HYDRAD .amr file reader
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/20/2021
// *
// ****

// Need to know whether electron mass density is included in the .amr file
#include "../../../../Radiation_Model/source/config.h"

// **** AMR FILE HEADER STRUCTURE ****

// Define the HYDRAD .amr file header structure
struct strAMRFileHeader {

	double fProfileTime, fFullLength;
	int iProfileNumber, iNumberOfCells;

};

// Define a type for the HYDRAD .amr file header structure
typedef strAMRFileHeader AMRFILEHEADER;
typedef strAMRFileHeader* PAMRFILEHEADER;

// **** AMR FILE CLASS ****

// Define the HYDRAD .amr file class
class CAMRFile {

	private:

	char szAMRFilename[256];
	int iMaxRL;
#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE 
	int iNumberOfColumns = 7;
#else // NLTE_CHROMOSPHERE
	int iNumberOfColumns = 6;
#endif // NLTE_CHROMOSPHERE
#else // OPTICALLY_THICK_RADIATION
	int iNumberOfColumns = 6;
#endif // OPTICALLY_THICK_RADIATION
	AMRFILEHEADER AMRFileHeader;
	double **ppfAMRQuantities = NULL;
	int **ppiUniqueIDStructure = NULL;

	// Memory allocation functions
	bool AllocateAMRQuantities( void );
	bool AllocateUniqueIDStructure( void );
	// Function to generate the unique ID structure for cell connectivity
	bool GenerateUniqueIDStructure( void );

	// Open and read an existing .amr file
	bool ReadAMRFile( char *pszAMRFilename );

	// Function to free all allocated memory
	bool FreeAll( void );

	public:

	// Constructor
	CAMRFile( void );
	CAMRFile( char *pszAMRFilename, int iAMRMaxRL );

	// Destructor
	~CAMRFile( void );

	// Create a new .amr file using the specified mesh definition file
	bool DefineAMRFile( char *pszMeshDefinitionFilename );

	// Function to interpolate quantities from an existing .amr file into the current .amr file 
	bool InterpolateAMRFile( CAMRFile *pAMRFile );

	// Functions for retrieving .amr file information
	void GetHeaderData( PAMRFILEHEADER pAMRFileHeader );
	bool GetCellData( int iCellNumber, double *pfAMRQuantities );
	bool GetCellData( int iCellNumber, int iColumnNumber, double *pfAMRQuantity );

	// Functions for setting .amr file information
	void SetHeaderData( AMRFILEHEADER AMRFileHeader );
	bool SetCellData( int iCellNumber, double *pfAMRQuantities );
	bool SetCellData( int iCellNumber, int iColumnNumber, double fAMRQuantity );

	// File operation functions
	bool SaveAMRFile( char *pszAMRFilename );
};

// Define a type for the HYDRAD .amr file class
typedef CAMRFile* PAMRFILE;