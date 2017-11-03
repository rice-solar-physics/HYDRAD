// ****
// *
// * Loop class definition
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


#include "strand.h"


class CLoop {

    private:

    // The directory from which to read the data for the strands comprising the loop
    char szWorkingDirectory[256];

    // The number of strands comprising the loop and the range of profiles to use as individual strands
    int iNumStrands, iStrandRange[3];

    // Array of pointers to individual strands
    PPSTRAND ppStrand;

    // Function to initialise the Loop object
    void Initialise( char *pszWorkingDirectory, int iFrom, int iTo, int iStep, double fdH );

    // Function to free all allocated memory
    void FreeAll( void );

    public:

    // Constructor
    CLoop( char *pszWorkingDirectory, int iFrom, int iTo, int iStep, double fdH );
        
    // Destructor
    ~CLoop( void );

    // Function to the return a string pointing to the working directory name
    char *pGetWorkingDirectory( void );

    // Functions to return the number of strands and the range of profiles used as individual strands
    int GetNumStrands( void );
    void GetStrandRange( int *piStrandRange );

    // Function to return the timestamp for a specified strand
    double GetTimeStamp( int iStrand );

    // Function to return the number of grid cells for a specified strand
    int GetNumCells( int iStrand );

    // Function to return the diameter of a specified strand
    double GetDiameter( int iStrand );

    // Function to return the length of a specified strand
    double GetLength( int iStrand );

    // Function to return the physical properties of a specified grid cell in a specified strand
    void GetPHYData( int iStrand, int iCell, PPHYDATA pStrandPHYData );

    // Function to return the number of elements for which the non-equilibrium ionisation state exists for a specified strand
    int GetNumNEQElements( int iStrand );

    // Function to return the atomic number of the specified element in the list for a specified strand
    int GetAtomicNumber( int iStrand, int iElement );

    // Function to return the non-equilibrium ion population fraction of the specified ion / element for a specified cell / strand
    double GetNEQ( int iStrand, int iCell, int iElement, int iIon );

};

typedef CLoop* PLOOP;
