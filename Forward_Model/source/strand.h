// ****
// *
// * Strand class definition
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


struct SPHYData {

    // The curvi-linear coordinate along the magnetic field and the cell width
    double fs, fds;

    // The bulk velocity, number density, electron and ion temperatures
    double fv, fn, fTe, fTi;

    // The bulk velocity projected onto the line-of-sight
    double fvp;

};

typedef SPHYData PHYDATA;
typedef SPHYData* PPHYDATA;
typedef SPHYData** PPPHYDATA;


class CStrand {

    private:

    // The time stamp from the numerical experiment
    double fTimeStamp;

    // The strand diameter and length
    double fSD, fLength;

    // The number of grid cells
    int iNumCells;

    // Pointer to an array of structures containing the physical properties of the strand as a function of position in field-aligned coordinates
    PPHYDATA pPHYData;

    // The number of non-equilibrium (NEQ) elements and their atomic numbers
    int iNumNEQElements, *piNEQ_Z;

    // The non-equilibrium population fractions: pppNEQ[X][Y][Z], where X = cell, Y = element and Z = ion
    double ***pppfNEQ;

    // Function to initialise the Strand object
    void Initialise( char *pszWorkingDirectory, int iProfileNumber, double fH );

    // Function to free all allocated memory
    void FreeAll( void );

    public:

    // Constructor
    CStrand( char *pszWorkingDirectory, int iProfileNumber, double fH );
        
    // Destructor
    ~CStrand( void );

    // Function to return the number of grid cells
    int GetNumCells( void );

    // Function to return the timestamp
    double GetTimeStamp( void );

    // Function to return the diameter of the strand
    double GetDiameter( void );

    // Function to return the length of the strand
    double GetLength( void );

    // Function to return the physical properties of a specified grid cell
    void GetPHYData( int iCell, PPHYDATA pStrandPHYData );

    // Function to return the number of elements for which the non-equilibrium ionisation state exists
    int GetNumNEQElements( void );

    // Function to return the atomic number of the specified element in the list for a specified strand
    int GetAtomicNumber( int iElement );

    // Function to return the non-equilibrium ion population fraction of the specified ion / element for a specified cell
    double GetNEQ( int iCell, int iElement, int iIon );

};

typedef CStrand* PSTRAND;
typedef CStrand** PPSTRAND;
