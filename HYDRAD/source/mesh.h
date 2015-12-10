// ****
// *
// * Class definition of the adaptive mesh
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 05/11/2012
// *
// ****


#include "eqns.h"


// **** ADAPTIVE MESH CLASS ****

// Define the adaptive mesh class
class CAdaptiveMesh : private CEquations {

    private:
	
    // The mesh time elapsed and the time-step taken
    double mesh_time, mesh_delta_t;

    // The file number for the next set of output profiles
    int iFileNumber;

    void CreateInitialMesh( void );

    void ZeroCellProperties( PCELLPROPERTIES pCellProperties );

    void Adapt( void );

    void FreeCurrentRow( void );
    void FreePreviousRow( void );

    void Solve( void );
    void Integrate( void );

    void WriteToFile( void );

    public:

    // Constructor
    CAdaptiveMesh( void );
	
    // Destructor
    ~CAdaptiveMesh( void );
	
};

// Define a type for the adaptive mesh class
typedef CAdaptiveMesh* PMESH;
