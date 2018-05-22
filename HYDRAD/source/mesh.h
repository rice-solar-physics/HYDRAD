// ****
// *
// * Class definition of the adaptive mesh
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 05/22/2018
// *
// ****

#include <time.h>

#include "eqns.h"


// **** ADAPTIVE MESH CLASS ****

// Define the adaptive mesh class
class CAdaptiveMesh : private CEquations {

    private:
	
	// The run (wall) and mesh (model) time elapsed, and the current time-step
    double run_time, mesh_time, mesh_delta_t;

    // The file number for the next set of output profiles
    int iFileNumber;

    void CreateInitialMesh( void );

    void ZeroCellProperties( PCELLPROPERTIES pCellProperties );

    void Adapt( void );

    void FreeCurrentRow( void );
    void FreePreviousRow( void );

    void Solve( void );
    void Integrate( void );

	void ShowProgress( clock_t *ptimer );
    void WriteToFile( void );

    public:

    // Constructor
    CAdaptiveMesh( void );
	
    // Destructor
    ~CAdaptiveMesh( void );
	
};

// Define a type for the adaptive mesh class
typedef CAdaptiveMesh* PMESH;
