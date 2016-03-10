// ****
// *
// * Class definition of the time-dependent hydrodynamic
// * equations, inherited by the adaptive mesh class
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/03/2016
// *
// ****


#include "cell.h"


// **** EQUATIONS CLASS ****

// Define the equations class
class CEquations {

    private:

    double **ppGravity;
    int igdp;

#ifdef USE_TABULATED_CROSS_SECTION
    double **ppCrossSection;
    int icsdp;
#endif // USE_TABULATED_CROSS_SECTION

    double lower_radiation_temperature_boundary;

    void Initialise( void );
    void FreeAll( void );

    // Function for finding the gravitational acceleration from the look-up table
    double CalculateGravity( double s );

#ifdef USE_TABULATED_CROSS_SECTION
    // Function for finding the cross-section from the look-up table
    double CalculateCrossSection( double s );
#endif // USE_TABULATED_CROSS_SECTION

    // Function for finding the smallest time-scale
    void GetSmallestTimeScale( double *delta_t, int iFirstStep );

#ifdef USE_KINETIC_MODEL
    int iNumCells;

    // Functions for the Spitzer-Harm part of the solution
    // Tabulated values from tables I and II (for Z = 1) in Spitzer & Harm, 1953, Phys. Rev., 89, 977
    double SH_Table[51][3];
    void Get_SH_Table( void );
#endif // USE_KINETIC_MODEL

    public:

    // User specifiable and loop parameters
    PARAMETERS Params;

    // Pointer to the heating model
    PHEAT pHeat;

    // Pointers to the radiation models
    PRADIATION pRadiation, pRadiation2;

#ifdef OPTICALLY_THICK_RADIATION
    // Pointers to the ions for which optically-thick radiative emission
    // will be calculated
    POPTICALLYTHICKION pHI, pMgII, pCaII;

    // Pointer to the centre of the row at the current time (approx. the loop apex)
    PCELL pCentreOfCurrentRow;
#endif // OPTICALLY_THICK_RADIATION

    // Pointer to the left-most cell at the previous time (the start of the previous row)
    PCELL pStartOfPreviousRow;

    // Pointer to the left-most and right-most cells at the current time (the start and end of the current row)
    PCELL pStartOfCurrentRow, pEndOfCurrentRow;

    // Pointer to the active cell
    PCELL pActiveCell;

    // Constructor
    CEquations( void );

    // Destructor
    ~CEquations( void );

    // Function for calculating physical quantities
    void CalculatePhysicalQuantities( void );

    // Function for evaluating the terms of the equations
    void EvaluateTerms( double current_time, double *delta_t, int iFirstStep );

    // Functions to provide a 2nd order accurate numerical integration of
    // the system of equations
    void Half_Time_Step( PCELLPROPERTIES pNewCellProperties, double delta_t );
    void Full_Time_Step( PCELLPROPERTIES CellProperties, double delta_t );

#ifdef USE_KINETIC_MODEL
    // Pointer to an indexed list of cells
    PCELL *ppCellList;

    void CountCells( void );
    void CreateIndexedCellList( void );

    void CalculateKineticModel( int iFirstStep );
    void CalculateNonMaxDFN( void );
#endif // USE_KINETIC_MODEL

};
