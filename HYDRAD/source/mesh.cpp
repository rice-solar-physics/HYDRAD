// ****
// *
// * Function bodies for the class definition of the adaptive mesh
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/11/2018
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mesh.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"

#ifdef OPENMP
	#include <omp.h>
#endif // OPENMP


// **** ADAPTIVE MESH CLASS ****

// Constructor
CAdaptiveMesh::CAdaptiveMesh( void )
{
// Create the initial mesh given the steady-state profiles
CreateInitialMesh();

// Solve the equations
Solve();
}

// Destructor
CAdaptiveMesh::~CAdaptiveMesh( void )
{
// Free all of the memory allocated within the mesh object
FreeCurrentRow();
}

void CAdaptiveMesh::ZeroCellProperties( PCELLPROPERTIES pCellProperties )
{
memset( pCellProperties, 0, sizeof(CELLPROPERTIES) );
}

void CAdaptiveMesh::CreateInitialMesh( void )
{
FILE *pFile;
PCELL pPreviousCell = NULL;
CELLPROPERTIES CellProperties;
int i, j;

printf( "\nProcessing the initial conditions...\n" );

#ifdef NON_EQUILIBRIUM_RADIATION
PCELL pNextActiveCell;
FILE *pIonFile = NULL;
char szIonFilename[32];
double fTemp;
#endif // NON_EQUILIBRIUM_RADIATION

// ******************************************************************************
// *                                                                            *
// *    CREATE THE COMPUTATIONAL MESH                                           *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    INITIALISE THE MESH                                                     *
// ******************************************************************************

// Reset the initial cell properties to zero before the initial mesh is created
ZeroCellProperties( &CellProperties );

// Open the initial density, momentum density and thermal + kinetic energy density profiles
pFile = fopen( Params.Profiles, "r" );

ReadDouble( pFile, &mesh_time );
fscanf( pFile, "%i", &iFileNumber );

ReadDouble( pFile, &Params.L );

fscanf( pFile, "%i", &Params.iNumberOfCells );

// Create the initial mesh using values from the user-specified .amr file
for( i = 0; i < Params.iNumberOfCells; i++ )
{
    // Fill in the cell properties structure

    ReadDouble( pFile, &(CellProperties.s[1]) );
    ReadDouble( pFile, &(CellProperties.cell_width) );

    CellProperties.s[0] = CellProperties.s[1] - ( 0.5 * CellProperties.cell_width );
    CellProperties.s[2] = CellProperties.s[1] + ( 0.5 * CellProperties.cell_width );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
    ReadDouble( pFile, &(CellProperties.rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
    ReadDouble( pFile, &(CellProperties.rho[1]) );
    ReadDouble( pFile, &(CellProperties.rho_v[1]) );
    ReadDouble( pFile, &(CellProperties.TE_KE[1][ELECTRON]) );
    ReadDouble( pFile, &(CellProperties.TE_KE[1][HYDROGEN]) );

    // The number of refinement levels in the .amr file MUST match MAX_REFINEMENT_LEVEL
    fscanf( pFile, "%i", &(CellProperties.iRefinementLevel) );
    for( j=1; j<=MAX_REFINEMENT_LEVEL; j++ )
	fscanf( pFile, "%i", &(CellProperties.iUniqueID[j]) );

#ifdef NON_EQUILIBRIUM_RADIATION
    // Reset the pointer to the ion population fractions object
    CellProperties.pIonFrac = NULL;
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
    CellProperties.pKinetic = new CKinetic();
#endif // USE_KINETIC_MODEL

    // Create a new cell
    pActiveCell = new CAdaptiveMeshCell( &CellProperties );

    // Set the RIGHT pointer to NULL
    pActiveCell->SetPointer( RIGHT, NULL );

    // If this is the first cell
    if( !i )
    {
        // Need to know the address of the left-most cell at each time t
	// in order to derive the row of cells at t + delta_t
	pStartOfCurrentRow = pActiveCell;
	pStartOfPreviousRow = NULL;		

	// Set the LEFT pointer of the left-most cell to NULL
	pActiveCell->SetPointer( LEFT, NULL );
    }
    // If this is not the left-most cell set the RIGHT pointer of the previous cell
    // to point at the new cell and the LEFT pointer of the new cell to point at
    // the previous cell
    else
    {
	pPreviousCell->SetPointer( RIGHT, pActiveCell );
	pActiveCell->SetPointer( LEFT, pPreviousCell );
    }
	
    // Keep track of the cell pointer so that the appropriate LEFT and RIGHT pointers can be set
    // next time around
    pPreviousCell = pActiveCell;
}

fclose( pFile );

// ******************************************************************************
// *    ADAPT THE INITIAL MESH                                                  *
// ******************************************************************************

#ifdef ADAPT
if( mesh_time == 0.0 )
    Adapt();
#endif // ADAPT

#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
	CreateIndexedCellList();
#endif // OPENMP || USE_KINETIC_MODEL

#if defined (OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE)
InitialiseRadiativeRates();
CalculateInitialPhysicalQuantities();
#else // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
CalculatePhysicalQuantities();
#endif // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE

// ******************************************************************************
// *    INITIALISE THE IONISATION STATE                                         *
// ******************************************************************************

#ifdef NON_EQUILIBRIUM_RADIATION
// If the mesh time is greater than 0.0 then the nonequilibrium ion population corresponding to
// the current time should be used
if( mesh_time > 0.0 )
{
    sprintf( szIonFilename, "Results/profile%i.ine", iFileNumber );
    pIonFile = fopen( szIonFilename, "r" );
}

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    // Initialise the fractional populations of the ions
    CellProperties.pIonFrac = new CIonFrac( NULL, (char *)"Radiation_Model/config/elements_neq.cfg", pRadiation );

    // If the mesh time is greater than 0.0 then get the nonequilibrium ion population corresponding to
    // the current time and position
    if( mesh_time > 0.0 )
    {
        ReadDouble( pIonFile, &fTemp );
        CellProperties.pIonFrac->ReadAllIonFracFromFile( pIonFile );
    }
    else
    {
	// Make sure the ionisation balance EXACTLY matches equilibrium at t = 0 seconds
#ifdef DENSITY_DEPENDENT_RATES
	CellProperties.pIonFrac->ResetAllIonFrac( log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]) );
#else // DENSITY_DEPENDENT_RATES
	CellProperties.pIonFrac->ResetAllIonFrac( log10(CellProperties.T[ELECTRON]) );
#endif // DENSITY_DEPENDENT_RATES
    }

    pActiveCell->UpdateCellProperties( &CellProperties );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

if( mesh_time > 0.0 )
    fclose( pIonFile );
#endif // NON_EQUILIBRIUM_RADIATION

// ******************************************************************************
// *    OUTPUT THE INITIAL STATE                                                *
// ******************************************************************************

EvaluateTerms( mesh_time, &mesh_delta_t, TRUE );
WriteToFile();
}

#ifdef ADAPT
void CAdaptiveMesh::Adapt( void )
{
PCELL pNextActiveCell, pFarLeftCell, pLeftCell, pRightCell, pFarRightCell, pNewCell[2];
CELLPROPERTIES FarLeftCellProperties, LeftCellProperties, CellProperties, RightCellProperties, FarRightCellProperties, NewCellProperties[2];
double drho = 0.0, dTE_KEe = 0.0, dTE_KEh = 0.0, drho_e = 0.0, x[6], y[6];
int iProlonged, iRestricted, j;

#ifdef NON_EQUILIBRIUM_RADIATION
double **ppIonFrac[6];
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef LINEAR_RESTRICTION
#else // LINEAR_RESTRICTION
double error;
#endif // LINEAR_RESTRICTION

#ifdef ENFORCE_CONSERVATION
double fWeight, temp1, temp2, temp3, temp4;
#endif// ENFORCE_CONSERVATION

// ******************************************************************************
// *                                                                            *
// *    PROLONGATION                                                            *
// *                                                                            *
// ******************************************************************************

// Look for cells that can be merged

do {

    iProlonged = FALSE;

    pNextActiveCell = pStartOfCurrentRow;

    while( pNextActiveCell->pGetPointer( RIGHT ) )
    {
        pActiveCell = pNextActiveCell;
	pActiveCell->GetCellProperties( &CellProperties );

	pRightCell = pActiveCell->pGetPointer( RIGHT );
	pRightCell->GetCellProperties( &RightCellProperties );

	// Only neighbouring cells derived from the same parent cell can be merged
	if( CellProperties.iRefinementLevel > 0 &&
            CellProperties.iRefinementLevel == RightCellProperties.iRefinementLevel &&
            CellProperties.iUniqueID[CellProperties.iRefinementLevel] == RightCellProperties.iUniqueID[RightCellProperties.iRefinementLevel] )
	{
#ifdef REFINE_ON_DENSITY
            drho = 1.0 - ( min( CellProperties.rho[1], RightCellProperties.rho[1] ) / max( CellProperties.rho[1], RightCellProperties.rho[1] ) );
#endif // REFINE_ON_DENSITY
#ifdef REFINE_ON_ELECTRON_ENERGY
            dTE_KEe = 1.0 - ( min( CellProperties.TE_KE[1][ELECTRON], RightCellProperties.TE_KE[1][ELECTRON] ) / max( CellProperties.TE_KE[1][ELECTRON], RightCellProperties.TE_KE[1][ELECTRON] ) );
#endif // REFINE_ON_ELECTRON_ENERGY
#ifdef REFINE_ON_HYDROGEN_ENERGY
            dTE_KEh = 1.0 - ( min( CellProperties.TE_KE[1][HYDROGEN], RightCellProperties.TE_KE[1][HYDROGEN] ) / max( CellProperties.TE_KE[1][HYDROGEN], RightCellProperties.TE_KE[1][HYDROGEN] ) );
#endif // REFINE_ON_HYDROGEN_ENERGY
#if defined (OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE)
#ifdef REFINE_ON_ELECTRON_DENSITY
            drho_e = 1.0 - ( min( CellProperties.rho_e, RightCellProperties.rho_e ) / max( CellProperties.rho_e, RightCellProperties.rho_e ) );
#endif // REFINE_ON_ELECTRON_DENSITY
#endif // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
            if( drho < MIN_FRAC_DIFF && dTE_KEe < MIN_FRAC_DIFF && dTE_KEh < MIN_FRAC_DIFF && drho_e < MIN_FRAC_DIFF )
            {
                iProlonged = TRUE;
				Params.iNumberOfCells--;
				
		ZeroCellProperties( &(NewCellProperties[0]) );

#ifdef NON_EQUILIBRIUM_RADIATION
		// Create a new ionfrac object for the new cell
                if( CellProperties.pIonFrac )
                    NewCellProperties[0].pIonFrac = new CIonFrac( CellProperties.pIonFrac, 0, pRadiation );
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
		// Create a new kinetic object for the new cell
		NewCellProperties[0].pKinetic = new CKinetic();
#endif // USE_KINETIC_MODEL

		NewCellProperties[0].iRefinementLevel = CellProperties.iRefinementLevel - 1;

		for( j=1; j<=NewCellProperties[0].iRefinementLevel; j++ )
                    NewCellProperties[0].iUniqueID[j] = CellProperties.iUniqueID[j];

		NewCellProperties[0].cell_width = CellProperties.cell_width + RightCellProperties.cell_width;

		NewCellProperties[0].s[0] = CellProperties.s[0];
		NewCellProperties[0].s[1] = CellProperties.s[2];
		NewCellProperties[0].s[2] = RightCellProperties.s[2];

// ******************************************************************************
// *    ENSURE MASS, MOMENTUM AND ENERGY IS CONSERVED                           *
// ******************************************************************************

		// Averaging the mass density, momentum density and energy density of the two grid cells to calculate the corresponding quantities for the merged cell guarantees conservation

		NewCellProperties[0].rho[1] = ( CellProperties.rho[1] + RightCellProperties.rho[1] ) / 2.0;
		
		pLeftCell = pActiveCell->pGetPointer( LEFT );
		pFarRightCell = pRightCell->pGetPointer( RIGHT );
		if( pLeftCell && pFarRightCell )
		{
		    NewCellProperties[0].rho_v[1] = ( CellProperties.rho_v[1] + RightCellProperties.rho_v[1] ) / 2.0;
		}
		else
		{
// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************
		    NewCellProperties[0].rho_v[1] = 0.0;
		}

		for( j=0; j<SPECIES; j++ )
                    NewCellProperties[0].TE_KE[1][j] = ( CellProperties.TE_KE[1][j] + RightCellProperties.TE_KE[1][j] ) / 2.0;

#ifdef NON_EQUILIBRIUM_RADIATION
                if( NewCellProperties[0].pIonFrac )
                {
                    x[1] = CellProperties.s[1];
                    x[2] = RightCellProperties.s[1];
                    ppIonFrac[1] = CellProperties.pIonFrac->ppGetIonFrac();
                    ppIonFrac[2] = RightCellProperties.pIonFrac->ppGetIonFrac();
                    NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[0].s[1] );
                }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
		NewCellProperties[0].rho_e = ( CellProperties.rho_e + RightCellProperties.rho_e ) / 2.0;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

		pNewCell[0] = new CAdaptiveMeshCell( &(NewCellProperties[0]) );

		if( pLeftCell )
		{
                    pLeftCell->SetPointer( RIGHT, pNewCell[0] );
                    pNewCell[0]->SetPointer( LEFT, pLeftCell );
		}
		else
		{
                    pNewCell[0]->SetPointer( LEFT, NULL );
                    pStartOfCurrentRow = pNewCell[0];
		}

		if( pFarRightCell )
		{
                    pNewCell[0]->SetPointer( RIGHT, pFarRightCell );
                    pFarRightCell->SetPointer( LEFT, pNewCell[0] );
                    pNextActiveCell = pFarRightCell;
		}
		else
		{
                    pNewCell[0]->SetPointer( RIGHT, NULL );
                    pNextActiveCell = pNewCell[0];
		}

#ifdef NON_EQUILIBRIUM_RADIATION
                if( CellProperties.pIonFrac && RightCellProperties.pIonFrac )
                {
                    delete CellProperties.pIonFrac;
                    delete RightCellProperties.pIonFrac;
                }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
		delete CellProperties.pKinetic;
		delete RightCellProperties.pKinetic;
#endif // USE_KINETIC_MODEL

		delete pActiveCell;
		delete pRightCell;
            }
            else
            {
                pNextActiveCell = pRightCell;
            }
	}
	else
	{
            pNextActiveCell = pRightCell;
	}
    }
} while ( iProlonged );

// ******************************************************************************
// *                                                                            *
// *    RESTRICTION                                                             *
// *                                                                            *
// ******************************************************************************

// Look for cells that need to be refined

do {

    iRestricted = FALSE;

    pNextActiveCell = pStartOfCurrentRow;

    while( pNextActiveCell->pGetPointer( RIGHT ) )
    {
        pActiveCell = pNextActiveCell;
	pActiveCell->GetCellProperties( &CellProperties );

	pRightCell = pActiveCell->pGetPointer( RIGHT );
	pRightCell->GetCellProperties( &RightCellProperties );
#ifdef REFINE_ON_DENSITY
	drho = 1.0 - ( min( CellProperties.rho[1], RightCellProperties.rho[1] ) / max( CellProperties.rho[1], RightCellProperties.rho[1] ) );
#endif // REFINE_ON_DENSITY
#ifdef REFINE_ON_ELECTRON_ENERGY
	dTE_KEe = 1.0 - ( min( CellProperties.TE_KE[1][ELECTRON], RightCellProperties.TE_KE[1][ELECTRON] ) / max( CellProperties.TE_KE[1][ELECTRON], RightCellProperties.TE_KE[1][ELECTRON] ) );
#endif // REFINE_ON_ELECTRON_ENERGY
#ifdef REFINE_ON_HYDROGEN_ENERGY
	dTE_KEh = 1.0 - ( min( CellProperties.TE_KE[1][HYDROGEN], RightCellProperties.TE_KE[1][HYDROGEN] ) / max( CellProperties.TE_KE[1][HYDROGEN], RightCellProperties.TE_KE[1][HYDROGEN] ) );
#endif // REFINE_ON_HYDROGEN_ENERGY
#if defined (OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE)
#ifdef REFINE_ON_ELECTRON_DENSITY
            drho_e = 1.0 - ( min( CellProperties.rho_e, RightCellProperties.rho_e ) / max( CellProperties.rho_e, RightCellProperties.rho_e ) );
#endif // REFINE_ON_ELECTRON_DENSITY
#endif // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
	if( ( drho > MAX_FRAC_DIFF || dTE_KEe > MAX_FRAC_DIFF || dTE_KEh > MAX_FRAC_DIFF || drho_e > MAX_FRAC_DIFF || abs( CellProperties.iRefinementLevel - RightCellProperties.iRefinementLevel ) > 1 ) && ( CellProperties.iRefinementLevel < MAX_REFINEMENT_LEVEL || RightCellProperties.iRefinementLevel < MAX_REFINEMENT_LEVEL ) )
	{
            iRestricted = TRUE;
			Params.iNumberOfCells++;
			
            if( CellProperties.iRefinementLevel <= RightCellProperties.iRefinementLevel )
            {
                // Refine the current cell

                ZeroCellProperties( &(NewCellProperties[0]) );
		ZeroCellProperties( &(NewCellProperties[1]) );

#ifdef NON_EQUILIBRIUM_RADIATION
                if( CellProperties.pIonFrac )
                {
                    // Create a new ionfrac object for each new cell
                    NewCellProperties[0].pIonFrac = new CIonFrac( CellProperties.pIonFrac, 0, pRadiation );
                    NewCellProperties[1].pIonFrac = new CIonFrac( CellProperties.pIonFrac, 0, pRadiation );
                }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
		// Create a new kinetic object for each new cell
		NewCellProperties[0].pKinetic = new CKinetic();
		NewCellProperties[1].pKinetic = new CKinetic();
#endif // USE_KINETIC_MODEL

		NewCellProperties[0].iRefinementLevel = CellProperties.iRefinementLevel + 1;
		NewCellProperties[1].iRefinementLevel = NewCellProperties[0].iRefinementLevel;

		for( j=1; j<=CellProperties.iRefinementLevel; j++ )
		{
                    NewCellProperties[0].iUniqueID[j] = CellProperties.iUniqueID[j];
                    NewCellProperties[1].iUniqueID[j] = CellProperties.iUniqueID[j];
		}

		NewCellProperties[0].cell_width = 0.5 * CellProperties.cell_width;
		NewCellProperties[1].cell_width = NewCellProperties[0].cell_width;

		NewCellProperties[0].s[0] = CellProperties.s[0];
		NewCellProperties[0].s[1] = NewCellProperties[0].s[0] + ( 0.5 * NewCellProperties[0].cell_width );
		NewCellProperties[0].s[2] = CellProperties.s[1];
		NewCellProperties[1].s[0] = CellProperties.s[1];
		NewCellProperties[1].s[1] = NewCellProperties[1].s[0] + ( 0.5 * NewCellProperties[0].cell_width );
		NewCellProperties[1].s[2] = CellProperties.s[2];
		
		pLeftCell = pActiveCell->pGetPointer( LEFT );

		if( pLeftCell )
		{
                    pLeftCell->GetCellProperties( &LeftCellProperties );

                    do {

                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] = rand();

                    } while ( NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == LeftCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] ||
                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == RightCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] );

                    NewCellProperties[1].iUniqueID[NewCellProperties[1].iRefinementLevel] = NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel];

                    pFarLeftCell = pLeftCell->pGetPointer( LEFT );
                    pFarRightCell = pRightCell->pGetPointer( RIGHT );
			
                    if( pFarLeftCell && pFarRightCell )
                    {
			pFarLeftCell->GetCellProperties( &FarLeftCellProperties );
			pFarRightCell->GetCellProperties( &FarRightCellProperties );

			// NEW CELLS

			x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];
			x[5] = FarRightCellProperties.s[1];

			y[1] = FarLeftCellProperties.rho[1];
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
			y[5] = FarRightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]), &error );
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[1] = FarLeftCellProperties.rho_v[1];
			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
			y[5] = FarRightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]) );
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]), &error );
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[1] = FarLeftCellProperties.TE_KE[1][j];
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
                            y[5] = FarRightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                            LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
#else
                            FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]), &error );
                            FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac && NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[1] = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[5] = FarRightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 2, NewCellProperties[0].s[1] );
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[2]), &(ppIonFrac[2]), 2, NewCellProperties[1].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 4, NewCellProperties[0].s[1] );
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 4, NewCellProperties[1].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION
                    

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = FarLeftCellProperties.rho_e;
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
			y[5] = FarRightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
                        LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
                        LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#else
                        FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e), &error );
                        FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
		    }
                    else if( !pFarLeftCell )
                    {
// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// LEFT-HAND BOUNDARY
                        
                        // Implement hydrostatic boundary conditions in the left-most of the two new cells

                        x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];

                        y[1] = LeftCellProperties.rho[1];
                        y[2] = CellProperties.rho[1];
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );

                        NewCellProperties[0].rho_v[1] = 0.0;

                        for( j=0; j<SPECIES; j++ )
                        {
                            y[1] = LeftCellProperties.TE_KE[1][j];
                            y[2] = CellProperties.TE_KE[1][j];
                            LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                        }

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac )
                        {
                            ppIonFrac[1] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = CellProperties.pIonFrac->ppGetIonFrac();
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[0].s[1] );
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = LeftCellProperties.rho_e;
                        y[2] = CellProperties.rho_e;
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

			pFarRightCell->GetCellProperties( &FarRightCellProperties );

			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];
			x[5] = FarRightCellProperties.s[1];
						
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
			y[5] = FarRightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
			y[5] = FarRightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
                            y[5] = FarRightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
#else
                            FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[5] = FarRightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[2]), &(ppIonFrac[2]), 2, NewCellProperties[1].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 4, NewCellProperties[1].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
			y[5] = FarRightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
                    }
                    else if( !pFarRightCell )
                    {
			pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

			x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];

			y[1] = FarLeftCellProperties.rho[1];
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[1] = FarLeftCellProperties.rho_v[1];
			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[1] = FarLeftCellProperties.TE_KE[1][j];
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
#else
                            FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac )
                        {
                            ppIonFrac[1] = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 2, NewCellProperties[0].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 4, NewCellProperties[0].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = FarLeftCellProperties.rho_e;
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// RIGHT-HAND BOUNDARY

                        // Implement hydrostatic boundary conditions in the right-most of the two new cells

                        x[1] = CellProperties.s[1];
			x[2] = RightCellProperties.s[1];

                        y[1] = CellProperties.rho[1];
                        y[2] = RightCellProperties.rho[1];
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );

                        NewCellProperties[1].rho_v[1] = 0.0;

                        for( j=0; j<SPECIES; j++ )
                        {
                            y[1] = CellProperties.TE_KE[1][j];
                            y[2] = RightCellProperties.TE_KE[1][j];
                            LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
                        }

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[1] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[1].s[1] );
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = CellProperties.rho_e;
                        y[2] = RightCellProperties.rho_e;
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
                    }
		}
		else
		{
                    do {

                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] = rand();

                    } while ( NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == RightCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] );

                    NewCellProperties[1].iUniqueID[NewCellProperties[1].iRefinementLevel] = NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel];

// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// LEFT-HAND BOUNDARY

                    // NEW CELLS

                    x[1] = CellProperties.s[1];
		    x[2] = RightCellProperties.s[1];

                    y[1] = CellProperties.rho[1];
                    y[2] = RightCellProperties.rho[1];
                    LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
                    LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );

                    NewCellProperties[0].rho_v[1] = 0.0;
                    NewCellProperties[1].rho_v[1] = 0.0;

                    for( j=0; j<SPECIES; j++ )
                    {
                        y[1] = CellProperties.TE_KE[1][j];
                        y[2] = RightCellProperties.TE_KE[1][j];
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
                    }

#ifdef NON_EQUILIBRIUM_RADIATION
                    if( NewCellProperties[0].pIonFrac && NewCellProperties[1].pIonFrac )
                    {
                        ppIonFrac[1] = CellProperties.pIonFrac->ppGetIonFrac();
                        ppIonFrac[2] = RightCellProperties.pIonFrac->ppGetIonFrac();
                        NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[0].s[1] );
                        NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[1].s[1] );
                    }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
		    y[1] = CellProperties.rho_e;
                    y[2] = RightCellProperties.rho_e;
                    LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
                    LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
		}

// ******************************************************************************
// *    ENSURE MASS AND ENERGY IS CONSERVED BY THE RESTRICTION OPERATOR         *
// ******************************************************************************
#ifdef ENFORCE_CONSERVATION
		// A correction is applied to the mass and energy densities in the event that the sum of the integrated quantities in the new cells does not equal the integrated quantity in the original cell

		fWeight = ( 2.0 * CellProperties.rho[1] ) / ( NewCellProperties[0].rho[1] + NewCellProperties[1].rho[1] );
		NewCellProperties[0].rho[1] *= fWeight;
		NewCellProperties[1].rho[1] *= fWeight;

		if( NewCellProperties[0].rho_v[1] || NewCellProperties[1].rho_v[1] )
		{
                    temp1 = ( 2.0 * CellProperties.rho_v[1] ) - ( NewCellProperties[0].rho_v[1] + NewCellProperties[1].rho_v[1] );
                    temp2 = fabs( NewCellProperties[0].rho_v[1] );
                    temp3 = fabs( NewCellProperties[1].rho_v[1] );
                    temp4 = temp1 / ( temp2 + temp3 );
                    NewCellProperties[0].rho_v[1] += ( temp2 * temp4 );
                    NewCellProperties[1].rho_v[1] += ( temp3 * temp4 );
		}

		for( j=0; j<SPECIES; j++ )
		{
		    fWeight = ( 2.0 * CellProperties.TE_KE[1][j] ) / ( NewCellProperties[0].TE_KE[1][j] + NewCellProperties[1].TE_KE[1][j] );
		    NewCellProperties[0].TE_KE[1][j] *= fWeight;
		    NewCellProperties[1].TE_KE[1][j] *= fWeight;
		}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
		fWeight = ( 2.0 * CellProperties.rho_e ) / ( NewCellProperties[0].rho_e + NewCellProperties[1].rho_e );
		NewCellProperties[0].rho_e *= fWeight;
		NewCellProperties[1].rho_e *= fWeight;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // ENFORCE_CONSERVATION

                pNewCell[0] = new CAdaptiveMeshCell( &(NewCellProperties[0]) );
		pNewCell[1] = new CAdaptiveMeshCell( &(NewCellProperties[1]) );
                
		if( pLeftCell )
		{
                    pLeftCell->SetPointer( RIGHT, pNewCell[0] );
                    pNewCell[0]->SetPointer( LEFT, pLeftCell );
		}
		else
		{
                    pNewCell[0]->SetPointer( LEFT, NULL );
                    pStartOfCurrentRow = pNewCell[0];
		}

		pNewCell[0]->SetPointer( RIGHT, pNewCell[1] );
		pNewCell[1]->SetPointer( LEFT, pNewCell[0] );

		pNewCell[1]->SetPointer( RIGHT, pRightCell );
		pRightCell->SetPointer( LEFT, pNewCell[1] );

#ifdef NON_EQUILIBRIUM_RADIATION
                if( CellProperties.pIonFrac )
                    delete CellProperties.pIonFrac;
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
        	delete CellProperties.pKinetic;
#endif // USE_KINETIC_MODEL

		delete pActiveCell;

		pNextActiveCell = pNewCell[1];
            }
            else
            {
                // Refine the right-hand cell

		// First shift the left-hand and active cells to the right
		pLeftCell = pActiveCell;
		pLeftCell->GetCellProperties( &LeftCellProperties );
		pActiveCell = pRightCell;
		pActiveCell->GetCellProperties( &CellProperties );

		ZeroCellProperties( &(NewCellProperties[0]) );
		ZeroCellProperties( &(NewCellProperties[1]) );

#ifdef NON_EQUILIBRIUM_RADIATION
                if( CellProperties.pIonFrac )
                {
                    // Create a new ionfrac object for each new cell
                    NewCellProperties[0].pIonFrac = new CIonFrac( CellProperties.pIonFrac, 0, pRadiation );
                    NewCellProperties[1].pIonFrac = new CIonFrac( CellProperties.pIonFrac, 0, pRadiation );
                }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
		// Create a new kinetic object for each new cell
		NewCellProperties[0].pKinetic = new CKinetic();
		NewCellProperties[1].pKinetic = new CKinetic();
#endif // USE_KINETIC_MODEL

		NewCellProperties[0].iRefinementLevel = CellProperties.iRefinementLevel + 1;
		NewCellProperties[1].iRefinementLevel = NewCellProperties[0].iRefinementLevel;

		for( j=1; j<=CellProperties.iRefinementLevel; j++ )
		{
                    NewCellProperties[0].iUniqueID[j] = CellProperties.iUniqueID[j];
                    NewCellProperties[1].iUniqueID[j] = CellProperties.iUniqueID[j];
		}

		NewCellProperties[0].cell_width = 0.5 * CellProperties.cell_width;
		NewCellProperties[1].cell_width = NewCellProperties[0].cell_width;

		NewCellProperties[0].s[0] = CellProperties.s[0];
		NewCellProperties[0].s[1] = NewCellProperties[0].s[0] + ( 0.5 * NewCellProperties[0].cell_width );
		NewCellProperties[0].s[2] = CellProperties.s[1];
		NewCellProperties[1].s[0] = CellProperties.s[1];
		NewCellProperties[1].s[1] = NewCellProperties[1].s[0] + ( 0.5 * NewCellProperties[0].cell_width );
		NewCellProperties[1].s[2] = CellProperties.s[2];

		pRightCell = pActiveCell->pGetPointer( RIGHT );

                if( pRightCell )
		{
                    pRightCell->GetCellProperties( &RightCellProperties );

                    do {

                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] = rand();

                    } while ( NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == LeftCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] ||
                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == RightCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] );

                    NewCellProperties[1].iUniqueID[NewCellProperties[1].iRefinementLevel] = NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel];

                    pFarLeftCell = pLeftCell->pGetPointer( LEFT );
                    pFarRightCell = pRightCell->pGetPointer( RIGHT );
			
                    if( pFarLeftCell && pFarRightCell )
                    {
			pFarLeftCell->GetCellProperties( &FarLeftCellProperties );
			pFarRightCell->GetCellProperties( &FarRightCellProperties );

			// NEW CELLS

			x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];
			x[5] = FarRightCellProperties.s[1];

			y[1] = FarLeftCellProperties.rho[1];
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
			y[5] = FarRightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]), &error );
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[1] = FarLeftCellProperties.rho_v[1];
			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
			y[5] = FarRightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]) );
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]), &error );
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[1] = FarLeftCellProperties.TE_KE[1][j];
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
                            y[5] = FarRightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                            LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
#else
                            FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]), &error );
                            FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac && NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[1] = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[5] = FarRightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 2, NewCellProperties[0].s[1] );
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[2]), &(ppIonFrac[2]), 2, NewCellProperties[1].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 4, NewCellProperties[0].s[1] );
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 4, NewCellProperties[1].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = FarLeftCellProperties.rho_e;
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
			y[5] = FarRightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e), &error );
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
                    }
                    else if( !pFarLeftCell )
                    {
// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// LEFT-HAND BOUNDARY

                        // Implement hydrostatic boundary conditions in the left-most of the two new cells

                        x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];

                        y[1] = LeftCellProperties.rho[1];
                        y[2] = CellProperties.rho[1];
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );

                        NewCellProperties[0].rho_v[1] = 0.0;

                        for( j=0; j<SPECIES; j++ )
                        {
                            y[1] = LeftCellProperties.TE_KE[1][j];
                            y[2] = CellProperties.TE_KE[1][j];
                            LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                        }

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac )
                        {
                            ppIonFrac[1] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = CellProperties.pIonFrac->ppGetIonFrac();
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[0].s[1] );
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
                        y[1] = LeftCellProperties.rho_e;
                        y[2] = CellProperties.rho_e;
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

			pFarRightCell->GetCellProperties( &FarRightCellProperties );

			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];
			x[5] = FarRightCellProperties.s[1];
						
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
			y[5] = FarRightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
			y[5] = FarRightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
                            y[5] = FarRightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
#else
                            FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[5] = FarRightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[2]), &(ppIonFrac[2]), 2, NewCellProperties[1].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 4, NewCellProperties[1].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
			y[5] = FarRightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[2]), &(y[2]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#else
			FitPolynomial4( &(x[1]), &(y[1]), NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
                    }
                    else if( !pFarRightCell )
                    {
			pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

			x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			x[3] = CellProperties.s[1];
			x[4] = RightCellProperties.s[1];

                        y[1] = FarLeftCellProperties.rho[1];
			y[2] = LeftCellProperties.rho[1];
			y[3] = CellProperties.rho[1];
			y[4] = RightCellProperties.rho[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]), &error );
#endif // LINEAR_RESTRICTION

			y[1] = FarLeftCellProperties.rho_v[1];
			y[2] = LeftCellProperties.rho_v[1];
			y[3] = CellProperties.rho_v[1];
			y[4] = RightCellProperties.rho_v[1];
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_v[1]), &error );
#endif // LINEAR_RESTRICTION

			for( j=0; j<SPECIES; j++ )
			{
                            y[1] = FarLeftCellProperties.TE_KE[1][j];
                            y[2] = LeftCellProperties.TE_KE[1][j];
                            y[3] = CellProperties.TE_KE[1][j];
                            y[4] = RightCellProperties.TE_KE[1][j];
#ifdef LINEAR_RESTRICTION
                            LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
#else
                            FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]), &error );
#endif // LINEAR_RESTRICTION
			}

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[0].pIonFrac )
                        {
                            ppIonFrac[1] = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[3] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[4] = RightCellProperties.pIonFrac->ppGetIonFrac();
#ifdef LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( &(x[1]), &(ppIonFrac[1]), 2, NewCellProperties[0].s[1] );
#else // LINEAR_RESTRICTION
                            NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 4, NewCellProperties[0].s[1] );
#endif // LINEAR_RESTRICTION
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE

			y[1] = FarLeftCellProperties.rho_e;
			y[2] = LeftCellProperties.rho_e;
			y[3] = CellProperties.rho_e;
			y[4] = RightCellProperties.rho_e;
#ifdef LINEAR_RESTRICTION
			LinearFit( &(x[1]), &(y[1]), NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
#else
			FitPolynomial4( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e), &error );
#endif // LINEAR_RESTRICTION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// RIGHT-HAND BOUNDARY

                        // Implement hydrostatic boundary conditions in the right-most of the two new cells

                        x[1] = CellProperties.s[1];
			x[2] = RightCellProperties.s[1];

                        y[1] = CellProperties.rho[1];
                        y[2] = RightCellProperties.rho[1];
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );

                        NewCellProperties[1].rho_v[1] = 0.0;

                        for( j=0; j<SPECIES; j++ )
                        {
                            y[1] = CellProperties.TE_KE[1][j];
                            y[2] = RightCellProperties.TE_KE[1][j];
                            LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
                        }

#ifdef NON_EQUILIBRIUM_RADIATION
                        if( NewCellProperties[1].pIonFrac )
                        {
                            ppIonFrac[1] = CellProperties.pIonFrac->ppGetIonFrac();
                            ppIonFrac[2] = RightCellProperties.pIonFrac->ppGetIonFrac();
                            NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[1].s[1] );
                        }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
			y[1] = CellProperties.rho_e;
                        y[2] = RightCellProperties.rho_e;
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
                    }
		}
		else
		{
                    do {

                        NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] = rand();

                    } while ( NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel] == LeftCellProperties.iUniqueID[NewCellProperties[0].iRefinementLevel] );

                    NewCellProperties[1].iUniqueID[NewCellProperties[1].iRefinementLevel] = NewCellProperties[0].iUniqueID[NewCellProperties[0].iRefinementLevel];

// ******************************************************************************
// *    IMPLEMENT HYDROSTATIC BOUNDARY CONDITIONS                               *
// ******************************************************************************

// RIGHT-HAND BOUNDARY
                    
                    // NEW CELLS
                    
		    x[1] = LeftCellProperties.s[1];
		    x[2] = CellProperties.s[1];

                    y[1] = LeftCellProperties.rho[1];
                    y[2] = CellProperties.rho[1];
                    LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho[1]) );
                    LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho[1]) );

                    NewCellProperties[0].rho_v[1] = 0.0;
                    NewCellProperties[1].rho_v[1] = 0.0;

                    for( j=0; j<SPECIES; j++ )
                    {
                        y[1] = LeftCellProperties.TE_KE[1][j];
                        y[2] = CellProperties.TE_KE[1][j];
                        LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].TE_KE[1][j]) );
                        LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].TE_KE[1][j]) );
                    }

#ifdef NON_EQUILIBRIUM_RADIATION
                    if( NewCellProperties[0].pIonFrac && NewCellProperties[1].pIonFrac )
                    {
                        ppIonFrac[1] = LeftCellProperties.pIonFrac->ppGetIonFrac();
                        ppIonFrac[2] = CellProperties.pIonFrac->ppGetIonFrac();
                        NewCellProperties[0].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[0].s[1] );
                        NewCellProperties[1].pIonFrac->InterpolateAllIonFrac( x, ppIonFrac, 2, NewCellProperties[1].s[1] );
                    }
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
                    y[1] = LeftCellProperties.rho_e;
                    y[2] = CellProperties.rho_e;
                    LinearFit( x, y, NewCellProperties[0].s[1], &(NewCellProperties[0].rho_e) );
                    LinearFit( x, y, NewCellProperties[1].s[1], &(NewCellProperties[1].rho_e) );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
		}

// ******************************************************************************
// *    ENSURE MASS AND ENERGY IS CONSERVED BY THE RESTRICTION OPERATOR         *
// ******************************************************************************
#ifdef ENFORCE_CONSERVATION
		// A correction is applied to the mass and energy densities in the event that the sum of the integrated quantities in the new cells does not equal the integrated quantity in the original cell

		fWeight = ( 2.0 * CellProperties.rho[1] ) / ( NewCellProperties[0].rho[1] + NewCellProperties[1].rho[1] );
		NewCellProperties[0].rho[1] *= fWeight;
		NewCellProperties[1].rho[1] *= fWeight;

		if( NewCellProperties[0].rho_v[1] || NewCellProperties[1].rho_v[1] )
		{
                    temp1 = ( 2.0 * CellProperties.rho_v[1] ) - ( NewCellProperties[0].rho_v[1] + NewCellProperties[1].rho_v[1] );
                    temp2 = fabs( NewCellProperties[0].rho_v[1] );
                    temp3 = fabs( NewCellProperties[1].rho_v[1] );
                    temp4 = temp1 / ( temp2 + temp3 );
                    NewCellProperties[0].rho_v[1] += ( temp2 * temp4 );
                    NewCellProperties[1].rho_v[1] += ( temp3 * temp4 );
		}

		for( j=0; j<SPECIES; j++ )
		{
		    fWeight = ( 2.0 * CellProperties.TE_KE[1][j] ) / ( NewCellProperties[0].TE_KE[1][j] + NewCellProperties[1].TE_KE[1][j] );
		    NewCellProperties[0].TE_KE[1][j] *= fWeight;
		    NewCellProperties[1].TE_KE[1][j] *= fWeight;
		}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
		fWeight = ( 2.0 * CellProperties.rho_e ) / ( NewCellProperties[0].rho_e + NewCellProperties[1].rho_e );
		NewCellProperties[0].rho_e *= fWeight;
		NewCellProperties[1].rho_e *= fWeight;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // ENFORCE_CONSERVATION

		pNewCell[0] = new CAdaptiveMeshCell( &(NewCellProperties[0]) );
		pNewCell[1] = new CAdaptiveMeshCell( &(NewCellProperties[1]) );
                    
		pLeftCell->SetPointer( RIGHT, pNewCell[0] );
		pNewCell[0]->SetPointer( LEFT, pLeftCell );

		pNewCell[0]->SetPointer( RIGHT, pNewCell[1] );
		pNewCell[1]->SetPointer( LEFT, pNewCell[0] );
			
		if( pRightCell )
		{
                    pNewCell[1]->SetPointer( RIGHT, pRightCell );
                    pRightCell->SetPointer( LEFT, pNewCell[1] );
		}
		else
		{
                    pNewCell[1]->SetPointer( RIGHT, NULL );
		}

#ifdef NON_EQUILIBRIUM_RADIATION
                if( CellProperties.pIonFrac )
                    delete CellProperties.pIonFrac;
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
		delete CellProperties.pKinetic;
#endif // USE_KINETIC_MODEL

		delete pActiveCell;

		pNextActiveCell = pNewCell[1];
            }
	}
	else
	{
            pNextActiveCell = pRightCell;
	}
    }

} while( iRestricted );
}
#endif // ADAPT

void CAdaptiveMesh::Solve( void )
{
clock_t timer[3];
double fNextOutputTime;
int iOutputStepCount = 0;
#ifdef ADAPT
	int iAdaptiveStepCount = 0;
	bool bAdapt;
#endif // ADAPT

printf( "\nSolving...\n\n" );

fNextOutputTime = mesh_time + Params.OutputPeriod;

// Now time-step the mesh
timer[0] = timer[1] = clock();
while( mesh_time <= Params.Duration )
{
    iOutputStepCount++;
    if( iOutputStepCount == OUTPUT_EVERY_N_TIME_STEPS )
    {
		timer[2] =clock();
		ShowProgress( timer );
		iOutputStepCount = 0;
		timer[1] = clock();
    }

#ifdef ADAPT
	iAdaptiveStepCount++;
	bAdapt = false;
	if( iAdaptiveStepCount == ADAPT_EVERY_N_TIME_STEPS )
	{
		bAdapt = true;
		iAdaptiveStepCount = 0;
	}
	Integrate( bAdapt );
#else // ADAPT
	Integrate();
#endif // ADAPT
	
    // If the output period has elapsed since the last output then write the profiles to a file
    if( mesh_time >= fNextOutputTime )
    {
        WriteToFile();
		fNextOutputTime += Params.OutputPeriod;
    }
}

printf( "\nDone!\n\n" );
}

#ifdef ADAPT
void CAdaptiveMesh::Integrate( bool bAdapt )
#else // ADAPT
void CAdaptiveMesh::Integrate( void )
#endif // ADAPT
{
PCELL pNextActiveCell, pPreviousCell = NULL, pNewCell;
CELLPROPERTIES CellProperties, NewCellProperties;
double delta_t;

#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
int iCell = 0;
#endif // OPENMP || USE_KINETIC_MODEL

// The integration is 2nd order in time and uses the derivatives calculated at half the time-step
delta_t = 0.5 * mesh_delta_t;

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &NewCellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
    // Create a new ionfrac object for the new cell
    NewCellProperties.pIonFrac = new CIonFrac( NewCellProperties.pIonFrac, 0, pRadiation );
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_KINETIC_MODEL
    NewCellProperties.pKinetic = new CKinetic();
#endif // USE_KINETIC_MODEL

    Half_Time_Step( &NewCellProperties, delta_t );

    // Create a new mesh cell at the current s but one time step advanced
    pNewCell = new CAdaptiveMeshCell( &NewCellProperties );

    // Set the RIGHT pointer to NULL
    pNewCell->SetPointer( RIGHT, NULL );
		
    // Set the BOTTOM pointer to the active cell (the cell from which the new cell has been derived)
    pNewCell->SetPointer( BOTTOM, pActiveCell );

    // If this is the left-most cell then set the LEFT pointer to NULL
    if( !( pActiveCell->pGetPointer( LEFT ) ) )
    {
        pNewCell->SetPointer( LEFT, NULL );
	
	// Also, set the pointer to the left-most cell at the current time to this pointer and change
	// the pointer to the start of the previous row accordingly
	pStartOfPreviousRow = pStartOfCurrentRow;
	pStartOfCurrentRow = pNewCell;
    }
    // This cell is not the left-most and so it will have a LEFT pointer
    else
    {
        // Set the RIGHT pointer of the previous cell to the new cell
	pPreviousCell->SetPointer( RIGHT, pNewCell );
	pNewCell->SetPointer( LEFT, pPreviousCell );
    }

    pPreviousCell = pNewCell;
    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );

#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
    ppCellList[iCell] = pNewCell;
    iCell++;
#endif // OPENMP || USE_KINETIC_MODEL
}

// Calculate the physical quantities in the new cells
CalculatePhysicalQuantities();

// Evaluate the terms of the equations
EvaluateTerms( mesh_time + delta_t, &delta_t, FALSE );

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    Full_Time_Step( &CellProperties, mesh_delta_t );

    pActiveCell->UpdateCellProperties( &CellProperties );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

#ifdef ADAPT
// Adapt the mesh
if( bAdapt )
{
	Adapt();
#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
	CreateIndexedCellList();
#endif // OPENMP || USE_KINETIC_MODEL
}
#endif // ADAPT

// Calculate the physical quantities
CalculatePhysicalQuantities();

// Evaluate the terms of the equations
EvaluateTerms( mesh_time + mesh_delta_t, &delta_t, TRUE );

// Limit time-step increases to be within a specified fractional amount of the previous time-step
if( delta_t > mesh_delta_t )
{
    if( delta_t / mesh_delta_t > TIME_STEP_INCREASE_LIMIT )
        delta_t = TIME_STEP_INCREASE_LIMIT * mesh_delta_t;
}

// The mesh has been time-stepped
mesh_time += mesh_delta_t;
mesh_delta_t = delta_t;

// Free the previous row in order to conserve memory
FreePreviousRow();
}

void CAdaptiveMesh::FreeCurrentRow( void )
{
#if defined(NON_EQUILIBRIUM_RADIATION) || defined(USE_KINETIC_MODEL)
	CELLPROPERTIES CellProperties;
#endif // NON_EQUILIBRIUM_RADIATION || USE_KINETIC_MODEL
PCELL pNextActiveCell;

pNextActiveCell = pStartOfCurrentRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );

#if defined(NON_EQUILIBRIUM_RADIATION) || defined(USE_KINETIC_MODEL)
	pActiveCell->GetCellProperties( &CellProperties );
#ifdef NON_EQUILIBRIUM_RADIATION
	delete CellProperties.pIonFrac;
#endif // NON_EQUILIBRIUM_RADIATION
#ifdef USE_KINETIC_MODEL
    delete CellProperties.pKinetic;
#endif // USE_KINETIC_MODEL
#endif // NON_EQUILIBRIUM_RADIATION || USE_KINETIC_MODEL

    delete pActiveCell;
}
}

void CAdaptiveMesh::FreePreviousRow( void )
{
#if defined(NON_EQUILIBRIUM_RADIATION) || defined(USE_KINETIC_MODEL)
	CELLPROPERTIES CellProperties;
#endif // NON_EQUILIBRIUM_RADIATION || USE_KINETIC_MODEL
PCELL pNextActiveCell;

pNextActiveCell = pStartOfPreviousRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );

#if defined(NON_EQUILIBRIUM_RADIATION) || defined(USE_KINETIC_MODEL)
	pActiveCell->GetCellProperties( &CellProperties );
#ifdef NON_EQUILIBRIUM_RADIATION
	delete CellProperties.pIonFrac;
#endif // NON_EQUILIBRIUM_RADIATION
#ifdef USE_KINETIC_MODEL
    delete CellProperties.pKinetic;
#endif // USE_KINETIC_MODEL
#endif // NON_EQUILIBRIUM_RADIATION || USE_KINETIC_MODEL

    delete pActiveCell;
}
}

void CAdaptiveMesh::ShowProgress( clock_t *ptimer )
{
	void FindTemporalUnits( double *pfTime, int *piUnits );
	double fTime[4];
	int iUnits[4], i;
	char cUnitLabel[6] = {'s','m','h','d','w','y'};

#ifdef OPENMP
	int iDivisor;
	#pragma omp parallel
		iDivisor = CLOCKS_PER_SEC * omp_get_num_threads();
	fTime[0] =  ((double)(ptimer[2]-ptimer[1])) / iDivisor;
	fTime[1] = ((double)(clock() - ptimer[0])) / iDivisor;
#else // OPENMP
	fTime[0] =  ((double)(ptimer[2]-ptimer[1])) / CLOCKS_PER_SEC;
	fTime[1] = ((double)(clock() - ptimer[0])) / CLOCKS_PER_SEC;
#endif // OPENMP	
	fTime[2] = fTime[1] / mesh_time;
	fTime[3] = ( Params.Duration - mesh_time ) * fTime[2];
	// Convert the times into convenient units given their magnitudes
	for( i=0; i<4; i++ )
		FindTemporalUnits( &(fTime[i]), &(iUnits[i]) );

	printf( "model-time elapsed = %.4e s; dt = %.4e s\n", mesh_time, mesh_delta_t );	
	printf( "wall-time elapsed = %.4e %c; ~ wall-time remaining = %.4e %c\n", fTime[1], cUnitLabel[iUnits[1]], fTime[3], cUnitLabel[iUnits[3]] );
	printf( "\twall-time/%i steps = %.4e %c; wall-time/second = %.4e %c\n",  OUTPUT_EVERY_N_TIME_STEPS, fTime[0], cUnitLabel[iUnits[0]], fTime[2], cUnitLabel[iUnits[2]] );
	printf( "\trefinement level = %i/%i; total grid cells = %i\n\n", GetMaxRL(), MAX_REFINEMENT_LEVEL, Params.iNumberOfCells );
	/*
	printf( "model-time elapsed = %.4e s; wall-time elapsed = %.4e %c\n", mesh_time, fTime[1], cUnitLabel[iUnits[1]] );
	printf( "\tdt = %.4e s; refinement level = %i/%i; total grid cells = %i;\n", mesh_delta_t, GetMaxRL(), MAX_REFINEMENT_LEVEL, GetNumCells() );
    	printf( "\twall-time/%i steps = %.4e %c; ",  OUTPUT_EVERY_N_TIME_STEPS, fTime[0], cUnitLabel[iUnits[0]] );
	printf( "wall-time/second = %.4e %c; ", fTime[2], cUnitLabel[iUnits[2]] );
	printf( "~ wall-time remaining = %.4e %c\n\n", fTime[3], cUnitLabel[iUnits[3]] );
	*/
}
		void FindTemporalUnits( double *pfTime, int *piUnits )
		{
			// The time is initally in units of seconds
			*piUnits = 0;	// seconds
			if( *pfTime < 60.0 ) return;

			// Convert the time into minutes
			*piUnits = 1;	// minutes
			*pfTime = *pfTime / 60.0;
			if( *pfTime < 60.0 ) return;
			
			// Convert the time into hours
			*piUnits = 2; // hours
			*pfTime = *pfTime / 60.0;
			if( *pfTime < 24.0 ) return;
			
			// Convert the time into days
			*piUnits = 3;	// days
			*pfTime = *pfTime / 24.0;
			if( *pfTime < 7.0 ) return;
			
			// Convert the time into weeks
			*piUnits = 4;	// weeks
			*pfTime = *pfTime / 7.0;
			if( *pfTime < 52.0 ) return;
			
			// Convert the time into years
			*piUnits = 5;	// years
			*pfTime = *pfTime / 52.0;
			return;
		}

void CAdaptiveMesh::WriteToFile( void )
{
FILE *pAMRFile;
char szAMRFilename[256];

#ifdef WRITE_FILE_PHYSICAL
FILE *pPhysicalFile;
char szPhysicalFilename[256];
sprintf( szPhysicalFilename, "Results/profile%i.phy", iFileNumber );
pPhysicalFile = fopen( szPhysicalFilename, "w" );
#endif // WRITE_FILE_PHYSICAL

#ifdef NON_EQUILIBRIUM_RADIATION
#ifdef WRITE_FILE_ION_POPULATIONS
FILE *pNEqIonFile;
char szNEqIonFilename[256];
sprintf( szNEqIonFilename, "Results/profile%i.ine", iFileNumber );
pNEqIonFile = fopen( szNEqIonFilename, "w" );
#endif // WRITE_FILE_ION_POPULATIONS
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef WRITE_FILE_SCALES
FILE *pScaleFile;
char szScaleFilename[256];
sprintf( szScaleFilename, "Results/profile%i.scl", iFileNumber );
pScaleFile = fopen( szScaleFilename, "w" );
#endif // WRITE_FILE_SCALES

#ifdef WRITE_FILE_TERMS
FILE *pTermsFile;
char szTermsFilename[256];
int iTerm;
sprintf( szTermsFilename, "Results/profile%i.trm", iFileNumber );
pTermsFile = fopen( szTermsFilename, "w" );
#endif // WRITE_FILE_TERMS

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
#ifdef WRITE_FILE_HSTATE
FILE *pHStateFile;
char szHStateFilename[256];
sprintf( szHStateFilename, "Results/profile%i.Hstate", iFileNumber );
pHStateFile = fopen( szHStateFilename, "w" );
#endif // WRITE_FILE_HSTATE
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;

int j;

// Write data into the files

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

#ifdef WRITE_FILE_PHYSICAL
    fprintf( pPhysicalFile, "%.8e\t%.8e\t%.8e", CellProperties.s[1], CellProperties.v[1], CellProperties.Cs );

    for( j=0; j<SPECIES; j++ )
        fprintf( pPhysicalFile, "\t%.8e", CellProperties.n[j] );

    for( j=0; j<SPECIES; j++ )
        fprintf( pPhysicalFile, "\t%.8e", CellProperties.P[1][j] );

    for( j=0; j<SPECIES; j++ )
        fprintf( pPhysicalFile, "\t%.8e", CellProperties.T[j] );

    for( j=0; j<SPECIES; j++ )
        fprintf( pPhysicalFile, "\t%.8e", CellProperties.Fc[1][j] );

    fprintf( pPhysicalFile, "\n" );
#endif // WRITE_FILE_PHYSICAL

#ifdef NON_EQUILIBRIUM_RADIATION
#ifdef WRITE_FILE_ION_POPULATIONS
    fprintf( pNEqIonFile, "%.8e", CellProperties.s[1] );
    CellProperties.pIonFrac->WriteAllIonFracToFile( pNEqIonFile );
#endif // WRITE_FILE_ION_POPULATIONS
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef WRITE_FILE_SCALES
    fprintf( pScaleFile, "%.8e\t%.8e\t%.8e", CellProperties.s[1], CellProperties.cell_width, CellProperties.advection_delta_t );

    for( j=0; j<SPECIES; j++ )
        fprintf( pScaleFile, "\t%.8e", CellProperties.conduction_delta_t[j] );

    fprintf( pScaleFile, "\t%.8e", CellProperties.viscosity_delta_t );
    fprintf( pScaleFile, "\t%.8e", CellProperties.collision_delta_t );
    fprintf( pScaleFile, "\t%.8e", CellProperties.radiation_delta_t );
#ifdef NON_EQUILIBRIUM_RADIATION
    fprintf( pScaleFile, "\t%.8e", CellProperties.atomic_delta_t );
#endif // NON_EQUILIBRIUM_RADIATION
    fprintf( pScaleFile, "\n" );
#endif // WRITE_FILE_SCALES

#ifdef WRITE_FILE_TERMS
    fprintf( pTermsFile, "%.8e\n", CellProperties.s[1] );

    fprintf( pTermsFile, "%.8e", CellProperties.drhobydt );
    for( iTerm=0; iTerm<MASS_TERMS; iTerm++ )
        fprintf( pTermsFile, "\t%.8e", CellProperties.rho_term[iTerm] );

    fprintf( pTermsFile, "\n" );

    fprintf( pTermsFile, "%.8e", CellProperties.drho_vbydt );
    for( iTerm=0; iTerm<MOMENTUM_TERMS; iTerm++ )
        fprintf( pTermsFile, "\t%.8e", CellProperties.rho_v_term[iTerm] );

    fprintf( pTermsFile, "\n" );

    for( j=0; j<SPECIES; j++ )
    {
        fprintf( pTermsFile, "%.8e", CellProperties.dTE_KEbydt[j] );
	for( iTerm=0; iTerm<ENERGY_TERMS; iTerm++ )
            fprintf( pTermsFile, "\t%.8e", CellProperties.TE_KE_term[iTerm][j] );

	fprintf( pTermsFile, "\n" );
    }
#endif // WRITE_FILE_TERMS

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
#ifdef WRITE_FILE_HSTATE
    fprintf( pHStateFile, "%.8e", CellProperties.s[1] );
    for( j=0; j<6; j++ )
	fprintf( pHStateFile, "\t%.8e", CellProperties.Hstate[j] );
    fprintf( pHStateFile, "\n" );
#endif // WRITE_FILE_HSTATE
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

#ifdef WRITE_FILE_PHYSICAL
fclose( pPhysicalFile );
#endif // WRITE_FILE_PHYSICAL

#ifdef NON_EQUILIBRIUM_RADIATION
#ifdef WRITE_FILE_ION_POPULATIONS
fclose( pNEqIonFile );
#endif // WRITE_FILE_ION_POPULATIONS
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef WRITE_FILE_SCALES
fclose( pScaleFile );
#endif // WRITE_FILE_SCALES

#ifdef WRITE_FILE_TERMS
fclose( pTermsFile );
#endif // WRITE_FILE_TERMS

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
#ifdef WRITE_FILE_HSTATE
fclose( pHStateFile );
#endif // WRITE_FILE_HSTATE
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

// Create the .amr file so that the simulation can be continued from the current output if necessary

sprintf( szAMRFilename, "Results/profile%i.amr", iFileNumber );

pAMRFile = fopen( szAMRFilename, "w" );

fprintf( pAMRFile, "%g\n%i\n%.16e\n%i\n", mesh_time, iFileNumber, Params.L, Params.iNumberOfCells );

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    fprintf( pAMRFile, "%.16e\t%.16e\t", CellProperties.s[1], CellProperties.cell_width );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
    fprintf( pAMRFile, "%.16e\t", CellProperties.rho_e );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

    fprintf( pAMRFile, "%.16e\t%.16e", CellProperties.rho[1], CellProperties.rho_v[1] );

    for( j=0; j<SPECIES; j++ )
        fprintf( pAMRFile, "\t%.16e", CellProperties.TE_KE[1][j] );

    fprintf( pAMRFile, "\t%i", CellProperties.iRefinementLevel );
    for( j=1; j<=MAX_REFINEMENT_LEVEL; j++ )
        fprintf( pAMRFile, "\t%i", CellProperties.iUniqueID[j] );

    fprintf( pAMRFile, "\n" );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

fclose( pAMRFile );

#ifdef USE_KINETIC_MODEL
FILE *pDFNFile;
char szDFNFilename[256];
double *pupsilon, *pMaxDFN_ee, *pNonMaxDFN;

// Write the distribution functions to a file
sprintf( szDFNFilename, "Results/profile%i.dfn", iFileNumber );

pDFNFile = fopen( szDFNFilename, "w" );

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    // Write the grid cell locations along the top row of the file
    fprintf( pDFNFile, "%.8e\t", CellProperties.s[1] );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// Write the thermal speeds along the next row of the file
pNextActiveCell = pStartOfCurrentRow;
pActiveCell = pNextActiveCell;
pActiveCell->GetCellProperties( &CellProperties );
pupsilon = CellProperties.pKinetic->Get_pupsilon();
fprintf( pDFNFile, "\n" );
for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
    fprintf( pDFNFile, "%.8e\t", pupsilon[j] );

pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
fprintf( pDFNFile, "\n" );
for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
    fprintf( pDFNFile, "%.8e\t", log10( pMaxDFN_ee[j] ) );

pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
pActiveCell = pNextActiveCell;
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
fprintf( pDFNFile, "\n" );
for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
    fprintf( pDFNFile, "%.8e\t", log10( pMaxDFN_ee[j] ) );

pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );
	
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    fprintf( pDFNFile, "\n" );
    for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
        fprintf( pDFNFile, "%.8e\t", log10( pNonMaxDFN[j] ) );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

pActiveCell = pNextActiveCell;
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
fprintf( pDFNFile, "\n" );
for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
    fprintf( pDFNFile, "%.8e\t", log10( pMaxDFN_ee[j] ) );

pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
pActiveCell = pNextActiveCell;
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
fprintf( pDFNFile, "\n" );
for( j=0; j<DISTRIBUTION_DATA_POINTS; j++ )
    fprintf( pDFNFile, "%.8e\t", log10( pMaxDFN_ee[j] ) );

fclose( pDFNFile );
#endif // USE_KINETIC_MODEL

iFileNumber++;
}
