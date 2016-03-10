// ****
// *
// * Function bodies for the class definition of the 
// * time-dependent hydrodynamic equations, inherited by the 
// * adaptive mesh class
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/03/2016
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "eqns.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"


double fLogLambda_ei( double Te, double Ti, double n )
{
double limit;

Te *= 1.29199251618e-4;
Ti *= 1.29199251618e-4;

limit = Ti * ( ELECTRON_MASS / AVERAGE_PARTICLE_MASS );

if( limit < Te && Te < 10.0 )
{
    return 23.0 - log( sqrt( n ) * pow( Te, (-1.5) ) );
}
else if( limit < 10.0 && 10.0 < Te )
{
    return 24.0 - log( sqrt( n ) * ( 1.0 / Te ) );
}
else // if( Te < limit )
{
    return 30.0 - log( sqrt( n ) * pow( Ti, (-1.5) ) );
}
}

double fLogLambda_ii( double Ti, double n )
{
double term1;

// Assume that proton-proton collisions dominate ion-ion interactions
// sqrt(2.0) * Z^3 = 1.4142135623730950488016887242097 * 1.0 * 1.0 * 1.0
term1 = 1.4142135623730950488016887242097 * ( sqrt( n ) / pow( Ti, 1.5 ) );
return 23.0 - log( term1 );
}


// Constructor
CEquations::CEquations( void )
{
Initialise();
}

// Destructor
CEquations::~CEquations( void )
{
FreeAll();
}

void CEquations::Initialise( void )
{
FILE *pFile;
double fTemp;
int i, iTemp;

#ifdef USE_KINETIC_MODEL
ppCellList = NULL;
#endif // USE_KINETIC_MODEL

// Get the HYDRAD initial conditions
pFile = fopen( "HYDRAD/config/HYDRAD.cfg", "r" );
// Get the initial profiles
fscanf( pFile, "%s", Params.Profiles );
// Get the gravity look-up table filename
fscanf( pFile, "%s", Params.GravityFilename );
#ifdef USE_TABULATED_CROSS_SECTION
// Get the cross-section look-up table filename
fscanf( pFile, "%s", Params.CrossSectionFilename );
#endif // USE_TABULATED_CROSS_SECTION
// Get the duration
ReadDouble( pFile, &Params.Duration );
// Get the output period
ReadDouble( pFile, &Params.OutputPeriod );
fclose( pFile );

// Get the loop length from the profiles file
pFile = fopen( Params.Profiles, "r" );
ReadDouble( pFile, &fTemp );
fscanf( pFile, "%i", &iTemp );
ReadDouble( pFile, &Params.L );
fclose( pFile );

#ifdef USE_KINETIC_MODEL
// Get the tabulated values from tables I and II (for Z = 1) in Spitzer & Harm, 1953, Phys. Rev., 89, 977
Get_SH_Table();
#endif // USE_KINETIC_MODEL

// Initialise the gravity look-up table
pFile = fopen( Params.GravityFilename, "r" );
fscanf( pFile, "%i", &igdp );
ppGravity = (double**)malloc( igdp * sizeof( double ) );
for( i=0; i<igdp; i++ )
{
    ppGravity[i] = (double*)malloc( 2 * sizeof( double ) );
    ReadDouble( pFile, &(ppGravity[i][0]) );
    ReadDouble( pFile, &(ppGravity[i][1]) );
}
fclose( pFile );

#ifdef USE_TABULATED_CROSS_SECTION
// Initialise the cross-section look-up table
pFile = fopen( Params.CrossSectionFilename, "r" );
fscanf( pFile, "%i", &icsdp );
ppCrossSection = (double**)malloc( icsdp * sizeof( double ) );
for( i=0; i<icsdp; i++ )
{
    ppCrossSection[i] = (double*)malloc( 2 * sizeof( double ) );
    ReadDouble( pFile, &(ppCrossSection[i][0]) );
    ReadDouble( pFile, &(ppCrossSection[i][1]) );
}
fclose( pFile );
#endif // USE_TABULATED_CROSS_SECTION

// Create the heating object and set the lower radiation temperature boundary
pHeat = new CHeat( (char *)"Heating_Model/config/heating_model.cfg", Params.L );

// Create the radiation objects
pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_neq.cfg" );
pRadiation2 = new CRadiation( (char *)"Radiation_Model/config/elements_eq.cfg" );
lower_radiation_temperature_boundary = MINIMUM_RADIATION_TEMPERATURE + ZERO_OVER_TEMPERATURE_INTERVAL;

#ifdef OPTICALLY_THICK_RADIATION
// Create the optically-thick ion objects
pHI = new COpticallyThickIon( 1, (char *)"h_1", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
pMgII = new COpticallyThickIon( 12, (char *)"mg_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
pCaII = new COpticallyThickIon( 20, (char *)"ca_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
#endif // OPTICALLY_THICK_RADIATION
}

void CEquations::FreeAll( void )
{
int i;

#ifdef USE_KINETIC_MODEL
free( ppCellList );
#endif // USE_KINETIC_MODEL

#ifdef USE_TABULATED_CROSS_SECTION
// Free the memory allocated to the cross-section look-up table
for( i=0; i<icsdp; i++ )
    free( ppCrossSection[i]);
free( ppCrossSection );
#endif // USE_TABULATED_CROSS_SECTION

// Free the memory allocated to the gravity look-up table
for( i=0; i<igdp; i++ )
    free( ppGravity[i]);
free( ppGravity );

// Delete the heating object
delete pHeat;

// Delete the radiation objects
delete pRadiation;
delete pRadiation2;

#ifdef OPTICALLY_THICK_RADIATION
// Delete the optically-thick ion objects
delete pHI;
delete pMgII;
delete pCaII;
#endif // OPTICALLY_THICK_RADIATION
}

void CEquations::CalculatePhysicalQuantities( void )
{
PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;

// Variables used for calculation of thermal conduction time-scale
double Kappa[SPECIES];
Kappa[ELECTRON] = SPITZER_ELECTRON_CONDUCTIVITY;
Kappa[HYDROGEN] = SPITZER_ION_CONDUCTIVITY;

// General variables
double term1, term2;
int j;

#ifdef OPTICALLY_THICK_RADIATION
double fHI;
#endif // OPTICALLY_THICK_RADIATION

pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

#ifdef OPTICALLY_THICK_RADIATION
    // Locate the apex cell from which the column densities will be calculated along each leg
    if( CellProperties.s[1] <= Params.L / 2.0 )
        pCentreOfCurrentRow = pActiveCell;
#endif // OPTICALLY_THICK_RADIATION

// ******************************************************************************
// *                                                                            *
// *    CALCULATE THE PHYSICAL QUANTITIES                                       *
// *                                                                            *
// ******************************************************************************

    CellProperties.n[HYDROGEN] = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;

    CellProperties.v[1] = CellProperties.rho_v[1] / CellProperties.rho[1];

    term1 = CellProperties.TE_KE[1][HYDROGEN] / CellProperties.rho[1];
    term2 = 0.5 * CellProperties.v[1] * CellProperties.v[1];
    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    CellProperties.T[HYDROGEN] = (4.830917874e15) * AVERAGE_PARTICLE_MASS * ( term1 - term2 );

#ifdef OPTICALLY_THICK_RADIATION
    if( CellProperties.T[HYDROGEN] < OPTICALLY_THICK_TEMPERATURE )
    {
        // 1.44e-4 provides some negligible number of electrons to avoid division by zero errors
        fHI = pHI->GetIonFrac( log10(CellProperties.T[HYDROGEN]) );
        // 1.000144 = 1.0 + 1.44e-4
        CellProperties.n[ELECTRON] = ( 1.000144 - fHI ) * CellProperties.n[HYDROGEN];
    }
    else
    {
        // 1.000144 = 1.0 + 1.44e-4
        CellProperties.n[ELECTRON] = 1.000144 * CellProperties.n[HYDROGEN];
    }
#else // OPTICALLY_THICK_RADIATION
    CellProperties.n[ELECTRON] = CellProperties.n[HYDROGEN];
#endif // OPTICALLY_THICK_RADIATION
    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    CellProperties.T[ELECTRON] = ( (4.830917874e15) * CellProperties.TE_KE[1][ELECTRON] ) / CellProperties.n[ELECTRON] ;

    // Temperature limiter
    for( j=0; j<SPECIES; j++ )
    {
        if( CellProperties.T[j] < MINIMUM_TEMPERATURE )
	{
            CellProperties.T[j] = MINIMUM_TEMPERATURE;
            // BOLTZMANN_CONSTANT / GAMMA_MINUS_ONE = 2.07e-16
            CellProperties.TE_KE[1][j] = (2.07e-16) * CellProperties.n[j] * CellProperties.T[j];

            if( j == HYDROGEN )
                CellProperties.TE_KE[1][j] += 0.5 * CellProperties.rho_v[1] * CellProperties.v[1];
        }
    }

    for( j=0; j<SPECIES; j++ )
    {
        CellProperties.P[1][j] = BOLTZMANN_CONSTANT * CellProperties.n[j] * CellProperties.T[j];
	CellProperties.TE_KE_P[1][j] = CellProperties.TE_KE[1][j] + CellProperties.P[1][j];
    }

    term1 = ( GAMMA * ( CellProperties.P[1][ELECTRON] + CellProperties.P[1][HYDROGEN] ) ) / CellProperties.rho[1];
    CellProperties.Cs = pow( term1, 0.5 );
    CellProperties.M = fabs( CellProperties.v[1] / CellProperties.Cs );

// ******************************************************************************
// *                                                                            *
// *    CALCULATE THE TIMESCALES                                                *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    COLLISIONS                                                              *
// ******************************************************************************

#ifdef FORCE_SINGLE_FLUID
    CellProperties.nu_ie = 1.0 / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
#else // FORCE_SINGLE_FLUID
    // Calculate the collision frequency between the electrons and ions
    // 4.820055089755540 * ELECTRON_MASS = 4.391070186767300e-27
    // For derivation of 4.820055089755540 refer to research diary entry: Monday 6th February 2007 (p26)
    CellProperties.nu_ie = ( ( (4.391070186767300e-27) / AVERAGE_PARTICLE_MASS ) * CellProperties.n[ELECTRON] * fLogLambda_ei( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] ) ) / pow( CellProperties.T[ELECTRON], 1.5 );
#endif // FORCE_SINGLE_FLUID
    // Calculated in EvaluateTerms when the collisional term is known
    CellProperties.collision_delta_t = Params.Duration;

// ******************************************************************************
// *    ADVECTION                                                               *
// ******************************************************************************

    // The time-step is calculated by the CFL condition
    CellProperties.advection_delta_t = SAFETY_ADVECTION * ( CellProperties.cell_width / ( fabs( CellProperties.v[1] ) + CellProperties.Cs ) );

// ******************************************************************************
// *    THERMAL CONDUCTION                                                      *
// ******************************************************************************

    // The time-step is calculated by: delta_t = ( n * k_B ) * ( cell width )^2 / ( 2.0 * coefficient of thermal conduction )
    // 0.5 * BOLTZMANN_CONSTANT = 6.9e-17
    term1 = SAFETY_CONDUCTION * (6.9e-17) * CellProperties.cell_width * CellProperties.cell_width;

#ifdef USE_KINETIC_MODEL
    CellProperties.conduction_delta_t[ELECTRON] = Params.Duration;
    for( j=1; j<SPECIES; j++ )
#else // USE_KINETIC_MODEL
    for( j=0; j<SPECIES; j++ )
#endif // USE_KINETIC_MODEL
    {
#ifdef OPTICALLY_THICK_RADIATION
        if( j == HYDROGEN && CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
            CellProperties.conduction_delta_t[j] = ( term1 * CellProperties.n[j] ) / ( ( Kappa[j] + pHI->Getkappa_0( log10(CellProperties.T[ELECTRON]) ) ) * pow( CellProperties.T[j], 2.5 ) );
        else
#endif // OPTICALLY_THICK_RADIATION
            CellProperties.conduction_delta_t[j] = ( term1 * CellProperties.n[j] ) / ( Kappa[j] * pow( CellProperties.T[j], 2.5 ) );
    }

// ******************************************************************************
// *    RADIATION                                                               *
// ******************************************************************************

    // Calculated in EvaluateTerms when the radiative term is known
    CellProperties.radiation_delta_t = Params.Duration;

// ******************************************************************************
// *    DYNAMIC VISCOSITY                                                       *
// ******************************************************************************

    // The time-step is calculated by: delta_t = ( rho ) * ( cell width )^2 / ( 2.0 * (4/3) * coefficient of dynamic viscosity )
    CellProperties.eta = DYNAMIC_VISCOSITY * ( pow( CellProperties.T[HYDROGEN], 2.5 ) / fLogLambda_ii( CellProperties.T[HYDROGEN], CellProperties.n[HYDROGEN] ) );
    CellProperties.viscosity_delta_t = SAFETY_VISCOSITY * ( ( CellProperties.rho[1] * CellProperties.cell_width * CellProperties.cell_width ) / ( 2.6666666666666666666666666666667 * CellProperties.eta ) );

    pActiveCell->UpdateCellProperties( &CellProperties );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}
}

void CEquations::EvaluateTerms( double current_time, double *delta_t, int iFirstStep )
{
PCELL pNextActiveCell, pFarLeftCell, pLeftCell, pRightCell;
CELLPROPERTIES CellProperties, FarLeftCellProperties, LeftCellProperties, RightCellProperties;

#ifdef NON_EQUILIBRIUM_RADIATION
PCELL pFarRightCell;
CELLPROPERTIES FarRightCellProperties;

// Variables used for time-dependent ionisation balance calculation
double **ppni0, **ppni1, **ppni2, **ppni3, **ppni4, ps[5], **ppdnibydt;
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
double fHI_c, frho_c, cell_width_cos_theta;
#endif // OPTICALLY_THICK_RADIATION

// Variables used by advective flux transport algorithm
double Q1, Q2, Q3, QT;

// Variables used by thermal and viscous flux transport algorithms
double T[3][SPECIES], Kappa[SPECIES], max_flux_coeff[SPECIES], Fc_max, v[2], n, P;
#ifdef NUMERICAL_VISCOSITY
double rho_v[2];
#endif // NUMERICAL_VISCOSITY

Kappa[ELECTRON] = SPITZER_ELECTRON_CONDUCTIVITY;
Kappa[HYDROGEN] = SPITZER_ION_CONDUCTIVITY;

max_flux_coeff[ELECTRON] = 1.5 / SQRT_ELECTRON_MASS; 
max_flux_coeff[HYDROGEN] = 1.5 / SQRT_AVERAGE_PARTICLE_MASS;

// Variables used for interpolation
double x[5], y[5], UpperValue, LowerValue, error;

// General variables
double term1, term2;
int j;

#ifdef USE_KINETIC_MODEL
CalculateKineticModel( iFirstStep );
#endif // USE_KINETIC_MODEL

#ifdef USE_TABULATED_CROSS_SECTION
double fCrossSection[3], fCellVolume;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *                                                                            *
// *    CALCULATE THE CELL-INTERFACE TERMS                                      *
// *                                                                            *
// ******************************************************************************

pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );

while( pNextActiveCell->pGetPointer( RIGHT ) )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    pLeftCell = pActiveCell->pGetPointer( LEFT );
    pLeftCell->GetCellProperties( &LeftCellProperties );

    pRightCell = pActiveCell->pGetPointer( RIGHT );
    pRightCell->GetCellProperties( &RightCellProperties );

    pFarLeftCell = pLeftCell->pGetPointer( LEFT );
    pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

// ******************************************************************************
// *    ADVECTIVE FLUX TRANSPORT ALGORITHM                                      *
// ******************************************************************************

    x[1] = LeftCellProperties.s[1];
    x[2] = CellProperties.s[1];
    y[1] = LeftCellProperties.v[1];
    y[2] = CellProperties.v[1];
    LinearFit( x, y, CellProperties.s[0], &(CellProperties.v[0]) );
    LeftCellProperties.v[2] = CellProperties.v[0];

    if( CellProperties.v[0] > 0.0 )
    {
// CALCULATE THE DENSITY
	 
        x[1] = FarLeftCellProperties.s[1];
	x[2] = LeftCellProperties.s[1];
	y[1] = FarLeftCellProperties.rho[1];
	y[2] = LeftCellProperties.rho[1];
	LinearFit( x, y, CellProperties.s[0], &Q1 );

	x[1] = LeftCellProperties.s[1];
	x[2] = CellProperties.s[1];
	y[1] = LeftCellProperties.rho[1];
	y[2] = CellProperties.rho[1];
	LinearFit( x, y, CellProperties.s[0], &Q2 );

        Q3 = LeftCellProperties.rho[1];

	if( CellProperties.rho[1] <= LeftCellProperties.rho[1] )
	{
	    QT = max( Q1, Q2 );
	    if( Q3 < QT )
	        CellProperties.rho[0] = Q3;
	    else
	        CellProperties.rho[0] = QT;
	}
	else
	{
	    QT = min( Q1, Q2 );
	    if( Q3 > QT )
	        CellProperties.rho[0] = Q3;
	    else
	        CellProperties.rho[0] = QT;
	}
	        
	LeftCellProperties.rho[2] = CellProperties.rho[0];
		
// CALCULATE THE MOMENTUM

        x[1] = FarLeftCellProperties.s[1];
	x[2] = LeftCellProperties.s[1];
	y[1] = FarLeftCellProperties.rho_v[1];
	y[2] = LeftCellProperties.rho_v[1];
	LinearFit( x, y, CellProperties.s[0], &Q1 );

	x[1] = LeftCellProperties.s[1];
	x[2] = CellProperties.s[1];
	y[1] = LeftCellProperties.rho_v[1];
	y[2] = CellProperties.rho_v[1];
	LinearFit( x, y, CellProperties.s[0], &Q2 );

        Q3 = LeftCellProperties.rho_v[1];

	if( CellProperties.rho_v[1] <= LeftCellProperties.rho_v[1] )
	{
	    QT = max( Q1, Q2 );
	    if( Q3 < QT )
	        CellProperties.rho_v[0] = Q3;
	    else
	        CellProperties.rho_v[0] = QT;
	}
	else
	{
	    QT = min( Q1, Q2 );
	    if( Q3 > QT )
	        CellProperties.rho_v[0] = Q3;
	    else
	        CellProperties.rho_v[0] = QT;
	}
	        
	LeftCellProperties.rho_v[2] = CellProperties.rho_v[0];

// CALCULATE THE ENERGY
		
	for( j=0; j<SPECIES; j++ )
	{
            x[1] = FarLeftCellProperties.s[1];
            x[2] = LeftCellProperties.s[1];
            y[1] = FarLeftCellProperties.TE_KE_P[1][j];
            y[2] = LeftCellProperties.TE_KE_P[1][j];
            LinearFit( x, y, CellProperties.s[0], &Q1 );

            x[1] = LeftCellProperties.s[1];
            x[2] = CellProperties.s[1];
            y[1] = LeftCellProperties.TE_KE_P[1][j];
            y[2] = CellProperties.TE_KE_P[1][j];
            LinearFit( x, y, CellProperties.s[0], &Q2 );

            Q3 = LeftCellProperties.TE_KE_P[1][j];

            if( CellProperties.TE_KE_P[1][j] <= LeftCellProperties.TE_KE_P[1][j] )
            {
                QT = max( Q1, Q2 );
		if( Q3 < QT )
                    CellProperties.TE_KE_P[0][j] = Q3;
		else
		    CellProperties.TE_KE_P[0][j] = QT;
            }
            else
            {
		QT = min( Q1, Q2 );
		if( Q3 > QT )
		    CellProperties.TE_KE_P[0][j] = Q3;
		else
		    CellProperties.TE_KE_P[0][j] = QT;
            }
	        
            LeftCellProperties.TE_KE_P[2][j] = CellProperties.TE_KE_P[0][j];
	}
    }
    else
    {
// CALCULATE THE DENSITY
	        
        x[1] = CellProperties.s[1];
	x[2] = RightCellProperties.s[1];
	y[1] = CellProperties.rho[1];
	y[2] = RightCellProperties.rho[1];
	LinearFit( x, y, CellProperties.s[0], &Q1 );

	x[1] = LeftCellProperties.s[1];
	x[2] = CellProperties.s[1];
	y[1] = LeftCellProperties.rho[1];
	y[2] = CellProperties.rho[1];
	LinearFit( x, y, CellProperties.s[0], &Q2 );

        Q3 = CellProperties.rho[1];

	if( CellProperties.rho[1] <= LeftCellProperties.rho[1] )
	{
            QT = min( Q1, Q2 );
	    if( Q3 > QT )
	        CellProperties.rho[0] = Q3;
	    else
	        CellProperties.rho[0] = QT;
	}
	else
	{
	    QT = max( Q1, Q2 );
	    if( Q3 < QT )
	        CellProperties.rho[0] = Q3;
	    else
	        CellProperties.rho[0] = QT;
	}
	        
	LeftCellProperties.rho[2] = CellProperties.rho[0];
		
// CALCULATE THE MOMENTUM
		
	x[1] = CellProperties.s[1];
	x[2] = RightCellProperties.s[1];
	y[1] = CellProperties.rho_v[1];
	y[2] = RightCellProperties.rho_v[1];
	LinearFit( x, y, CellProperties.s[0], &Q1 );

	x[1] = LeftCellProperties.s[1];
	x[2] = CellProperties.s[1];
	y[1] = LeftCellProperties.rho_v[1];
	y[2] = CellProperties.rho_v[1];
	LinearFit( x, y, CellProperties.s[0], &Q2 );

        Q3 = CellProperties.rho_v[1];

	if( CellProperties.rho_v[1] <= LeftCellProperties.rho_v[1] )
	{
	    QT = min( Q1, Q2 );
	    if( Q3 > QT )
	        CellProperties.rho_v[0] = Q3;
	    else
	        CellProperties.rho_v[0] = QT;
	}
	else
	{
	    QT = max( Q1, Q2 );
	    if( Q3 < QT )
	        CellProperties.rho_v[0] = Q3;
	    else
	        CellProperties.rho_v[0] = QT;
	}
	        
	LeftCellProperties.rho_v[2] = CellProperties.rho_v[0];

// CALCULATE THE ENERGY
		
        for( j=0; j<SPECIES; j++ )
	{
            x[1] = CellProperties.s[1];
            x[2] = RightCellProperties.s[1];
            y[1] = CellProperties.TE_KE_P[1][j];
            y[2] = RightCellProperties.TE_KE_P[1][j];
            LinearFit( x, y, CellProperties.s[0], &Q1 );

            x[1] = LeftCellProperties.s[1];
            x[2] = CellProperties.s[1];
            y[1] = LeftCellProperties.TE_KE_P[1][j];
            y[2] = CellProperties.TE_KE_P[1][j];
            LinearFit( x, y, CellProperties.s[0], &Q2 );

            Q3 = CellProperties.TE_KE_P[1][j];

            if( CellProperties.TE_KE_P[1][j] <= LeftCellProperties.TE_KE_P[1][j] )
            {
		QT = min( Q1, Q2 );
		if( Q3 > QT )
		    CellProperties.TE_KE_P[0][j] = Q3;
		else
                    CellProperties.TE_KE_P[0][j] = QT;
            }
            else
            {
		QT = max( Q1, Q2 );
		if( Q3 < QT )
		    CellProperties.TE_KE_P[0][j] = Q3;
		else
	        CellProperties.TE_KE_P[0][j] = QT;
            }
	        
            LeftCellProperties.TE_KE_P[2][j] = CellProperties.TE_KE_P[0][j];
	}
    }

// ******************************************************************************
// *    THERMAL FLUX TRANSPORT ALGORITHM                                        *
// ******************************************************************************

    x[1] = FarLeftCellProperties.s[1];
    x[2] = LeftCellProperties.s[1];
    x[3] = CellProperties.s[1];
    x[4] = RightCellProperties.s[1];

    // Only need to calculate the diffusive terms on the first time-step
    if( iFirstStep )
    {
#ifdef USE_KINETIC_MODEL
	j = ELECTRON;

	y[1] = FarLeftCellProperties.Fc[1][j];
	y[2] = LeftCellProperties.Fc[1][j];
	y[3] = CellProperties.Fc[1][j];
	y[4] = RightCellProperties.Fc[1][j];
		
	FitPolynomial4( x, y, CellProperties.s[0], &(CellProperties.Fc[0][j]), &error );
	LeftCellProperties.Fc[2][j] = CellProperties.Fc[0][j];

	for( j=1; j<SPECIES; j++ )
#else // USE_KINETIC_MODEL
	// Calculate the conducted heat fluxes for the different species
	for( j=0; j<SPECIES; j++ )
#endif // USE_KINETIC_MODEL
	{
            y[1] = FarLeftCellProperties.T[j];
            y[2] = LeftCellProperties.T[j];
            y[3] = CellProperties.T[j];
            y[4] = RightCellProperties.T[j];

            FitPolynomial4( x, y, ( CellProperties.s[1] - CellProperties.cell_width ), &(T[0][j]), &error );
            FitPolynomial4( x, y, CellProperties.s[0], &(T[1][j]), &error );
            T[2][j] = CellProperties.T[j];

            // Calculate the conducted heat flux at the left boundary
#ifdef OPTICALLY_THICK_RADIATION
            if( j == HYDROGEN && T[1][ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
                CellProperties.Fc[0][j] = - ( Kappa[j] + pHI->Getkappa_0( log10(T[1][ELECTRON]) ) ) * pow( T[1][j], 2.5 ) * ( ( T[2][j] - T[0][j] ) / CellProperties.cell_width );
            else
#endif // OPTICALLY_THICK_RADIATION
                CellProperties.Fc[0][j] = - Kappa[j] * pow( T[1][j], 2.5 ) * ( ( T[2][j] - T[0][j] ) / CellProperties.cell_width );

            // Estimate the maximum conducted heat flux (treats n as approximately constant across cell)
            // BOLTZMANN_CONSTANT^1.5 = 1.621132937e-24
            Fc_max = (1.621132937e-24) * max_flux_coeff[j] * CellProperties.n[j] * pow( T[1][j], 1.5 );

            term1 = CellProperties.Fc[0][j] * Fc_max;
            term2 = ( CellProperties.Fc[0][j] * CellProperties.Fc[0][j] ) + ( Fc_max * Fc_max );
            CellProperties.Fc[0][j] =  term1 / sqrt( term2 );

            LeftCellProperties.Fc[2][j] = CellProperties.Fc[0][j];
		
            if( pLeftCell->pGetPointer( LEFT )->pGetPointer( LEFT ) )
                LeftCellProperties.Fc[1][j] = 0.5 * ( LeftCellProperties.Fc[0][j] + LeftCellProperties.Fc[2][j] );
	}

// ******************************************************************************
// *    VISCOUS FLUX TRANSPORT ALGORITHM                                        *
// ******************************************************************************

#ifdef USE_KINETIC_MODEL
	j = ELECTRON;
	y[1] = FarLeftCellProperties.T[j];
	y[2] = LeftCellProperties.T[j];
	y[3] = CellProperties.T[j];
	y[4] = RightCellProperties.T[j];
	FitPolynomial4( x, y, CellProperties.s[0], &(T[1][j]), &error );
#endif // USE_KINETIC_MODEL

	j = HYDROGEN;
	y[1] = FarLeftCellProperties.n[j];
	y[2] = LeftCellProperties.n[j];
	y[3] = CellProperties.n[j];
	y[4] = RightCellProperties.n[j];
	FitPolynomial4( x, y, CellProperties.s[0], &n, &error );

	y[1] = FarLeftCellProperties.v[1];
	y[2] = LeftCellProperties.v[1];
	y[3] = CellProperties.v[1];
	y[4] = RightCellProperties.v[1];
	FitPolynomial4( x, y, ( CellProperties.s[1] - CellProperties.cell_width ), &(v[0]), &error );

	v[1] = CellProperties.v[1];

	// Calculate the viscosity flux at the left boundary
	term1 = DYNAMIC_VISCOSITY * ( pow( T[1][j], 2.5 ) / fLogLambda_ii( T[1][j], n ) );
	CellProperties.Feta[0] = term1 * 1.3333333333333333333333333333333 * ( ( v[1] - v[0] ) / CellProperties.cell_width );

	// The viscous flux forms the anisotropic part of the pressure tensor and therefore cannot exceed the total pressure
	// Hence, it must be limited to the value of the total pressure
	P = BOLTZMANN_CONSTANT * n * T[1][j];
	term1 = CellProperties.Feta[0] * P;
	term2 = ( CellProperties.Feta[0] * CellProperties.Feta[0] ) + ( P * P );
	CellProperties.Feta[0] =  term1 / sqrt( term2 );

	LeftCellProperties.Feta[2] = CellProperties.Feta[0];

	if( pLeftCell->pGetPointer( LEFT )->pGetPointer( LEFT ) )
	{
            LeftCellProperties.Feta[1] = 0.5 * ( LeftCellProperties.Feta[0] + LeftCellProperties.Feta[2] );
            LeftCellProperties.dvbyds = LeftCellProperties.Feta[1] / LeftCellProperties.eta;
	}

// ******************************************************************************
// *    NUMERICAL VISCOSITY                                                     *
// ******************************************************************************

#ifdef NUMERICAL_VISCOSITY
	y[1] = FarLeftCellProperties.rho_v[1];
	y[2] = LeftCellProperties.rho_v[1];
	y[3] = CellProperties.rho_v[1];
	y[4] = RightCellProperties.rho_v[1];
	FitPolynomial4( x, y, ( CellProperties.s[1] - CellProperties.cell_width ), &(rho_v[0]), &error );
	rho_v[1] = CellProperties.rho_v[1];

	term1 = ( CellProperties.cell_width * CellProperties.cell_width ) / ( 2.0 * RELATIVE_VISCOUS_TIME_SCALE * ( CellProperties.advection_delta_t / SAFETY_ADVECTION ) );
	CellProperties.Fnumerical[0] = term1 * ( ( rho_v[1] - rho_v[0] ) / CellProperties.cell_width );

        LeftCellProperties.Fnumerical[2] = CellProperties.Fnumerical[0];
#endif // NUMERICAL_VISCOSITY
    }

// ******************************************************************************
// *    CELL BOUNDARY PRESSURES                                                 *
// ******************************************************************************

    for( j=0; j<SPECIES; j++ )
    {
	y[1] = FarLeftCellProperties.P[1][j];
	y[2] = LeftCellProperties.P[1][j];
	y[3] = CellProperties.P[1][j];
	y[4] = RightCellProperties.P[1][j];
        FitPolynomial4( x, y, CellProperties.s[0], &(CellProperties.P[0][j]), &error );
	LeftCellProperties.P[2][j] = CellProperties.P[0][j];
    }

    pActiveCell->UpdateCellProperties( &CellProperties );
    pLeftCell->UpdateCellProperties( &LeftCellProperties );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// ******************************************************************************
// *                                                                            *
// *    CALCULATE THE CELL-CENTERED TERMS                                       *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    COLUMN NUMBER AND MASS DENSITIES                                        *
// ******************************************************************************

#ifdef OPTICALLY_THICK_RADIATION
// Left-hand leg of the loop
fHI_c = 0.0;
frho_c = 0.0;
pNextActiveCell = pCentreOfCurrentRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    {
        cell_width_cos_theta = CellProperties.cell_width * fabs( cos( ( _PI_ * CellProperties.s[1] ) / Params.L ) );
        // Calculate the neutral hydrogen column number density
        fHI_c += ( CellProperties.n[HYDROGEN] - ( CellProperties.n[ELECTRON] / 1.000144 ) ) * cell_width_cos_theta;
        // Calculate the mass column density
        frho_c += CellProperties.rho[1] * cell_width_cos_theta;

        CellProperties.HI_c = fHI_c;
        CellProperties.rho_c = frho_c;

        pActiveCell->UpdateCellProperties( &CellProperties );
    }
    pNextActiveCell = pActiveCell->pGetPointer( LEFT );
}

// Right-hand leg of the loop
fHI_c = 0.0;
frho_c = 0.0;
pNextActiveCell = pCentreOfCurrentRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    {
        cell_width_cos_theta = CellProperties.cell_width * fabs( cos( ( _PI_ * CellProperties.s[1] ) / Params.L ) );
        // Calculate the neutral hydrogen column number density
        fHI_c += ( CellProperties.n[HYDROGEN] - ( CellProperties.n[ELECTRON] / 1.000144 ) ) * cell_width_cos_theta;
        // Calculate the mass column density
        frho_c += CellProperties.rho[1] * cell_width_cos_theta;

        CellProperties.HI_c = fHI_c;
        CellProperties.rho_c = frho_c;

        pActiveCell->UpdateCellProperties( &CellProperties );
    }
    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}
#endif // OPTICALLY_THICK_RADIATION

// ******************************************************************************
// *    TERMS OF THE CONVERSATION EQUATIONS                                     *
// ******************************************************************************

pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );

while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    pLeftCell = pActiveCell->pGetPointer( LEFT );
    pLeftCell->GetCellProperties( &LeftCellProperties );

    pRightCell = pActiveCell->pGetPointer( RIGHT );
    pRightCell->GetCellProperties( &RightCellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
    pFarLeftCell = pLeftCell->pGetPointer( LEFT );
    pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

    pFarRightCell = pRightCell->pGetPointer( RIGHT );
    pFarRightCell->GetCellProperties( &FarRightCellProperties );
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef USE_TABULATED_CROSS_SECTION
    fCrossSection[0] = CalculateCrossSection( CellProperties.s[0] );
    fCrossSection[1] = CalculateCrossSection( CellProperties.s[1] );
    fCrossSection[2] = CalculateCrossSection( CellProperties.s[2] );
    fCellVolume = fCrossSection[1] * CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *                                                                            *
// *    MASS CONSERVATION                                                       *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    ADVECTION                                                               *
// ******************************************************************************

#ifdef USE_TABULATED_CROSS_SECTION
    LowerValue = CellProperties.rho[0] * CellProperties.v[0] * fCrossSection[0];
    UpperValue = CellProperties.rho[2] * CellProperties.v[2] * fCrossSection[2];
    CellProperties.rho_term[0] = - ( UpperValue - LowerValue ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
    LowerValue = CellProperties.rho[0] * CellProperties.v[0];
    UpperValue = CellProperties.rho[2] * CellProperties.v[2];
    CellProperties.rho_term[0] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

    CellProperties.drhobydt = CellProperties.rho_term[0];

// ******************************************************************************
// *                                                                            *
// *    MOMENTUM CONSERVATION                                                   *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    ADVECTION                                                               *
// ******************************************************************************

#ifdef USE_TABULATED_CROSS_SECTION
    LowerValue = CellProperties.rho_v[0] * CellProperties.v[0] * fCrossSection[0];
    UpperValue = CellProperties.rho_v[2] * CellProperties.v[2] * fCrossSection[2];
    CellProperties.rho_v_term[0] = - ( UpperValue - LowerValue ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
    LowerValue = CellProperties.rho_v[0] * CellProperties.v[0];
    UpperValue = CellProperties.rho_v[2] * CellProperties.v[2];
    CellProperties.rho_v_term[0] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *    PRESSURE GRADIENT                                                       *
// ******************************************************************************

    UpperValue = CellProperties.P[2][ELECTRON] + CellProperties.P[2][HYDROGEN];
    LowerValue = CellProperties.P[0][ELECTRON] + CellProperties.P[0][HYDROGEN];
    CellProperties.rho_v_term[1] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;

// ******************************************************************************
// *    GRAVITY                                                                 *
// ******************************************************************************

    CellProperties.rho_v_term[2] = CellProperties.rho[1] * CalculateGravity( CellProperties.s[1] );

    // Terms that must be integrated to first order in time only, otherwise they're unconditionally unstable
    if( iFirstStep )
    {
// ******************************************************************************
// *    VISCOUS STRESS                                                          *
// ******************************************************************************

#ifdef USE_TABULATED_CROSS_SECTION
	CellProperties.rho_v_term[3] = ( ( CellProperties.Feta[2] * fCrossSection[2] ) - ( CellProperties.Feta[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
	CellProperties.rho_v_term[3] = ( CellProperties.Feta[2] - CellProperties.Feta[0] ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *    NUMERICAL VISCOSITY                                                     *
// ******************************************************************************

#ifdef NUMERICAL_VISCOSITY
	// Numerical viscosity is used to stabilise the solutions as in the Lax scheme
#ifdef USE_TABULATED_CROSS_SECTION
	CellProperties.rho_v_term[4] = ( ( CellProperties.Fnumerical[2] * fCrossSection[2] ) - ( CellProperties.Fnumerical[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
	CellProperties.rho_v_term[4] = ( CellProperties.Fnumerical[2] - CellProperties.Fnumerical[0] ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION
#endif // NUMERICAL_VISCOSITY
    }

    CellProperties.drho_vbydt = CellProperties.rho_v_term[0] + CellProperties.rho_v_term[1] + CellProperties.rho_v_term[2] + CellProperties.rho_v_term[3] + CellProperties.rho_v_term[4];

// ******************************************************************************
// *                                                                            *
// *    ENERGY CONSERVATION                                                     *
// *                                                                            *
// ******************************************************************************

// ******************************************************************************
// *    ADVECTION                                                               *
// ******************************************************************************

#ifdef USE_TABULATED_CROSS_SECTION
    for( j=0; j<SPECIES; j++ )
    {
	LowerValue = CellProperties.TE_KE_P[0][j] * CellProperties.v[0] * fCrossSection[0];
	UpperValue = CellProperties.TE_KE_P[2][j] * CellProperties.v[2] * fCrossSection[2];
	CellProperties.TE_KE_term[0][j] = - ( UpperValue - LowerValue ) / fCellVolume;
    }
#else // USE_TABULATED_CROSS_SECTION
    for( j=0; j<SPECIES; j++ )
    {
	LowerValue = CellProperties.TE_KE_P[0][j] * CellProperties.v[0];
	UpperValue = CellProperties.TE_KE_P[2][j] * CellProperties.v[2];
	CellProperties.TE_KE_term[0][j] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
    }
#endif // USE_TABULATED_CROSS_SECTION

    // Terms that must be integrated to first order in time only, otherwise they're unconditionally unstable
    if( iFirstStep )
    {
// ******************************************************************************
// *    THERMAL CONDUCTION                                                      *
// ******************************************************************************

#ifdef USE_TABULATED_CROSS_SECTION
	for( j=0; j<SPECIES; j++ )
            CellProperties.TE_KE_term[1][j] = - ( ( CellProperties.Fc[2][j] * fCrossSection[2] ) - ( CellProperties.Fc[0][j] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
	for( j=0; j<SPECIES; j++ )
            CellProperties.TE_KE_term[1][j] = - ( CellProperties.Fc[2][j] - CellProperties.Fc[0][j] ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *    VISCOUS STRESS                                                          *
// ******************************************************************************

        // Heating due to the viscous stress and work done on (by) the flow by (on) the viscous stress
#ifdef USE_TABULATED_CROSS_SECTION
	CellProperties.TE_KE_term[7][HYDROGEN] = ( ( CellProperties.Feta[2] * CellProperties.v[2] * fCrossSection[2] ) - ( CellProperties.Feta[0] * CellProperties.v[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
	CellProperties.TE_KE_term[7][HYDROGEN] = ( ( CellProperties.Feta[2] * CellProperties.v[2] ) - ( CellProperties.Feta[0] * CellProperties.v[0] ) ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION

// ******************************************************************************
// *    NUMERICAL VISCOSITY                                                     *
// ******************************************************************************

#ifdef NUMERICAL_VISCOSITY
	// Numerical viscosity is used to stabilise the solutions as in the Lax scheme
#ifdef USE_TABULATED_CROSS_SECTION
	CellProperties.TE_KE_term[8][HYDROGEN] = ( ( CellProperties.Fnumerical[2] * CellProperties.v[2] * fCrossSection[2] ) - ( CellProperties.Fnumerical[0] * CellProperties.v[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_TABULATED_CROSS_SECTION
	CellProperties.TE_KE_term[8][HYDROGEN] = ( ( CellProperties.Fnumerical[2] * CellProperties.v[2] ) - ( CellProperties.Fnumerical[0] * CellProperties.v[0] ) ) / CellProperties.cell_width;
#endif // USE_TABULATED_CROSS_SECTION
#endif // NUMERICAL_VISCOSITY
    }

// ******************************************************************************
// *    GRAVITY                                                                 *
// ******************************************************************************

    CellProperties.TE_KE_term[2][HYDROGEN] = CellProperties.rho_v_term[2] * CellProperties.v[1];

// ******************************************************************************
// *    COLLISIONS                                                              *
// ******************************************************************************

    // The collision frequency is calculated using n_e, therefore the collisional coupling depends on (n_e)(n_H)
    CellProperties.TE_KE_term[3][ELECTRON] = (2.07e-16) * CellProperties.n[HYDROGEN] * CellProperties.nu_ie * ( CellProperties.T[HYDROGEN] - CellProperties.T[ELECTRON] );
    // Set the collisional timescale to 1% of the timescale for order of magnitude changes in the electron energy due to collisional energy exchange
    CellProperties.collision_delta_t = ( 0.01 * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[3][ELECTRON] );
    // If the collisional timescale is less than the minimum specified collisional timescale then scale the rate of energy exchange so that tiny timesteps can be avoided
    if( CellProperties.collision_delta_t < MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE )
        CellProperties.TE_KE_term[3][ELECTRON] *= CellProperties.collision_delta_t / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;

    CellProperties.TE_KE_term[3][HYDROGEN] = - CellProperties.TE_KE_term[3][ELECTRON];

// ******************************************************************************
// *    HEATING                                                                 *
// ******************************************************************************

    CellProperties.TE_KE_term[4][HEATED_SPECIES] = pHeat->CalculateHeating( CellProperties.s[1], current_time );
   
// ******************************************************************************
// *    RADIATION                                                               *
// ******************************************************************************

CellProperties.TE_KE_term[5][ELECTRON] = - SMALLEST_DOUBLE;
#ifdef OPTICALLY_THICK_RADIATION
    if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    {
        CellProperties.TE_KE_term[4][HEATED_SPECIES] += pHeat->CalculateVALHeating( log10( CellProperties.rho_c ) );
        CellProperties.TE_KE_term[5][ELECTRON] -= pHI->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10((4e-14)*CellProperties.HI_c), CellProperties.n[ELECTRON] * CellProperties.rho[1]);
        CellProperties.TE_KE_term[5][ELECTRON] -= pMgII->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10(CellProperties.rho_c), CellProperties.n[ELECTRON] * CellProperties.rho[1]);
        CellProperties.TE_KE_term[5][ELECTRON] -= pCaII->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10(CellProperties.rho_c), CellProperties.n[ELECTRON] * CellProperties.rho[1]);
        CellProperties.radiation_delta_t = ( SAFETY_RADIATION * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[5][ELECTRON] )
;
    }
    else
    {
	term1 = 1.0;
#else // OPTICALLY_THICK_RADIATION
    if( CellProperties.T[ELECTRON] < MINIMUM_RADIATION_TEMPERATURE )
    {
        // Provide some additional heating to the chromosphere if the temperature drops below the specified isothermal temperature
        CellProperties.TE_KE_term[5][ELECTRON] = ( ( ( MINIMUM_RADIATION_TEMPERATURE / CellProperties.T[ELECTRON] ) - 1.0 ) * CellProperties.TE_KE[1][ELECTRON] ) / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
	CellProperties.radiation_delta_t = MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
    }
    else
    {   
	term1 = 1.0;
	// Decrease the radiation smoothly to zero in the chromosphere
	if( CellProperties.T[ELECTRON] < lower_radiation_temperature_boundary )
	    term1 = ( CellProperties.T[ELECTRON] - MINIMUM_RADIATION_TEMPERATURE ) / ZERO_OVER_TEMPERATURE_INTERVAL;
#endif // OPTICALLY_THICK_RADIATION

#ifdef DECOUPLE_IONISATION_STATE_SOLVER
	#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
		printf( "nOPTICALLY_THICK_RADIATION; DECOUPLE_IONISATION_STATE_SOLVER; USE_POWER_LAW_RADIATIVE_LOSSES\n" );
		CellProperties.TE_KE_term[5][ELECTRON] -= term1 * pRadiation2->GetPowerLawRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) );
	#else // USE_POWER_LAW_RADIATIVE_LOSSES
		#ifdef NON_EQUILIBRIUM_RADIATION
			printf( "nOPTICALLY_THICK_RADIATION; DECOUPLE_IONISATION_STATE_SOLVER; nUSE_POWER_LAW_RADIATIVE_LOSSES; NON_EQUILIBRIUM_RADIATION\n" );
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) + pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) );
		#else // NON_EQUILIBRIUM_RADIATION
			printf( "nOPTICALLY_THICK_RADIATION; DECOUPLE_IONISATION_STATE_SOLVER; nUSE_POWER_LAW_RADIATIVE_LOSSES; nNON_EQUILIBRIUM_RADIATION\n" );
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) );
		#endif // NON_EQUILIBRIUM_RADIATION
	#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#else // DECOUPLE_IONISATION_STATE_SOLVER
	#ifdef NON_EQUILIBRIUM_RADIATION
		printf( "nOPTICALLY_THICK_RADIATION; nDECOUPLE_IONISATION_STATE_SOLVER; NON_EQUILIBRIUM_RADIATION\n" );
        	ppni2 = CellProperties.pIonFrac->ppGetIonFrac();
	        CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ), ppni2 ) + pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) );
	#else // NON_EQUILIBRIUM_RADIATION
		#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
			// printf( "nOPTICALLY_THICK_RADIATION; nDECOUPLE_IONISATION_STATE_SOLVER; nNON_EQUILIBRIUM_RADIATION; USE_POWER_LAW_RADIATIVE_LOSSES\n" );
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * pRadiation2->GetPowerLawRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) );
		#else // USE_POWER_LAW_RADIATIVE_LOSSES
			// printf( "nOPTICALLY_THICK_RADIATION; nDECOUPLE_IONISATION_STATE_SOLVER; nNON_EQUILIBRIUM_RADIATION; nUSE_POWER_LAW_RADIATIVE_LOSSES\n" );
        		CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ) ) );
		#endif // USE_POWER_LAW_RADIATIVE_LOSSES
	#endif // NON_EQUILIBRIUM_RADIATION
#endif // DECOUPLE_IONISATION_STATE_SOLVER
	    CellProperties.radiation_delta_t = ( SAFETY_RADIATION * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[5][ELECTRON] )
;
    }

// ******************************************************************************
// *    SMALL-SCALE ELECTRIC FIELDS                                             *
// ******************************************************************************

    // Derived from qnEv = v dP/ds
    // The term added to the electron energy equation is (e)nEv = v dPe/ds
    // The term added to the hydrogen energy equation is (-e)nEv = -v dPe/ds
    CellProperties.TE_KE_term[6][ELECTRON] = CellProperties.v[1] * ( ( CellProperties.P[2][ELECTRON] - CellProperties.P[0][ELECTRON] ) / CellProperties.cell_width );
    CellProperties.TE_KE_term[6][HYDROGEN] = -CellProperties.TE_KE_term[6][ELECTRON];

    for( j=0; j<SPECIES; j++ )
        CellProperties.dTE_KEbydt[j] = CellProperties.TE_KE_term[0][j] + CellProperties.TE_KE_term[1][j] + CellProperties.TE_KE_term[2][j] + CellProperties.TE_KE_term[3][j] + CellProperties.TE_KE_term[4][j] + CellProperties.TE_KE_term[5][j] + CellProperties.TE_KE_term[6][j] + CellProperties.TE_KE_term[7][j] + CellProperties.TE_KE_term[8][j];

// ******************************************************************************
// *                                                                            *
// *    TIME-DEPENDENT IONISATION                                               *
// *                                                                            *
// ******************************************************************************

#ifdef NON_EQUILIBRIUM_RADIATION
    ppni0 = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
    ppni1 = LeftCellProperties.pIonFrac->ppGetIonFrac();
    ppni2 = CellProperties.pIonFrac->ppGetIonFrac();
    ppni3 = RightCellProperties.pIonFrac->ppGetIonFrac();
    ppni4 = FarRightCellProperties.pIonFrac->ppGetIonFrac();

    ps[0] = FarLeftCellProperties.s[1];
    ps[1] = LeftCellProperties.s[1];
    ps[2] = CellProperties.s[1];
    ps[3] = RightCellProperties.s[1];
    ps[4] = FarRightCellProperties.s[1];

    ppdnibydt = CellProperties.pIonFrac->ppGetdnibydt();
    pRadiation->GetAlldnibydt( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ), ppni0, ppni1, ppni2, ppni3, ppni4, ps, CellProperties.s, CellProperties.v, CellProperties.cell_width, ppdnibydt, &(CellProperties.atomic_delta_t) );
#endif // NON_EQUILIBRIUM_RADIATION

    pActiveCell->UpdateCellProperties( &CellProperties );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// Find the smallest characteristic time-scale
GetSmallestTimeScale( delta_t, iFirstStep );
}

double CEquations::CalculateGravity( double s )
{
double x[3], y[3], g_parallel;
int i;

for( i=0; i<igdp; i++ )
{
    if( s < ppGravity[i][0] ) break;
}

if( i == 0 ) i = 1;
if( i == igdp ) i = igdp - 1;

x[1] = ppGravity[i-1][0];
x[2] = ppGravity[i][0];

y[1] = ppGravity[i-1][1];
y[2] = ppGravity[i][1];

LinearFit( x, y, s, &g_parallel );

return g_parallel;
}

#ifdef USE_TABULATED_CROSS_SECTION
double CEquations::CalculateCrossSection( double s )
{
double x[3], y[3], cross_section;
int i;

for( i=0; i<icsdp; i++ )
{
    if( s < ppCrossSection[i][0] ) break;
}

if( i == 0 ) i = 1;
if( i == icsdp ) i = icsdp - 1;

x[1] = ppCrossSection[i-1][0];
x[2] = ppCrossSection[i][0];

y[1] = ppCrossSection[i-1][1];
y[2] = ppCrossSection[i][1];

LinearFit( x, y, s, &cross_section );

return cross_section;
}
#endif // USE_TABULATED_CROSS_SECTION

void CEquations::GetSmallestTimeScale( double *delta_t, int iFirstStep )
{
if( !iFirstStep )
    return;
	
PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;
int j;

*delta_t = Params.Duration;

pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );

while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    // Advection timescale
    if( CellProperties.advection_delta_t < *delta_t )
	*delta_t = CellProperties.advection_delta_t;

    // Thermal conduction timescale
    for( j=0; j<SPECIES; j++ )
    {
	if( CellProperties.conduction_delta_t[j] < *delta_t )
            *delta_t = CellProperties.conduction_delta_t[j];
    }

    // Viscous timescale
    if( CellProperties.viscosity_delta_t < *delta_t )
	*delta_t = CellProperties.viscosity_delta_t;

    // Collisional timescale
    if( CellProperties.collision_delta_t >= MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE &&
        CellProperties.collision_delta_t < *delta_t )
        *delta_t = CellProperties.collision_delta_t;

    // Radiation timescale
    if( CellProperties.radiation_delta_t < *delta_t )
	*delta_t = CellProperties.radiation_delta_t;

#ifdef NON_EQUILIBRIUM_RADIATION
    // Ionisation and recombination timescale
    if( CellProperties.atomic_delta_t < *delta_t )
        *delta_t = CellProperties.atomic_delta_t;
#endif // NON_EQUILIBRIUM_RADIATION

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}
}

void CEquations::Half_Time_Step( PCELLPROPERTIES pNewCellProperties, double delta_t )
{
CELLPROPERTIES CellProperties;
int j;

pActiveCell->GetCellProperties( &CellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
// Integrate the ion fractions
pNewCellProperties->pIonFrac->IntegrateAllIonFrac( delta_t );
#endif // NON_EQUILIBRIUM_RADIATION

pNewCellProperties->rho[1] = CellProperties.rho[1] + ( delta_t * CellProperties.drhobydt );
pNewCellProperties->rho_v[1] = CellProperties.rho_v[1] + ( delta_t * CellProperties.drho_vbydt );

for( j=0; j<SPECIES; j++ )
    pNewCellProperties->TE_KE[1][j] = CellProperties.TE_KE[1][j] + ( delta_t * CellProperties.dTE_KEbydt[j] );
}

void CEquations::Full_Time_Step( PCELLPROPERTIES pCellProperties, double delta_t )
{
PCELL pBottomCell;
CELLPROPERTIES BottomCellProperties;
int j;

pBottomCell = pActiveCell->pGetPointer( BOTTOM );
pBottomCell->GetCellProperties( &BottomCellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
// Integrate the ion fractions using the previous set of values
pCellProperties->pIonFrac->CopyAllIonFrac( BottomCellProperties.pIonFrac );
pCellProperties->pIonFrac->IntegrateAllIonFrac( delta_t );
#endif // NON_EQUILIBRIUM_RADIATION

pCellProperties->rho[1] = BottomCellProperties.rho[1] + ( delta_t * pCellProperties->drhobydt );
pCellProperties->rho_v[1] = BottomCellProperties.rho_v[1] + ( delta_t * pCellProperties->drho_vbydt );

for( j=0; j<SPECIES; j++ )
    pCellProperties->TE_KE[1][j] = BottomCellProperties.TE_KE[1][j] + ( delta_t * pCellProperties->dTE_KEbydt[j] );
}

// ******************************************************************************
// *                                                                            *
// *    KINETIC COMPONENT FOR DISTRIBUTION FUNCTION CALCULATIONS                *
// *                                                                            *
// ******************************************************************************

#ifdef USE_KINETIC_MODEL
// Index labels for each column of the table in Spitzer & Harm, 1953, Phys. Rev., 89, 977
#define U							0		// Relative (to u_th) velocity values in column 0
#define X_E							1		// Coefficients when an electric field is present in column 1
#define X_T							2		// Coefficients when a temperature gradient is present in column 2

void CEquations::Get_SH_Table( void )
{
FILE *pFile;
int i;

pFile = fopen( "Kinetic_Model/config/spitzer_harm_table.txt" ,"r" );
for( i=0; i<51; i++ )
{
    ReadDouble( pFile, &(SH_Table[i][U]) );
    ReadDouble( pFile, &(SH_Table[i][X_E]) );
    ReadDouble( pFile, &(SH_Table[i][X_T]) );
}
fclose( pFile );
}

void CEquations::CountCells( void )
{
PCELL pNextActiveCell;

iNumCells = 0;

pNextActiveCell = pStartOfCurrentRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;

    iNumCells++;

    pNextActiveCell=pActiveCell->pGetPointer( RIGHT );
}
}

void CEquations::CreateIndexedCellList( void )
{
PCELL pNextActiveCell;
int iIndex;

// Free the previous indexed cell list and allocate sufficient memory for the next one
if( ppCellList )
    free( ppCellList );

ppCellList = (PCELL*)malloc( iNumCells * sizeof(PCELL) );

// Compile an indexed list of the cells
iIndex = 0;
pNextActiveCell = pStartOfCurrentRow;
while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;

    ppCellList[iIndex] = pActiveCell;
    iIndex++;
	
    pNextActiveCell=pActiveCell->pGetPointer( RIGHT );
}
}

void CEquations::CalculateKineticModel( int iFirstStep )
{
if( !iFirstStep )
    return;

PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;
double Tmax = 0.0;

// Calculate the physical quantities required for the kinetic model
pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    // Find the maximum temperature to determine the velocity range
    if( CellProperties.T[ELECTRON] > Tmax )
	Tmax = CellProperties.T[ELECTRON];

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// Calculate the quantities required for the kinetic model
pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    CellProperties.pKinetic->CalculateVelocityRange( CellProperties.T[ELECTRON], Tmax );
    CellProperties.pKinetic->CalculateCollisionalProperties( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );
    CellProperties.pKinetic->CalculateMaxDFN( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// Calculate the non-Maxwellian distribution function and its moments in each grid cell
CalculateNonMaxDFN();
}

void CEquations::CalculateNonMaxDFN( void )
{
// General variables
double *pupsilon, *pnu_ee, *pnu_ei, *plambda_ei, *pMaxDFN_ee, *pMaxDFN_ei, *pNonMaxDFN;
double dTbyds, K_SH, Kappa;
int iIndex, i, j, half_data_points = DISTRIBUTION_DATA_POINTS / 2;

// Variables for the BGK part of the solution
CELLPROPERTIES CellProperties;
double ds, previous_cell_width, *previous_pNonMaxDFN;
double term1, term2, term3;

// Variables for the Spitzer-Harm part of the solution
CELLPROPERTIES FarLeftCellProperties, LeftCellProperties, RightCellProperties, FarRightCellProperties;
double u_th, u_n;
double fScaleLength, E, lambda_ei, x_e, x_t, fp;
double x[6], y[6], fLowerValue, fUpperValue, fError;

// ******************************************************************************
// *    INITIALISE THE BOUNDARY CONDITIONS                                      *
// ******************************************************************************

iIndex = 0;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = 1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = iNumCells-2;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = iNumCells-1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

// ******************************************************************************
// *    BGK PART OF THE SOLUTION                                                *
// ******************************************************************************

// Solve for positive velocities

iIndex = 1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

for( iIndex=2; iIndex<iNumCells; iIndex++ )
{
    previous_cell_width = CellProperties.cell_width;
    previous_pNonMaxDFN = pNonMaxDFN;

    pActiveCell = ppCellList[iIndex];
    pActiveCell->GetCellProperties( &CellProperties );
    pupsilon = CellProperties.pKinetic->Get_pupsilon();
    pnu_ee = CellProperties.pKinetic->Get_pnu_ee();
    pnu_ei = CellProperties.pKinetic->Get_pnu_ei();
    pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
    pMaxDFN_ei = CellProperties.pKinetic->Get_pMaxDFN_ei();
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    ds = 0.5 * ( previous_cell_width + CellProperties.cell_width );

    for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
    {
	term1 = ( ( pnu_ee[i] * pMaxDFN_ee[i] ) + ( 4.0 * pnu_ei[i] * pMaxDFN_ei[i] ) ) * ds;
	term2 = pupsilon[i] * previous_pNonMaxDFN[i];
	term3 = pupsilon[i] + ( ( pnu_ee[i] + ( 4.0 * pnu_ei[i] ) ) * ds );

	pNonMaxDFN[i] = ( term1 + term2 ) / term3;
    }
}

// Solve for negative velocities

iIndex = iNumCells-2;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

for( iIndex=iNumCells-3; iIndex>=0; iIndex-- )
{
    previous_cell_width = CellProperties.cell_width;
    previous_pNonMaxDFN = pNonMaxDFN;

    pActiveCell = ppCellList[iIndex];
    pActiveCell->GetCellProperties( &CellProperties );
    pupsilon = CellProperties.pKinetic->Get_pupsilon();
    pnu_ee = CellProperties.pKinetic->Get_pnu_ee();
    pnu_ei = CellProperties.pKinetic->Get_pnu_ei();
    pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
    pMaxDFN_ei = CellProperties.pKinetic->Get_pMaxDFN_ei();
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    ds = 0.5 * ( previous_cell_width + CellProperties.cell_width );

    for( i=0; i<half_data_points; i++ )
    {
        term1 = ( ( pnu_ee[i] * pMaxDFN_ee[i] ) + ( 4.0 * pnu_ei[i] * pMaxDFN_ei[i] ) ) * ds;
	term2 = pupsilon[i] * previous_pNonMaxDFN[i];
	term3 = pupsilon[i] + ( ( pnu_ee[i] + ( 4.0 * pnu_ei[i] ) ) * ds );

	pNonMaxDFN[i] = ( term1 + term2 ) / term3;
    }
}

// ******************************************************************************
// *    SPITZER PART OF THE SOLUTION                                            *
// ******************************************************************************

for( iIndex=2; iIndex<iNumCells-2; iIndex++ )
{
    pActiveCell = ppCellList[iIndex];
    pActiveCell->GetCellProperties( &CellProperties );

    ppCellList[iIndex-2]->GetCellProperties( &FarLeftCellProperties );
    ppCellList[iIndex-1]->GetCellProperties( &LeftCellProperties );
    ppCellList[iIndex+1]->GetCellProperties( &RightCellProperties );
    ppCellList[iIndex+2]->GetCellProperties( &FarRightCellProperties );

    pupsilon = CellProperties.pKinetic->Get_pupsilon();
    plambda_ei = CellProperties.pKinetic->Get_plambda_ei();
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    u_th = CellProperties.pKinetic->Get_u_th();

    x[1] = FarLeftCellProperties.s[1];
    x[2] = LeftCellProperties.s[1];
    x[3] = CellProperties.s[1];
    x[4] = RightCellProperties.s[1];
    x[5] = FarRightCellProperties.s[1];

    y[1] = FarLeftCellProperties.T[ELECTRON];
    y[2] = LeftCellProperties.T[ELECTRON];
    y[3] = CellProperties.T[ELECTRON];
    y[4] = RightCellProperties.T[ELECTRON];
    y[5] = FarRightCellProperties.T[ELECTRON];

    FitPolynomial4( x, y, CellProperties.s[0], &fLowerValue, &fError );
    FitPolynomial4( &(x[1]), &(y[1]), CellProperties.s[2], &fUpperValue, &fError );

    dTbyds = ( fUpperValue - fLowerValue ) / CellProperties.cell_width;
    fScaleLength = KNUDSEN_NUMBER * fabs( CellProperties.T[ELECTRON] / dTbyds );
    E = -0.703 * ( BOLTZMANN_CONSTANT / ELECTRON_CHARGE ) * dTbyds;

    lambda_ei = fSHMeanFreePath( u_th, CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );

    // Solve for positive velocities

    for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
    {
        // Normalise the velocity to the thermal speed
	u_n = pupsilon[i] / u_th;

	if( plambda_ei[i] > fScaleLength || u_n > 3.20 )
            break;

	// Calculate X_E and X_T
        if( u_n > 0.10 )
	{
            // Find the data points to interpolate at the current velocity
            for( j=0; j<51; j++ )
            {
                if( u_n <= SH_Table[j][U] )
                    break;
            }

            x[1] = SH_Table[j-1][U];
            x[2] = SH_Table[j][U];

            y[1] = SH_Table[j-1][X_E];
            y[2] = SH_Table[j][X_E];
            LinearFit( x, y, u_n, &x_e );

            y[1] = SH_Table[j-1][X_T];
            y[2] = SH_Table[j][X_T];
            LinearFit( x, y, u_n, &x_t );
	}
	else
	{
            x_e = 0.0;
            x_t = 0.0;
	}

	term1 = -x_e * ( ( ELECTRON_CHARGE * E ) / ( BOLTZMANN_CONSTANT * CellProperties.T[ELECTRON] ) );
	term2 = 2.0 * x_t * ( 1.0 / CellProperties.T[ELECTRON] ) * dTbyds; 
	fp = lambda_ei * ( term1 + term2 );

	pNonMaxDFN[i] *= ( 1.0 + fp );
    }

    // Solve for negative velocities

    for( i=half_data_points-1; i>=0; i-- )
    {
        // Normalise the velocity to the thermal speed
	u_n = -pupsilon[i] / u_th;

	if( plambda_ei[i] > fScaleLength || u_n > 3.20 )
            break;

	// Calculate X_E and X_T
	if( u_n > 0.10 )
	{
            // Find the data points to interpolate at the current velocity
            for( j=0; j<51; j++ )
            {
                if( u_n <= SH_Table[j][U] )
                    break;
            }

            x[1] = SH_Table[j-1][U];
            x[2] = SH_Table[j][U];

            y[1] = SH_Table[j-1][X_E];
            y[2] = SH_Table[j][X_E];
            LinearFit( x, y, u_n, &x_e );

            y[1] = SH_Table[j-1][X_T];
            y[2] = SH_Table[j][X_T];
            LinearFit( x, y, u_n, &x_t );
	}
	else
	{
            x_e = 0.0;
            x_t = 0.0;
	}

	term1 = -x_e * ( ( ELECTRON_CHARGE * E ) / ( BOLTZMANN_CONSTANT * CellProperties.T[ELECTRON] ) );
	term2 = 2.0 * x_t * ( 1.0 / CellProperties.T[ELECTRON] ) * dTbyds; 
	fp = - lambda_ei * ( term1 + term2 );

	pNonMaxDFN[i] *= ( 1.0 + fp );
    }

    CellProperties.Fc[1][ELECTRON] = CellProperties.pKinetic->Get_Fc();

    // Calculate the thermal conduction time-scale
    K_SH = SPITZER_ELECTRON_CONDUCTIVITY * pow( CellProperties.T[ELECTRON], 2.5 );
    Kappa = fabs( CellProperties.Fc[1][ELECTRON] / dTbyds );
	
    if( Kappa > K_SH )
        Kappa = K_SH;

    // The time-step is calculated by: delta_t = ( n * k_B ) * ( cell width )^2 / ( 2.0 * coefficient of thermal conduction )
    // 0.5 * BOLTZMANN_CONSTANT = 6.9e-17
    term1 = SAFETY_CONDUCTION * (6.9e-17) * CellProperties.cell_width * CellProperties.cell_width;
    CellProperties.conduction_delta_t[ELECTRON] = ( term1 * CellProperties.n[ELECTRON] ) / Kappa;

    pActiveCell->UpdateCellProperties( &CellProperties );
}
}
#endif // USE_KINETIC_MODEL
