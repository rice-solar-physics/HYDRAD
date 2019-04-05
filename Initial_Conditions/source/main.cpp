// ****
// *
// * A hydrostatic code that calculates a set of initial conditions 
// * for the hydrodynamic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 04/05/2019
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "params.h"
#include "ode.h"
#include "misc.h"
#include "../../Radiation_Model/source/radiation.h"
#ifdef OPTICALLY_THICK_RADIATION
#include "../../Radiation_Model/source/OpticallyThick/OpticallyThickIon.h"
#ifdef NLTE_CHROMOSPHERE
#include <iostream>
#include <fstream>
#include <vector>
#include "../../Radiation_Model/source/OpticallyThick/RadiativeRates.h"
using namespace std;
#endif // NLTE_CHROMOSPHERE
#include  "../../Resources/source/fitpoly.h"
#endif // OPTICALLY_THICK_RADIATION
#include "../../Resources/source/file.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/Utils/regPoly/regpoly.h"


int main(void)
{
double FindHeatingRange( double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients );
double RefineSolution( double Log_10H0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients );
int CalculateSolution( double finalH0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients );
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, double *pfGravityCoefficients );
#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void RecalculateElectronDensity( PARAMETERS Params );
int AMR2PHY( PARAMETERS Params );
void RecalculateChromosphericHeating( PARAMETERS Params, int number_of_lines );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

PARAMETERS Params;
PRADIATION pRadiation;

FILE *pFile;
#ifdef USE_POLY_FIT_TO_GRAVITY
#else // USE_POLY_FIT_TO_GRAVITY
char szGravityFilename[512];
#endif // USE_POLY_FIT_TO_GRAVITY
double *pfGravityCoefficients;
int i;

double *s, *P, *T, *nH, *ne;
double Log_10H0, finalH0;

int iTRplusCoronaplusTRSteps, iTotalSteps;

printf( "\n\nCalculating initial hydrostatic conditions...\n\n" );

// Get the user-specified parameter values from the configuration file
GetConfigurationParameters( &Params );

// Initialise the radiative losses
pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_eq.cfg" );

// Initialise the gravitational geometry
#ifdef USE_POLY_FIT_TO_GRAVITY
pFile = fopen( POLY_FIT_TO_GRAVITY_FILE, "r" );
#else // USE_POLY_FIT_TO_GRAVITY
GenerateSemiCircularLoop( Params );
sprintf( szGravityFilename, "%s.gravity", Params.szOutputFilename );
pFile = fopen( szGravityFilename, "r" );
#endif // USE_POLY_FIT_TO_GRAVITY
pfGravityCoefficients = (double*)malloc( sizeof(double) * (POLY_ORDER+1) );
for( i=0; i<POLY_ORDER+1; i++ )
    ReadDouble( pFile, &(pfGravityCoefficients[i]) );
fclose( pFile );

// Allocate memory to store the hydrostatic profiles
s = (double*)malloc( sizeof(double) * MAX_CELLS );
P = (double*)malloc( sizeof(double) * MAX_CELLS );
T = (double*)malloc( sizeof(double) * MAX_CELLS );
nH = (double*)malloc( sizeof(double) * MAX_CELLS );
ne = (double*)malloc( sizeof(double) * MAX_CELLS );

#ifdef ISOTHERMAL
finalH0 = 0.0;
#else // ISOTHERMAL
Log_10H0 = FindHeatingRange( s, P, T, nH, ne, Params, pRadiation, pfGravityCoefficients );
finalH0 = RefineSolution( Log_10H0, s, P, T, nH, ne, Params, pRadiation, pfGravityCoefficients );
#endif // ISOTHERMAL
iTRplusCoronaplusTRSteps = CalculateSolution( finalH0, s, P, T, nH, ne, Params, pRadiation, pfGravityCoefficients );
iTotalSteps = AddChromospheres( iTRplusCoronaplusTRSteps, s, P, T, nH, ne, Params, pfGravityCoefficients );

printf( "Writing initial conditions file...\n\n" );

WriteAMRFile( iTotalSteps, s, T, nH, ne, Params );
WritePHYFile( iTotalSteps, s, T, nH, ne, Params );
WriteSOLFile( finalH0, Params );

// Free the memory allocated to the hydrostatic profiles
free( s );
free( P );
free( T );
free( nH );
free( ne );

// Free the memory allocated to the gravitational geometry
free( pfGravityCoefficients );

// Close the radiative losses
delete pRadiation;

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
// Recalculate the electron density using the model 6-level hydrogen atom and output a new .amr file
RecalculateElectronDensity( Params );
// Output a new .phy file
iTotalSteps = AMR2PHY( Params );
// Recalculate the chromospheric heating as a function of column mass because the optically-thick radiation will have changed due to the new electron density profile
RecalculateChromosphericHeating( Params, iTotalSteps );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

printf( "Done!\n\n" );

return 0;
}

double FindHeatingRange( double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, ne2, nH2;
double Log_10H0, H0, H, R;
double dPbyds, dFcbyds, dTbyds;
// **** NEW ****
double FracDiff_T, FracDiff_nH;
// **** NEW ****
int iStep;

max_ds = Params.Lfull / MIN_CELLS;

sL = Params.s0;
sR = Params.Lfull - sL;

// First find the lower- and upper-limits of the heating range
for( Log_10H0=Params.Log_10H0[0]; Log_10H0<=Params.Log_10H0[1]; Log_10H0+=Params.dLog_10H0 )
{
    iStep = 0;

    H0 = pow( 10.0, Log_10H0 );

    // Set the initial conditions
    Fc = 0.0;
    nH[iStep] = Params.n0;
#ifdef OPTICALLY_THICK_RADIATION
    // 1.000144 = 1.0 + 1.44e-4
    ne[iStep] = 1.000144 * nH[iStep];
#else // OPTICALLY_THICK_RADIATION
    ne[iStep] = nH[iStep];
#endif // OPTICALLY_THICK_RADIATION
    T[iStep] = Params.T0;
    P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

    s[iStep] = 0.0 + sL;
    ds = MIN_DS;

    for( ;; ) {
        do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

            // Calculate the heat input and the radiation
            H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), ne[iStep], nH[iStep] );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T[iStep] ), ne[iStep], nH[iStep] ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), ne[iStep], nH[iStep] ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc, T[iStep] );

            P2 = P[iStep] + ( dPbyds * (ds/2.0) );
            Fc2 = Fc + ( dFcbyds * (ds/2.0) );
            T2 = T[iStep] + ( dTbyds * (ds/2.0) );
            nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );
#ifdef OPTICALLY_THICK_RADIATION
    	    // 1.000144 = 1.0 + 1.44e-4
    	    ne2 = 1.000144 * nH2;
#else // OPTICALLY_THICK_RADIATION
            ne2 = nH2;
#endif // OPTICALLY_THICK_RADIATION

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, Params.Lfull, pfGravityCoefficients );

            // Calculate the heat input and the radiation
            H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T2 ), ne2, nH2 );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T2 ), ne2, nH2 ) + pRadiation->GetFreeFreeRad( log10( T2 ), ne2, nH2 ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************
// **** NEW ****
            T[iStep+1] = T[iStep] + ( dTbyds * ds );
            P[iStep+1] = P[iStep] + ( dPbyds * ds );
            nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );

            FracDiff_T = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
            FracDiff_nH = fabs( 1.0 - ( nH[iStep+1] / nH[iStep] ) );
            if( FracDiff_T > EPSILON || FracDiff_nH > EPSILON )
            {
		ds /= MAX_VARIATION;
		if( ds < MIN_DS )
		{
			ds = MIN_DS;
		}
            }

        } while ( ( FracDiff_T > EPSILON || FracDiff_nH > EPSILON ) && ds > MIN_DS );
// **** NEW ****
        s[iStep+1] = s[iStep] + ds;
        if( s[iStep+1] >= sR ) break;

        if( T[iStep+1] < Params.T0 ) break;

        Fc += dFcbyds * ds;
        P[iStep+1] = P[iStep] + ( dPbyds * ds );
        nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );
#ifdef OPTICALLY_THICK_RADIATION
        // 1.000144 = 1.0 + 1.44e-4
        ne[iStep+1] = 1.000144 * nH[iStep+1];
#else // OPTICALLY_THICK_RADIATION
        ne[iStep+1] = nH[iStep+1];
#endif // OPTICALLY_THICK_RADIATION

        ds *= MAX_VARIATION;
        if( ds > max_ds ) ds = max_ds;

        iStep++;
    }

    if( s[iStep+1] < sR ) break;
}

return Log_10H0;
}

double RefineSolution( double Log_10H0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, ne2, nH2;
double H0, H, R;
double dPbyds, dFcbyds, dTbyds;
// **** NEW ****
double FracDiff_T, FracDiff_nH;
// **** NEW ****
int iStep;

double H0lower, H0upper, dH0, finalH0;
double minFc = LARGEST_DOUBLE;

max_ds = Params.Lfull / MIN_CELLS;

sL = Params.s0;
sR = Params.Lfull - sL;

H0lower = pow( 10.0, (Log_10H0-Params.dLog_10H0) );
H0upper = pow( 10.0, (Log_10H0) );
dH0 = ( H0upper - H0lower ) / Params.Hintervals;

printf( "Peak heating range = %.4e -> %.4e erg cm^-3 s^-1\n\n", H0lower, H0upper );

finalH0 = H0upper;

// Now find the solution within the calculated heating range
for( H0=H0lower; H0<=H0upper; H0+=dH0 )
{
    iStep = 0;

    // Set the initial conditions
    Fc = 0.0;
    nH[iStep] = Params.n0;
#ifdef OPTICALLY_THICK_RADIATION
    // 1.000144 = 1.0 + 1.44e-4
    ne[iStep] = 1.000144 * nH[iStep];
#else // OPTICALLY_THICK_RADIATION
    ne[iStep] = nH[iStep];
#endif // OPTICALLY_THICK_RADIATION
    T[iStep] = Params.T0;
    P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

    s[iStep] = 0.0 + sL;
    ds = MIN_DS;

    for( ;; ) {
        do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients);

            // Calculate the heat input and the radiation
            H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), ne[iStep], nH[iStep] );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T[iStep] ), ne[iStep], nH[iStep] ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), ne[iStep], nH[iStep] ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc, T[iStep] );

            P2 = P[iStep] + ( dPbyds * (ds/2.0) );
            Fc2 = Fc + ( dFcbyds * (ds/2.0) );
            T2 = T[iStep] + ( dTbyds * (ds/2.0) );
            nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );
#ifdef OPTICALLY_THICK_RADIATION
    	    // 1.000144 = 1.0 + 1.44e-4
    	    ne2 = 1.000144 * nH2;
#else // OPTICALLY_THICK_RADIATION
            ne2 = nH2;
#endif // OPTICALLY_THICK_RADIATION

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, Params.Lfull, pfGravityCoefficients );

            // Calculate the heat input and the radiation
            H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T2 ), ne2, nH2 );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T2 ), ne2, nH2 ) + pRadiation->GetFreeFreeRad( log10( T2 ), ne2, nH2 ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************
// **** NEW ****
            T[iStep+1] = T[iStep] + ( dTbyds * ds );
            P[iStep+1] = P[iStep] + ( dPbyds * ds );
            nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );

            FracDiff_T = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
            FracDiff_nH = fabs( 1.0 - ( nH[iStep+1] / nH[iStep] ) );
            if( FracDiff_T > EPSILON || FracDiff_nH > EPSILON )
            {
		ds /= MAX_VARIATION;
		if( ds < MIN_DS )
		{
			ds = MIN_DS;
		}
            }

        } while ( ( FracDiff_T > EPSILON || FracDiff_nH > EPSILON ) && ds > MIN_DS );
// **** NEW ****
        s[iStep+1] = s[iStep] + ds;
        if( s[iStep+1] >= sR ) break;

        if( T[iStep+1] < Params.T0 ) break;

        Fc += dFcbyds * ds;
        P[iStep+1] = P[iStep] + ( dPbyds * ds );
        nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );
#ifdef OPTICALLY_THICK_RADIATION
        // 1.000144 = 1.0 + 1.44e-4
        ne[iStep+1] = 1.000144 * nH[iStep+1];
#else // OPTICALLY_THICK_RADIATION
        ne[iStep+1] = nH[iStep+1];
#endif // OPTICALLY_THICK_RADIATION

        ds *= MAX_VARIATION;
        if( ds > max_ds ) ds = max_ds;

        iStep++;
    }

    if( s[iStep+1] < sR ) break;

    if( ( s[iStep+1] >= sR ) && ( Fc > 0.0 ) && ( fabs(Fc) < fabs(minFc) ) )
    {
        minFc = Fc;
	finalH0 = H0;
    }
}

return finalH0;
}

int CalculateSolution( double finalH0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, double *pfGravityCoefficients )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, ne2, nH2;
double H0, H, R;
double dPbyds, dFcbyds, dTbyds;
// **** NEW ****
double FracDiff_T, FracDiff_nH;
// **** NEW ****
int iStep;

max_ds = Params.Lfull / MIN_CELLS;

sL = Params.s0;
sR = Params.Lfull - sL;

iStep = 0;

// Set the initial conditions
Fc = 0.0;
nH[iStep] = Params.n0;
#ifdef OPTICALLY_THICK_RADIATION
// 1.000144 = 1.0 + 1.44e-4
ne[iStep] = 1.000144 * nH[iStep];
#else // OPTICALLY_THICK_RADIATION
ne[iStep] = nH[iStep];
#endif // OPTICALLY_THICK_RADIATION
T[iStep] = Params.T0;
P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

s[iStep] = 0.0 + sL;
ds = MIN_DS;

H0 = finalH0;
printf( "Optimum peak heating rate = %.4e erg cm^-3 s^-1\n\n", H0 );

while( s[iStep] <= sR ) {
    do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

        // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
        R = 0.0;
#else // ISOTHERMAL
        H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
        R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), ne[iStep], nH[iStep] );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
        R = - ( pRadiation->GetRadiation( log10( T[iStep] ), ne[iStep], nH[iStep] ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), ne[iStep], nH[iStep] ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#endif // ISOTHERMAL

        // Get the heat flux gradient
        dFcbyds = CalcdFcbyds( H, R );

        // Get the temperature gradient
        dTbyds = CalcdTbyds( Fc, T[iStep] );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        Fc2 = Fc + ( dFcbyds * (ds/2.0) );
        T2 = T[iStep] + ( dTbyds * (ds/2.0) );
        nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );
#ifdef OPTICALLY_THICK_RADIATION
    	// 1.000144 = 1.0 + 1.44e-4
    	ne2 = 1.000144 * nH2;
#else // OPTICALLY_THICK_RADIATION
        ne2 = nH2;
#endif // OPTICALLY_THICK_RADIATION

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, Params.Lfull, pfGravityCoefficients );

        // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
        R = 0.0;
#else // ISOTHERMAL
        H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
        R = - pRadiation->GetPowerLawRad( log10( T2 ), ne2, nH2 );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
        R = - ( pRadiation->GetRadiation( log10( T2 ), ne2, nH2 ) + pRadiation->GetFreeFreeRad( log10( T2 ), ne2, nH2 ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#endif // ISOTHERMAL

        // Get the heat flux gradient
        dFcbyds = CalcdFcbyds( H, R );

        // Get the temperature gradient
        dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************
// **** NEW ****
            T[iStep+1] = T[iStep] + ( dTbyds * ds );
            P[iStep+1] = P[iStep] + ( dPbyds * ds );
            nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );

            FracDiff_T = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
            FracDiff_nH = fabs( 1.0 - ( nH[iStep+1] / nH[iStep] ) );
            if( FracDiff_T > EPSILON || FracDiff_nH > EPSILON )
            {
		ds /= MAX_VARIATION;
		if( ds < MIN_DS )
		{
			ds = MIN_DS;
		}
            }

        } while ( ( FracDiff_T > EPSILON || FracDiff_nH > EPSILON ) && ds > MIN_DS );
// **** NEW ****
    s[iStep+1] = s[iStep] + ds;

    Fc += dFcbyds * ds;
    P[iStep+1] = P[iStep] + ( dPbyds * ds );
    nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );
#ifdef OPTICALLY_THICK_RADIATION
    // 1.000144 = 1.0 + 1.44e-4
    ne[iStep+1] = 1.000144 * nH[iStep+1];
#else // OPTICALLY_THICK_RADIATION
    ne[iStep+1] = nH[iStep+1];
#endif // OPTICALLY_THICK_RADIATION

    ds *= MAX_VARIATION;
    if( ds > max_ds ) ds = max_ds;

    iStep++;
}

return iStep;
}

#ifdef OPTICALLY_THICK_RADIATION
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, double *pfGravityCoefficients )
{
double GetVALTemperature( double s, int iVALTemperatureDP, double **ppVALTemperature );

double ds, max_ds;
double P2, T2, nT, ne2, nH2;
double dPbyds;
double FracDiff;
int i, iStep, iCHRSteps, iTotalSteps;

// Get the temperature structure of the VAL atmosphere
FILE *pFile;
double **ppVALTemperature;
int iVALTemperatureDP;
pFile = fopen ( "Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.T", "r" );
// Get the number of data points in the file
fscanf( pFile, "%i", &iVALTemperatureDP );
// Allocate sufficient memory to hold the VAL atmosphere data
ppVALTemperature = (double**)malloc( sizeof(double) * 2 );
for( i=0; i<2; i++ )
    ppVALTemperature[i] = (double*)malloc( sizeof(double) * iVALTemperatureDP );
for( i=0; i<iVALTemperatureDP; i++ )
{
    // Array index [0][i] contain the spatial coordinate and [1][i] contain the volumetric heating rate
    ReadDouble( pFile, &(ppVALTemperature[0][i]) );
    ReadDouble( pFile, &(ppVALTemperature[1][i]) );
}
fclose( pFile );

// Create the optically-thick ion objects
POPTICALLYTHICKION pHI;
pHI = new COpticallyThickIon( 1, (char *)"h_1", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );

max_ds = Params.Lfull / MIN_CELLS;

// Left-hand chromosphere

// Shift the TR + corona + TR solution in the array to make room for the chromosphere
iStep = ( MAX_CELLS - iTRplusCoronaplusTRSteps ) / 2;
for( i=iTRplusCoronaplusTRSteps-1; i>=0; i-- )
{
    s[i+iStep] = s[i];
    P[i+iStep] = P[i];
    T[i+iStep] = T[i];
    nH[i+iStep] = nH[i];
    ne[i+iStep] = ne[i];
}

iCHRSteps = 0;

nH[iStep] = Params.n0;
// 1.000144 = 1.0 + 1.44e-4
ne[iStep] = 1.000144 * nH[iStep];
T[iStep] = Params.T0;
P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

ds = - MIN_DS;
s[iStep] = s[iStep+1] + ds;

for( ;; ) {
    do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = GetVALTemperature( s[iStep]+(ds/2.0), iVALTemperatureDP, ppVALTemperature );
        nT = P2 / ( BOLTZMANN_CONSTANT * T2 );
        nH2 = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T2) ) ) ) ) * nT;
        ne2 = nT - ( ( 1.0 - 1.44e-4 ) * nH2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep]+(ds/2.0), nH2, Params.Lfull, pfGravityCoefficients );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

        P[iStep-1] = P[iStep] + ( dPbyds * ds );
        FracDiff = fabs( 1.0 - ( P[iStep] / P[iStep-1] ) );
        if( FracDiff > EPSILON )
        {
            ds /= MAX_VARIATION;
            if( -ds < MIN_DS )
                ds = - MIN_DS;
        }

    } while ( FracDiff > EPSILON && -ds > MIN_DS );

    s[iStep-1] = s[iStep] + ds;
    if( s[iStep-1] < 0.0 ) break;

    T[iStep-1] = GetVALTemperature( s[iStep-1], iVALTemperatureDP, ppVALTemperature );
    nT = P[iStep-1] / ( BOLTZMANN_CONSTANT * T[iStep-1] );
    nH[iStep-1] = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T[iStep-1]) ) ) ) ) * nT;
    ne[iStep-1] = nT - ( ( 1.0 - 1.44e-4 ) * nH[iStep-1] );

    ds *= MAX_VARIATION;
    if( ds < -max_ds ) ds = - max_ds;

    iStep--;
    iCHRSteps++;
}

// Shift the chromosphere + TR + corona + TR solution in the array back to element zero
for( i=iStep; i<iStep+iCHRSteps+iTRplusCoronaplusTRSteps-1; i++ )
{
    s[i-iStep] = s[i];
    P[i-iStep] = P[i];
    T[i-iStep] = T[i];
    nH[i-iStep] = nH[i];
    ne[i-iStep] = ne[i];
}

// Right-hand chromosphere

iStep = iTRplusCoronaplusTRSteps + iCHRSteps - 1;

nH[iStep] = Params.n0;
// 1.000144 = 1.0 + 1.44e-4
ne[iStep] = 1.000144 * nH[iStep];
T[iStep] = Params.T0;
P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

ds = s[iStep-1] - s[iStep-2];
s[iStep] = s[iStep-1] + ds;

for( ;; ) {
    do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = GetVALTemperature( Params.Lfull - ( s[iStep]+(ds/2.0) ), iVALTemperatureDP, ppVALTemperature );
        nT = P2 / ( BOLTZMANN_CONSTANT * T2 );
        nH2 = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T2) ) ) ) ) * nT;
        ne2 = nT - ( ( 1.0 - 1.44e-4 ) * nH2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep]+(ds/2.0), nH2, Params.Lfull, pfGravityCoefficients );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

        P[iStep+1] = P[iStep] + ( dPbyds * ds );
        FracDiff = fabs( 1.0 - ( P[iStep+1] / P[iStep] ) );
        if( FracDiff > EPSILON )
        {
            ds /= MAX_VARIATION;
            if( ds < MIN_DS )
		ds = MIN_DS;
        }

    } while ( FracDiff > EPSILON && ds > MIN_DS );

    s[iStep+1] = s[iStep] + ds;
    if( s[iStep+1] > Params.Lfull ) break;

    T[iStep+1] = GetVALTemperature( Params.Lfull - s[iStep+1], iVALTemperatureDP, ppVALTemperature );
    nT = P[iStep+1] / ( BOLTZMANN_CONSTANT * T[iStep+1] );
    nH[iStep+1] = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T[iStep+1]) ) ) ) ) * nT;
    ne[iStep+1] = nT - ( ( 1.0 - 1.44e-4 ) * nH[iStep+1] );

    ds *= MAX_VARIATION;
    if( ds > max_ds ) ds = max_ds;

    iStep++;
}

iTotalSteps = iStep + 1;

// Delete the optically-thick ion object
delete pHI;

// Free the memory allocated to the VAL atmosphere
for( i=0; i<2; i++ )
    free( ppVALTemperature[i] );
free( ppVALTemperature );

return iTotalSteps;
}

double GetVALTemperature( double s, int iVALTemperatureDP, double **ppVALTemperature )
{
int i;
double x[5], y[5], fVALTemperature, error;

// Trap limits
if( s < ppVALTemperature[0][0] ) return ppVALTemperature[1][0];
else if( s > ppVALTemperature[0][iVALTemperatureDP-1] ) return ppVALTemperature[1][iVALTemperatureDP-1];

for( i=0; i<iVALTemperatureDP; i++ )
{
    if( s <= ppVALTemperature[0][i] ) break;
}

if( i < 2 ) i = 2;
else if( i > iVALTemperatureDP - 2 ) i = iVALTemperatureDP - 2;

x[1] = ppVALTemperature[0][i-2];
x[2] = ppVALTemperature[0][i-1];
x[3] = ppVALTemperature[0][i];
x[4] = ppVALTemperature[0][i+1];
y[1] = ppVALTemperature[1][i-2];
y[2] = ppVALTemperature[1][i-1];
y[3] = ppVALTemperature[1][i];
y[4] = ppVALTemperature[1][i+1];

FitPolynomial4( x, y, s, &fVALTemperature, &error );

return fVALTemperature;
}
#else // OPTICALLY_THICK_RADIATION
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, double *pfGravityCoefficients )
{
double ds, max_ds;
double P2, T2, nH2;
double dPbyds;
double FracDiff;
int i, iStep, iCHRSteps, iTotalSteps;

max_ds = Params.Lfull / MIN_CELLS;

// Left-hand chromosphere

// Shift the TR + corona + TR solution in the array to make room for the chromosphere
iStep = ( MAX_CELLS - iTRplusCoronaplusTRSteps ) / 2;
for( i=iTRplusCoronaplusTRSteps-1; i>=0; i-- )
{
    s[i+iStep] = s[i];
    P[i+iStep] = P[i];
    T[i+iStep] = T[i];
    nH[i+iStep] = nH[i];
    ne[i+iStep] = ne[i];
}

iCHRSteps = 0;

nH[iStep] = Params.n0;
ne[iStep] = nH[iStep];
T[iStep] = Params.T0;
P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

ds = - MIN_DS;
s[iStep] = s[iStep+1] + ds;

for( ;; ) {
    do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = T[iStep];
        nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, Params.Lfull, pfGravityCoefficients );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

        P[iStep-1] = P[iStep] + ( dPbyds * ds );
        FracDiff = fabs( 1.0 - ( P[iStep] / P[iStep-1] ) );
        if( FracDiff > EPSILON )
        {
            ds /= MAX_VARIATION;
            if( -ds < MIN_DS )
                ds = - MIN_DS;
        }

    } while ( FracDiff > EPSILON && -ds > MIN_DS );

    s[iStep-1] = s[iStep] + ds;
    if( s[iStep-1] < 0.0 ) break;

    T[iStep-1] = T[iStep];
    nH[iStep-1] = P[iStep-1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep-1] );
    ne[iStep-1] = nH[iStep-1];

    ds *= MAX_VARIATION;
    if( ds < -max_ds ) ds = - max_ds;

    iStep--;
    iCHRSteps++;
}

// Shift the chromosphere + TR + corona + TR solution in the array back to element zero
for( i=iStep; i<iStep+iCHRSteps+iTRplusCoronaplusTRSteps-1; i++ )
{
    s[i-iStep] = s[i];
    P[i-iStep] = P[i];
    T[i-iStep] = T[i];
    nH[i-iStep] = nH[i];
    ne[i-iStep] = ne[i];
}

// Right-hand chromosphere

iStep = iTRplusCoronaplusTRSteps + iCHRSteps - 1;

nH[iStep] = Params.n0;
ne[iStep] = nH[iStep];
T[iStep] = Params.T0;
P[iStep] = BOLTZMANN_CONSTANT * ( nH[iStep] + ne[iStep] ) * T[iStep];

ds = s[iStep-1] - s[iStep-2];
s[iStep] = s[iStep-1] + ds;

for( ;; ) {
    do {
// *****************************************************************************
// *    STEP 1                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], Params.Lfull, pfGravityCoefficients );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = T[iStep];
        nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, Params.Lfull, pfGravityCoefficients );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

        P[iStep+1] = P[iStep] + ( dPbyds * ds );
        FracDiff = fabs( 1.0 - ( P[iStep+1] / P[iStep] ) );
        if( FracDiff > EPSILON )
        {
            ds /= MAX_VARIATION;
            if( ds < MIN_DS )
		ds = MIN_DS;
        }

    } while ( FracDiff > EPSILON && ds > MIN_DS );

    s[iStep+1] = s[iStep] + ds;
    if( s[iStep+1] > Params.Lfull ) break;

    T[iStep+1] = T[iStep];
    nH[iStep+1] = P[iStep+1] / ( 2.0 * BOLTZMANN_CONSTANT * T[iStep+1] );
    ne[iStep+1] = nH[iStep+1];

    ds *= MAX_VARIATION;
    if( ds > max_ds ) ds = max_ds;

    iStep++;
}

iTotalSteps = iStep + 1;

return iTotalSteps;
}
#endif // OPTICALLY_THICK_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void RecalculateElectronDensity( PARAMETERS Params )
{
    void CalculateColumnMass( double fL, int iEntries, double *ps, double *prho, double *pMc );
    void GetZ_c( double TeZ_c, double *Z_c, int iEntries, double *ps, double *pTe );
    void GetMcZ_c( double *Z_c, double *McZ_c, int iEntries, double *ps, double *pMc );

    FILE *pAMRFile, *pTRANSFile, *pTERMSFile;
    FILE *pOUTPUTFile;
    char szAMRFilename[512], szBuffer[512];
    double *ps, *pds, *prho, *pne, *pnH, *pT, *pMc;
    double fEH;
    double *pTrt, *pnu0, *pTeZ_c, *pZ_c_LEFT, *pZ_c_RIGHT, *pMcZ_c_LEFT, *pMcZ_c_RIGHT;
    double *pterm1, *pterm2, *pH;
    double **ppTrad;
    double fL, Z_c[2], McZ_c[2], fA;
    double fPreviousIteration, fBuffer;
    int iEntries, iNBBT, iNBFT;
    int i, j;
    int iBuffer;
    
    PRADIATION pRadiationEQ, pRadiationNEQ;
    PRADIATIVERATES pRadiativeRates;
    double fBB_lu[6], fBB_ul[6], fBF[4], fFB[4], fColl_ex_lu[10], fColl_ex_ul[10], fColl_ion[5], fColl_rec[5];
    double fradT[10];
    double fZ[31];
    double fSum, fElement;
    double fY[6];
    int *piA, iNumElements, h, iIteration;

    printf( "Recalculating the chromospheric electron density\n\n" );

    // Open the .amr file
    sprintf( szAMRFilename, "%s.orig", Params.szOutputFilename );
    pAMRFile = fopen( szAMRFilename, "r");
    
	ReadDouble( pAMRFile, &fBuffer );		// Timestamp
	fscanf( pAMRFile, "%i", &iBuffer );		// File number
	ReadDouble( pAMRFile, &fL );			// Loop length
	fscanf( pAMRFile, "%i", &iEntries );		// Number of grid cells

	// Allocate memory to hold the position, hydrogen number density, and electron temperature profiles
	ps = (double*)malloc( sizeof(double) * iEntries );
	pds = (double*)malloc( sizeof(double) * iEntries );
	prho = (double*)malloc( sizeof(double) * iEntries );
	pne = (double*)malloc( sizeof(double) * iEntries );
	pnH = (double*)malloc( sizeof(double) * iEntries );
	pT = (double*)malloc( sizeof(double) * iEntries );
	pMc = (double*)malloc( sizeof(double) * iEntries );

	for( i=0; i<iEntries; i++ )
	{
	    // .amr file
	    ReadDouble( pAMRFile, &(ps[i]) );		// Position
	    ReadDouble( pAMRFile, &(pds[i]) );		// Cell width
	    ReadDouble( pAMRFile, &(prho[i]) );		// Mass density
	    ReadDouble( pAMRFile, &fBuffer );		// Momentum density
	    ReadDouble( pAMRFile, &fBuffer );		// Electron energy density
	    ReadDouble( pAMRFile, &fEH );		// Hydrogen energy density

	    fscanf( pAMRFile, "%i", &iBuffer );		// Refinement level
	    for( j=0; j<MAX_REFINEMENT_LEVEL; j++ )
		ReadDouble( pAMRFile, &fBuffer );	// Unique cell identifier

	    pnH[i] = prho[i] / AVERAGE_PARTICLE_MASS;
	    pne[i] = pnH[i];
	    pT[i] = fEH / ( 1.5 * BOLTZMANN_CONSTANT * pnH[i] );
	}

    fclose( pAMRFile );

    // Calculate the column mass density along the loop from the apex to each foot-point
    CalculateColumnMass( fL, iEntries, ps, prho, pMc );

    // Get the information about the transitions
    pTRANSFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/radiative_rates/transitions.txt", "r" );

	fscanf( pTRANSFile, "%i", &iNBBT );		// Number of bound-bound transitions
	fscanf( pTRANSFile, "%i", &iNBFT );		// Number of bound-free transitions

	// Allocate memory for the transition data
	pTrt = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pnu0 = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pTeZ_c  = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pZ_c_LEFT = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pZ_c_RIGHT = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pMcZ_c_LEFT = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pMcZ_c_RIGHT = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );

	// Column labels
	fscanf( pTRANSFile, "%s", szBuffer );
	fscanf( pTRANSFile, "%s", szBuffer );
	fscanf( pTRANSFile, "%s", szBuffer );
	fscanf( pTRANSFile, "%s", szBuffer );
	fscanf( pTRANSFile, "%s", szBuffer );

	for( i=0; i<iNBBT+iNBFT; i++ )
	{
	    fscanf( pTRANSFile, "%i", &iBuffer );	// i
	    fscanf( pTRANSFile, "%i", &iBuffer );	// j
	    ReadDouble( pTRANSFile, &(pTrt[i]) );	// T_rad^top
	    ReadDouble( pTRANSFile, &(pnu0[i]) );	// Transition rest frequency (wavelength)
	    ReadDouble( pTRANSFile, &(pTeZ_c[i]) );	// Temperature at Z_c
	}

    fclose( pTRANSFile );

    // Get the equation terms needed to calculate T_rad for each transition as a function of 's'
    pTERMSFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/radiative_rates/terms.txt", "r" );

	fscanf( pTERMSFile, "%i", &iNBBT );		// Number of bound-bound transitions
	fscanf( pTERMSFile, "%i", &iNBFT );		// Number of bound-free transitions

	// Allocate memory for the transition data
	pterm1 = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pterm2 = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
	pH = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );

	// Column labels
	fscanf( pTERMSFile, "%s", szBuffer );
	fscanf( pTERMSFile, "%s", szBuffer );
	fscanf( pTERMSFile, "%s", szBuffer );

	for( i=0; i<iNBBT+iNBFT; i++ )
	{
	    ReadDouble( pTERMSFile, &(pterm1[i]) );
	    ReadDouble( pTERMSFile, &(pterm2[i]) );
	    ReadDouble( pTERMSFile, &(pH[i]) );
	}

    fclose( pTERMSFile );

    // Calculate Z_c for each of the transitions
    for( i=0; i<iNBBT+iNBFT; i++ )
    {
    	GetZ_c( pTeZ_c[i], Z_c, iEntries, ps, pT );
        pZ_c_LEFT[i] = Z_c[0];
        pZ_c_RIGHT[i] = Z_c[1];
    }

    // Calculate Mc(Z_c) for each of the transitions
    for( i=0; i<iNBBT+iNBFT; i++ )
    {
        Z_c[0] = pZ_c_LEFT[i];
        Z_c[1] = pZ_c_RIGHT[i];
        GetMcZ_c( Z_c, McZ_c, iEntries, ps, pMc );
        pMcZ_c_LEFT[i] = McZ_c[0];
        pMcZ_c_RIGHT[i] = McZ_c[1];
    }

    // Calculate Trad for each transition in each grid cell
    // Allocate a 2D array of size ( # transitions ) x ( # grid cells )
    ppTrad = (double**)malloc( sizeof(double*) * (iNBBT+iNBFT) );
    for( i=0; i<iNBBT+iNBFT; i++ )
        ppTrad[i] = (double*)malloc( sizeof(double) * iEntries );

    for( i=0; i<iNBBT+iNBFT; i++ )
        for( j=0; j<iEntries; j++ )
        {
            if( ps[j] <= fL/2.0 )
            {
                if( ps[j] < pZ_c_LEFT[i] )
                    ppTrad[i][j] = pT[j];
                else
                {
                    fA = pterm1[i] + ( pterm2[i] * pow( (pMc[j]/pMcZ_c_LEFT[i]), pH[i] ) );
                    ppTrad[i][j] = ( (4.7979e-11) * pnu0[i] ) / log( (1.0/fA) + 1.0 );
                }
            }
            else
            {
                if( ps[j] > pZ_c_RIGHT[i] )
                    ppTrad[i][j] = pT[j];
                else
                {
                    fA = pterm1[i] + ( pterm2[i] * pow( (pMc[j]/pMcZ_c_RIGHT[i]), pH[i] ) );
                    ppTrad[i][j] = ( (4.7979e-11) * pnu0[i] ) / log( (1.0/fA) + 1.0 );
                }
            }
        }

    // Update the electron density profile

    pRadiationEQ = new CRadiation( (char *)"Radiation_Model/config/elements_eq.cfg" );
    pRadiationNEQ = new CRadiation( (char *)"Radiation_Model/config/elements_neq.cfg" );
    pRadiativeRates = new CRadiativeRates( (char *)"Radiation_Model/atomic_data/OpticallyThick/radiative_rates/rates_files.txt" );

    for( j=0; j<iEntries; j++ )
    {
	if( !(j%1000) )
		printf( "%.3g%%\n", 100.0*((double)j/(double)iEntries) );

        // Calculate the contribution to the electron density from elements other than hydrogen:
        fSum = 0.0;
        
	// Elements treated in equilibrium (except hydrogen)
	piA = pRadiationEQ->pGetAtomicNumbers( &iNumElements );
	for( i=0; i<iNumElements; i++ )
    	{
	    // Don't double count hydrogen
	    if( piA[i] > 1 )
	    {
	        fElement = 0.0;
			
	        pRadiationEQ->GetEquilIonFrac( piA[i], fZ, log10(pT[j]) );
	        for( h=1; h<piA[i]+1; h++ )
		    fElement += ((double)h) * fZ[h];

                fElement *= pRadiationEQ->GetAbundance( piA[i] );

                fSum += fElement;
	    }
        }

	// Elements treated out-of-equilibrium (except hydrogen)
	piA = pRadiationNEQ->pGetAtomicNumbers( &iNumElements );
	for( i=0; i<iNumElements; i++ )
    	{
	    // Don't double count hydrogen
	    if( piA[i] > 1 )
	    {
	        fElement = 0.0;
			
	        pRadiationNEQ->GetEquilIonFrac( piA[i], fZ, log10(pT[j]) );
	        for( h=1; h<piA[i]+1; h++ )
		    fElement += ((double)h) * fZ[h];

                fElement *= pRadiationNEQ->GetAbundance( piA[i] );

                fSum += fElement;
	    }
        }

	// Convert the radiation temperatures for each transition to log values
	for( i=0; i<iNBBT+iNBFT; i++ )
                fradT[i] = log10( ppTrad[i][j] );

        // Iterate until the solution converges
        for( iIteration=0; iIteration<=MAX_ITERATIONS; iIteration++ )
        {
	    fPreviousIteration = pne[j];

            pRadiativeRates->GetBoundBoundRates( fBB_lu, fBB_ul, &(fradT[0]) );
            pRadiativeRates->GetBoundFreeRates( fBF, &(fradT[6]) );
            pRadiativeRates->GetFreeBoundRates( fFB, &(fradT[6]), log10( pT[j] ), pne[j] );
            pRadiativeRates->GetCollisionalRatesRH( fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, log10( pT[j] ), pne[j] );
    
            pRadiativeRates->SolveHIIFraction( fY, fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, fBB_lu, fBB_ul, fBF, fFB );
                
            // Recalculate the electron density. To improve stability, iterate gradually to the converged solution rather than taking large steps
            pne[j] = pne[j] - CONVERGENCE_EPSILON * ( pne[j] - ( pnH[j] * ( fY[5] + fSum ) ) );
	    
	    // Check for convergence and exit the loop if the condition is met
            if( iIteration && (fabs(fPreviousIteration-pne[j])/pne[j]) < CONVERGENCE_CONDITION ) break;
        }
    }
    
    printf( "\nWriting new .amr file\n\n" );

    // Update the .amr file. The new .amr file contains a new column corresponding to the electron mass density
    // Open the .amr and .phy files
    pAMRFile = fopen( szAMRFilename, "r");
    pOUTPUTFile = fopen( Params.szOutputFilename, "w" );
    
	ReadDouble( pAMRFile, &fBuffer );	fprintf( pOUTPUTFile, "%g\r\n", fBuffer );		// Timestamp
	fscanf( pAMRFile, "%i", &iBuffer );	fprintf( pOUTPUTFile, "%i\r\n", iBuffer );		// File number
	ReadDouble( pAMRFile, &fL );		fprintf( pOUTPUTFile, "%g\r\n", fL );			// Loop length
	fscanf( pAMRFile, "%i", &iEntries );	fprintf( pOUTPUTFile, "%i\r\n", iEntries );		// Number of grid cells

	for( i=0; i<iEntries; i++ )
	{
	    // .amr file
	    ReadDouble( pAMRFile, &(ps[i]) );	fprintf( pOUTPUTFile, "%.16e\t", ps[i] );					// Position
	    ReadDouble( pAMRFile, &(pds[i]) );	fprintf( pOUTPUTFile, "%.16e\t", pds[i] );					// Cell width

	    ReadDouble( pAMRFile, &(prho[i]) );	fprintf( pOUTPUTFile, "%.16e\t", ELECTRON_MASS * pne[i] );			// Electron mass density
						fprintf( pOUTPUTFile, "%.16e\t", prho[i] );					// Hydrogen mass density

	    ReadDouble( pAMRFile, &fBuffer );	fprintf( pOUTPUTFile, "%.16e\t", fBuffer );					// Momentum density

	    ReadDouble( pAMRFile, &fBuffer );	fprintf( pOUTPUTFile, "%.16e\t", 1.5 * BOLTZMANN_CONSTANT * pne[i] * pT[i] );	// Electron energy density
	    ReadDouble( pAMRFile, &fBuffer );	fprintf( pOUTPUTFile, "%.16e\t", fBuffer );					// Hydrogen energy density

	    fscanf( pAMRFile, "%i", &iBuffer );	fprintf( pOUTPUTFile, "%i", iBuffer );			// Refinement level
	    for( j=0; j<MAX_REFINEMENT_LEVEL; j++ )
	    {
		ReadDouble( pAMRFile, &fBuffer ); fprintf( pOUTPUTFile, "\t%i", (int)fBuffer );		// Unique cell identifier
	    }
	    fprintf( pOUTPUTFile, "\r\n" );
	}

    fclose( pOUTPUTFile );
    fclose( pAMRFile );
    
    delete pRadiativeRates;
    delete pRadiationNEQ;
    delete pRadiationEQ;

    // Free the allocated memory
    for( i=0; i<iNBBT+iNBFT; i++ )
        free( ppTrad[i] );
    free( ppTrad );

    free( pH );
    free( pterm2 );
    free( pterm1 );

    free( pMcZ_c_RIGHT );
    free( pMcZ_c_LEFT );
    free( pZ_c_RIGHT );
    free( pZ_c_LEFT );
    free( pTeZ_c );
    free( pnu0 );
    free( pTrt );

    free( pMc );
    free( pT );
    free( pnH );
    free( pne );
    free( prho );
    free( pds );
    free( ps );
}

void CalculateColumnMass( double fL, int iEntries, double *ps, double *prho, double *pMc )
{
    int iApex, i;

    for( i=0; i<iEntries; i++ )
        if( ps[i] >= fL/2.0 ) break;

    // Note that this is the column mass density, but a ratio of column masses appears in the calculation and so the area factor cancels
    iApex = i;
    pMc[iApex] = prho[iApex] * (ps[iApex] - ps[iApex-1]);

    
    // Left-hand leg
    for( i=iApex-1; i>=0; i-- )
        pMc[i] = pMc[i+1] + ( prho[i] * (ps[i+1] - ps[i]) );

    
    // Right-hand leg
    for( i=iApex+1; i<iEntries; i++ )
        pMc[i] = pMc[i-1] + ( prho[i] * (ps[i] - ps[i-1]) );

}

void GetZ_c( double TeZ_c, double *Z_c, int iEntries, double *ps, double *pTe )
{
    int i;
    double x[3], y[3];

    // Left-hand chromosphere
    for( i=0; i<iEntries-1; i++ )
        if( ( TeZ_c <= pTe[i] && TeZ_c >= pTe[i+1] ) || ( TeZ_c >= pTe[i] && TeZ_c <= pTe[i+1] ) )
            break;

    x[1] = pTe[i];
    x[2] = pTe[i+1];
    y[1] = ps[i];
    y[2] = ps[i+1];
    LinearFit( x, y, TeZ_c, &(Z_c[0]) );

    // Right-hand chromosphere
    for( i=iEntries-2; i>=0; i-- )
        if( ( TeZ_c <= pTe[i] && TeZ_c >= pTe[i+1] ) || ( TeZ_c >= pTe[i] && TeZ_c <= pTe[i+1] ) )
            break;

    x[1] = pTe[i];
    x[2] = pTe[i+1];
    y[1] = ps[i];
    y[2] = ps[i+1];
    LinearFit( x, y, TeZ_c, &(Z_c[1]) );
}

void GetMcZ_c( double *Z_c, double *McZ_c, int iEntries, double *ps, double *pMc )
{
    int i;
    double x[3], y[3];

    // Left-hand chromosphere
    for( i=0; i<iEntries-1; i++ )
        if( Z_c[0] >= ps[i] && Z_c[0] <= ps[i+1] )
            break;

    x[1] = ps[i];
    x[2] = ps[i+1];
    y[1] = pMc[i];
    y[2] = pMc[i+1];
    LinearFit( x, y, Z_c[0], &(McZ_c[0]) );

    // Right-hand chromosphere
    for( i=iEntries-2; i>=0; i-- )
        if( Z_c[1] >= ps[i] && Z_c[1] <= ps[i+1] )
            break;

    x[1] = ps[i];
    x[2] = ps[i+1];
    y[1] = pMc[i];
    y[2] = pMc[i+1];
    LinearFit( x, y, Z_c[1], &(McZ_c[1]) );
}

int AMR2PHY( PARAMETERS Params )
{
    FILE *pAMRFile, *pPHYFile;
    char szPHYFilename[512];
    double *pfs, *pfds, frho_e, frho_H, frhov, fE_e, fE_H;
    double *pfv, *pfCs, *pfne, *pfnH, *pfPe, *pfPH, *pfTe, *pfTH, *pfF_ce, *pfF_cH;
    double flast_s, flast_ne, flast_nH, flast_Te;
    double fDiff_ne[2], fDiff_nH[2], fDiff_Te[2], fMAX_DS = Params.Lfull / MIN_CELLS;
    double delta = MAX_VARIATION - 1.0, kappa, s[2], T[2], x[4], y[4];
    double fBuffer;
    int iEntries, iTotalSteps;
    int i, j, iBuffer;
    
    printf( "Writing new .phy file\n\n" );

    // Open the .amr file
    pAMRFile = fopen( Params.szOutputFilename, "r");
    
	ReadDouble( pAMRFile, &fBuffer );		// Timestamp
	fscanf( pAMRFile, "%i", &iBuffer );		// File number
	ReadDouble( pAMRFile, &fBuffer );		// Loop length
	fscanf( pAMRFile, "%i", &iEntries );		// Number of grid cells

	// Allocate memory to hold the quantities to store in the .phy file
	pfs = (double*)malloc( sizeof(double) * iEntries );
	pfds = (double*)malloc( sizeof(double) * iEntries );
	pfv = (double*)malloc( sizeof(double) * iEntries );
	pfCs = (double*)malloc( sizeof(double) * iEntries );
	pfne = (double*)malloc( sizeof(double) * iEntries );
	pfnH = (double*)malloc( sizeof(double) * iEntries );
	pfPe = (double*)malloc( sizeof(double) * iEntries );
	pfPH = (double*)malloc( sizeof(double) * iEntries );
	pfTe = (double*)malloc( sizeof(double) * iEntries );
	pfTH = (double*)malloc( sizeof(double) * iEntries );
	pfF_ce = (double*)malloc( sizeof(double) * iEntries );
	pfF_cH = (double*)malloc( sizeof(double) * iEntries );

	for( i=0; i<iEntries; i++ )
	{
	    // .amr file
	    ReadDouble( pAMRFile, &(pfs[i]) );		// Position
	    ReadDouble( pAMRFile, &(pfds[i]) );		// Cell width
	    ReadDouble( pAMRFile, &frho_e );		// Electron mass density
	    ReadDouble( pAMRFile, &frho_H );		// Hydrogen mass density
	    ReadDouble( pAMRFile, &frhov );		// Momentum density
	    ReadDouble( pAMRFile, &fE_e );		// Electron energy density
	    ReadDouble( pAMRFile, &fE_H );		// Hydrogen energy density

	    fscanf( pAMRFile, "%i", &iBuffer );		// Refinement level
	    for( j=0; j<MAX_REFINEMENT_LEVEL; j++ )
		ReadDouble( pAMRFile, &fBuffer );	// Unique cell identifier

	    // Calculate the quantities to store in the .phy file
	    pfv[i] = frhov / frho_H;            
	    pfne[i] = frho_e / ELECTRON_MASS;
	    pfnH[i] = frho_H / AVERAGE_PARTICLE_MASS;
	    pfTe[i] = fE_e / ( 1.5 * BOLTZMANN_CONSTANT * pfne[i] );
	    pfTH[i] = ( fE_H - ( 0.5 * frhov * pfv[i] ) ) / ( 1.5 * BOLTZMANN_CONSTANT * pfnH[i] );
	    pfPe[i] = BOLTZMANN_CONSTANT * pfne[i] * pfTe[i];
    	    pfPH[i] = BOLTZMANN_CONSTANT * pfnH[i] * pfTH[i];
	    pfCs[i] = sqrt( ( GAMMA * ( pfPe[i] + pfPH[i] ) ) / frho_H );
	    pfF_ce[i] = 0.0;
	    pfF_cH[i] = 0.0;
	}

    fclose( pAMRFile );

    // Open the .phy file
    sprintf( szPHYFilename, "%s.phy", Params.szOutputFilename );
    pPHYFile = fopen( szPHYFilename, "w");
    
	i = 0; fprintf( pPHYFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", pfs[i], pfv[i], pfCs[i], pfne[i], pfnH[i], pfPe[i], pfPH[i], pfTe[i], pfTH[i], pfF_ce[i], pfF_cH[i] );
	i = 1; fprintf( pPHYFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", pfs[i], pfv[i], pfCs[i], pfne[i], pfnH[i], pfPe[i], pfPH[i], pfTe[i], pfTH[i], pfF_ce[i], pfF_cH[i] );
		flast_s = pfs[i]; flast_ne = pfne[i]; flast_nH = pfnH[i]; flast_Te = pfTe[i];

	iTotalSteps = 2;

        for( i=2; i<iEntries-2; i++ )
	{

	    fDiff_ne[0] = fabs( pfne[i] - flast_ne ) / flast_ne;
	    fDiff_ne[1] = fDiff_ne[0] * ( flast_ne / pfne[i] );
	    fDiff_nH[0] = fabs( pfnH[i] - flast_nH ) / flast_nH;
	    fDiff_nH[1] = fDiff_nH[0] * ( flast_nH / pfnH[i] );
	    fDiff_Te[0] = fabs( pfTe[i] - flast_Te ) / flast_Te;
	    fDiff_Te[1] = fDiff_Te[0] * ( flast_Te / pfTe[i] );

	    if( fDiff_ne[0] > delta || fDiff_ne[1] > delta || fDiff_nH[0] > delta || fDiff_nH[1] > delta || fDiff_Te[0] > delta || fDiff_Te[1] > delta || ( pfs[i] - flast_s ) >= fMAX_DS )
	    {
	        // Electron thermal conduction
	        kappa = SPITZER_ELECTRON_CONDUCTIVITY * pow( pfTe[i], 2.5 );
	        x[1] = pfs[i-1]; x[2] = pfs[i]; x[3] = pfs[i+1];
	        y[1] = pfTe[i-1]; y[2] = pfTe[i]; y[3] = pfTe[i+1];
	        s[0] = pfs[i] - ( 0.5 * pfds[i] ); s[1] = s[0] + pfds[i];
	        LinearFit( x, y, s[0], &(T[0]) );
	        LinearFit( &(x[1]), &(y[1]), s[1], &(T[1]) );
	        pfF_ce[i] = - kappa * ( ( T[1] - T[0] ) / pfds[i] );

	        // Hydrogen thermal conduction
	        kappa = SPITZER_ION_CONDUCTIVITY * pow( pfTH[i], 2.5 );
	        y[1] = pfTH[i-1]; y[2] = pfTH[i]; y[3] = pfTH[i+1];
	        LinearFit( x, y, s[0], &(T[0]) );
	        LinearFit( &(x[1]), &(y[1]), s[1], &(T[1]) );
	        pfF_cH[i] = - kappa * ( ( T[1] - T[0] ) / pfds[i] );

	        fprintf( pPHYFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", pfs[i], pfv[i], pfCs[i], pfne[i], pfnH[i], pfPe[i], pfPH[i], pfTe[i], pfTH[i], pfF_ce[i], pfF_cH[i] );
			flast_s = pfs[i]; flast_ne = pfne[i]; flast_nH = pfnH[i]; flast_Te = pfTe[i];

		iTotalSteps++;
	    }
	}
	
	i = iEntries - 2; fprintf( pPHYFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", pfs[i], pfv[i], pfCs[i], pfne[i], pfnH[i], pfPe[i], pfPH[i], pfTe[i], pfTH[i], pfF_ce[i], pfF_cH[i] );
	i = iEntries - 1; fprintf( pPHYFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", pfs[i], pfv[i], pfCs[i], pfne[i], pfnH[i], pfPe[i], pfPH[i], pfTe[i], pfTH[i], pfF_ce[i], pfF_cH[i] );

	iTotalSteps += 2;

    fclose( pPHYFile );

    // Free the allocated memory
    free( pfF_cH );
    free( pfF_ce );
    free( pfTH );
    free( pfTe );
    free( pfPH );
    free( pfPe );
    free( pfnH );
    free( pfne );
    free( pfCs );
    free( pfv );
    free( pfds );
    free( pfs );

    return iTotalSteps;
}

void RecalculateChromosphericHeating( PARAMETERS Params, int number_of_lines )
{
    char szPHYFilename[512];
    double s, v, Cs, n_e, n_H, P_e, P_H, T_e, T_H, Fce, FcH;
    double frho_c = 0., fHI_c =0., previous_s = Params.Lfull/2., cell_width_cos_theta = 0., fDensityDifference = 0., fRadiation=0.;
    double fSum, fElement, fZ[31];
    int *piA, iNumElements;
    int num_elements = 0;
    int i, j, k;

    PRADIATION pRadiation_EQ;
    pRadiation_EQ = new CRadiation( "Radiation_Model/config/elements_eq.cfg" );
#if defined (DECOUPLE_IONISATION_STATE_SOLVER) || defined (NON_EQUILIBRIUM_RADIATION)
    PRADIATION pRadiation_NEQ;
    pRadiation_NEQ = new CRadiation( "Radiation_Model/config/elements_neq.cfg" );
#endif // DECOUPLE_IONISATION_STATE_SOLVER || NON_EQUILIBRIUM_RADIATION

    COpticallyThickIon *pHI, *pMgII,*pCaII;   
    pHI = new COpticallyThickIon( 1, (char *)"h_1", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    pMgII = new COpticallyThickIon( 12, (char *)"mg_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    pCaII = new COpticallyThickIon( 20, (char *)"ca_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    
    // Create vectors to store the log of radiation and column density
    vector<double> fRad;
    vector<double> fRho;

    printf( "Recalculating the chromospheric heating\n\n" );

    // Reserve space for the maximum size of the arrays
    fRad.reserve(number_of_lines/2);
    fRho.reserve(number_of_lines/2);

    sprintf( szPHYFilename, "%s.phy", Params.szOutputFilename );    
    ifstream phy_infile( szPHYFilename );

    for( k=0;k<number_of_lines;++k )
    {
        // Read in the line and values on the line
        phy_infile >> s >> v >> Cs >> n_e >> n_H >> P_e >> P_H >> T_e >> T_H >> Fce >> FcH;
        
        // For positions in the second half of the loop, calculate the column densities
        if( k > number_of_lines/2 )
        {
            
            if( T_e < OPTICALLY_THICK_TEMPERATURE )
            {
                cell_width_cos_theta = fabs(s - previous_s) * fabs( cos( ( _PI_ * s ) / Params.Lfull ) );
                
                frho_c += n_H * AVERAGE_PARTICLE_MASS * cell_width_cos_theta;

#ifdef NLTE_CHROMOSPHERE
		fSum = 0.0;
		piA = pRadiation_EQ->pGetAtomicNumbers( &iNumElements );
		for( i=0; i<iNumElements; i++ )
		{
			// Don't double count hydrogen
			if( piA[i] > 1 )
			{
				fElement = 0.0;
			
				pRadiation_EQ->GetEquilIonFrac( piA[i], fZ, log10(T_e) );
				for( j=1; j<piA[i]+1; j++ )
					fElement += ((double)j) * fZ[j];

				fElement *= pRadiation_EQ->GetAbundance( piA[i] );

				fSum += fElement;
			}
		}
#if defined (DECOUPLE_IONISATION_STATE_SOLVER) || defined (NON_EQUILIBRIUM_RADIATION)
		piA = pRadiation_NEQ->pGetAtomicNumbers( &iNumElements );
		for( i=0; i<iNumElements; i++ )
		{
			// Don't double count hydrogen
			if( piA[i] > 1 )
			{
				fElement = 0.0;
			
				pRadiation_NEQ->GetEquilIonFrac( piA[i], fZ, log10(T_e) );
				for( j=1; j<piA[i]+1; j++ )
					fElement += ((double)j) * fZ[j];

				fElement *= pRadiation_NEQ->GetAbundance( piA[i] );

				fSum += fElement;
			}
		}
#endif // DECOUPLE_IONISATION_STATE_SOLVER || NON_EQUILIBRIUM_RADIATION

		fHI_c += ( n_H * ( 1.0 + fSum ) - n_e ) * cell_width_cos_theta;
#else // NLTE_CHROMOSPHERE
                // Calculate the neutral hydrogen column number density
                fDensityDifference = n_H - ( n_e / 1.000144 );
                if ( fabs(fDensityDifference)  < 1.0 )
                    fHI_c += 0.0;
                else
                    fHI_c += ( n_H - ( n_e / 1.000144 ) ) * cell_width_cos_theta;
#endif // NLTE_CHROMOSPHERE
                
                fRho.push_back(log10(frho_c));
                
                fRadiation = pHI->GetVolumetricLossRate( log10(T_e), log10((4e-14)*fHI_c), n_e * n_H * AVERAGE_PARTICLE_MASS) + pMgII->GetVolumetricLossRate( log10(T_e), log10(frho_c), n_e * n_H * AVERAGE_PARTICLE_MASS) + pCaII->GetVolumetricLossRate( log10(T_e), log10(frho_c), n_e * n_H * AVERAGE_PARTICLE_MASS);
             
                fRad.push_back(log10(fRadiation));
                
                ++num_elements;         
            }
        }
        
        previous_s = s;
    }
    
    phy_infile.close();

    delete pHI;
    delete pMgII;
    delete pCaII;

#if defined (DECOUPLE_IONISATION_STATE_SOLVER) || defined (NON_EQUILIBRIUM_RADIATION)
    delete pRadiation_NEQ;
#endif // DECOUPLE_IONISATION_STATE_SOLVER || NON_EQUILIBRIUM_RADIATION
    delete pRadiation_EQ;

    // Output the volumetric heating rate vs. column density
    ofstream outfile("Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat");
    outfile << num_elements << endl;
    outfile.precision(10);
    for( i=0; i<num_elements; ++i)
    {
        if( fRho[i] >= 0.0 ) fRad[i] = -300.0;
        outfile << fRho[i] << "\t" << fRad[i] << endl;
    }
    outfile.close();
}
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION