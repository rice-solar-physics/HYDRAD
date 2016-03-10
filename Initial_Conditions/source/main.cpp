// ****
// *
// * A hydrostatic code that calculates a set of initial conditions 
// * for the hydrodynamic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/02/2016
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
#include  "../../Resources/source/fitpoly.h"
#endif // OPTICALLY_THICK_RADIATION
#include "../../Resources/source/file.h"
#include "../../Resources/source/constants.h"


int main(void)
{
double FindHeatingRange( double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity );
double RefineSolution( double Log_10H0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity );
int CalculateSolution( double finalH0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity );
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, int igdp, double **ppGravity );

PARAMETERS Params;
PRADIATION pRadiation;

FILE *pFile;
#ifdef USE_TABULATED_GRAVITY
#else // USE_TABULATED_GRAVITY
char szGravityFilename[256];
#endif // USE_TABULATED_GRAVITY
double **ppGravity;
int i, igdp;

double *s, *P, *T, *nH, *ne;
double Log_10H0, finalH0;

int iTRplusCoronaplusTRSteps, iTotalSteps;

printf( "\n\nCalculating initial hydrostatic conditions...\n\n" );

// Get the user-specified parameter values from the configuration file
GetConfigurationParameters( &Params );

// Initialise the radiative losses
pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_eq.cfg" );

// Initialise the gravitational geometry
#ifdef USE_TABULATED_GRAVITY
pFile = fopen( TABULATED_GRAVITY_FILE, "r" );
#else // USE_TABULATED_GRAVITY
GenerateSemiCircularLoop( Params );
sprintf( szGravityFilename, "%s.gravity", Params.szOutputFilename );
pFile = fopen( szGravityFilename, "r" );
#endif // USE_TABULATED_GRAVITY
fscanf( pFile, "%i", &igdp );
ppGravity = (double**)malloc( igdp * sizeof( double ) );
for( i=0; i<igdp; i++ )
{
    ppGravity[i] = (double*)malloc( 2 * sizeof( double ) );
    ReadDouble( pFile, &(ppGravity[i][0]) );
    ReadDouble( pFile, &(ppGravity[i][1]) );
}
fclose( pFile );

// Allocate memory to store the hydrostatic profiles
s = (double*)malloc( sizeof(double) * MAX_CELLS );
P = (double*)malloc( sizeof(double) * MAX_CELLS );
T = (double*)malloc( sizeof(double) * MAX_CELLS );
nH = (double*)malloc( sizeof(double) * MAX_CELLS );
ne = (double*)malloc( sizeof(double) * MAX_CELLS );

Log_10H0 = FindHeatingRange( s, P, T, nH, ne, Params, pRadiation, igdp, ppGravity );
finalH0 = RefineSolution( Log_10H0, s, P, T, nH, ne, Params, pRadiation, igdp, ppGravity );
iTRplusCoronaplusTRSteps = CalculateSolution( finalH0, s, P, T, nH, ne, Params, pRadiation, igdp, ppGravity );
iTotalSteps = AddChromospheres( iTRplusCoronaplusTRSteps, s, P, T, nH, ne, Params, igdp, ppGravity );

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
for( i=0; i<igdp; i++ )
    free( ppGravity[i]);
free( ppGravity );

// Close the radiative losses
delete pRadiation;

printf( "Done!\n\n" );

return 0;
}

double FindHeatingRange( double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, nH2;
double Log_10H0, H0, H, R;
double dPbyds, dFcbyds, dTbyds;
double FracDiff;
int iStep;

max_ds = Params.Lfull / MIN_CELLS;

sL = Params.s0;
sR = Params.Lfull - sL;

// First find the lower- and upper-limits of the heating range
for( Log_10H0=Params.Log_10H0[0]; Log_10H0<=Params.Log_10H0[1]; Log_10H0+=Params.dLog_10H0 )
{
    iStep = 0;

#ifdef ISOTHERMAL
    H0 = 0.0;
#else // ISOTHERMAL
    H0 = pow( 10.0, Log_10H0 );
#endif // ISOTHERMAL

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
            dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

            // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
            H = 0.0;
            R = 0.0;
#else // ISOTHERMAL
            H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), log10( nH[iStep] ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T[iStep] ), log10( nH[iStep] ) ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), log10( nH[iStep] ) ) );
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

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, igdp, ppGravity );

            // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
            H = 0.0;
            R = 0.0;
#else // ISOTHERMAL
            H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T2 ), log10( nH2 ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T2 ), log10( nH2 ) ) + pRadiation->GetFreeFreeRad( log10( T2 ), log10( nH2 ) ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#endif // ISOTHERMAL

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

            T[iStep+1] = T[iStep] + ( dTbyds * ds );
            FracDiff = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
            if( FracDiff > EPSILON )
            {
		ds /= MAX_VARIATION;
		if( ds < MIN_DS )
		{
			ds = MIN_DS;
		}
            }
        } while ( FracDiff > EPSILON && ds > MIN_DS );

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

double RefineSolution( double Log_10H0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, nH2;
double H0, H, R;
double dPbyds, dFcbyds, dTbyds;
double FracDiff;
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
            dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

            // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
            H = 0.0;
            R = 0.0;
#else // ISOTHERMAL
            H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), log10( nH[iStep] ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T[iStep] ), log10( nH[iStep] ) ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), log10( nH[iStep] ) ) );
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

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

            // Get the pressure gradient
            dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, igdp, ppGravity );

            // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
            H = 0.0;
            R = 0.0;
#else // ISOTHERMAL
            H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
            R = - pRadiation->GetPowerLawRad( log10( T2 ), log10( nH2 ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
            R = - ( pRadiation->GetRadiation( log10( T2 ), log10( nH2 ) ) + pRadiation->GetFreeFreeRad( log10( T2 ), log10( nH2 ) ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#endif // ISOTHERMAL

            // Get the heat flux gradient
            dFcbyds = CalcdFcbyds( H, R );

            // Get the temperature gradient
            dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

            T[iStep+1] = T[iStep] + ( dTbyds * ds );
            FracDiff = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
            if( FracDiff > EPSILON )
            {
		ds /= MAX_VARIATION;
		if( ds < MIN_DS )
                    ds = MIN_DS;
            }

        } while ( FracDiff > EPSILON && ds > MIN_DS );

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

int CalculateSolution( double finalH0, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, PRADIATION pRadiation, int igdp, double **ppGravity )
{
double ds, max_ds, sL, sR;
double Fc, P2, Fc2, T2, nH2;
double H0, H, R;
double dPbyds, dFcbyds, dTbyds;
double FracDiff;
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
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

        // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
        H = 0.0;
        R = 0.0;
#else // ISOTHERMAL
        H = Eheat( s[iStep], H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
        R = - pRadiation->GetPowerLawRad( log10( T[iStep] ), log10( nH[iStep] ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
        R = - ( pRadiation->GetRadiation( log10( T[iStep] ), log10( nH[iStep] ) ) + pRadiation->GetFreeFreeRad( log10( T[iStep] ), log10( nH[iStep] ) ) );
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

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, igdp, ppGravity );

        // Calculate the heat input and the radiation
#ifdef ISOTHERMAL
        H = 0.0;
        R = 0.0;
#else // ISOTHERMAL
        H = Eheat( (s[iStep]+(ds/2.0)), H0, Params.sH0, Params.sH );
#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
        R = - pRadiation->GetPowerLawRad( log10( T2 ), log10( nH2 ) );
#else // USE_POWER_LAW_RADIATIVE_LOSSES
        R = - ( pRadiation->GetRadiation( log10( T2 ), log10( nH2 ) ) + pRadiation->GetFreeFreeRad( log10( T2 ), log10( nH2 ) ) );
#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#endif // ISOTHERMAL

        // Get the heat flux gradient
        dFcbyds = CalcdFcbyds( H, R );

        // Get the temperature gradient
        dTbyds = CalcdTbyds( Fc2, T2 );

// *****************************************************************************
// *    STEP 3                                                                 *
// *****************************************************************************

        T[iStep+1] = T[iStep] + ( dTbyds * ds );
        FracDiff = fabs( 1.0 - ( T[iStep+1] / T[iStep] ) );
        if( FracDiff > EPSILON )
        {
            ds /= MAX_VARIATION;
            if( ds < MIN_DS )
		ds = MIN_DS;
        }

    } while ( FracDiff > EPSILON && ds > MIN_DS );

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
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, int igdp, double **ppGravity )
{
double GetVALTemperature( double s, int iVALTemperatureDP, double **ppVALTemperature );

double ds, max_ds;
double P2, T2, nT, nH2, ne2;
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
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = GetVALTemperature( s[iStep]+(ds/2.0), iVALTemperatureDP, ppVALTemperature );
        nT = P2 / ( BOLTZMANN_CONSTANT * T2 );
        nH2 = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T2) ) ) ) ) * nT;
        ne2 = nT - ( ( 1.0 - 1.44e-4 ) * nH2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep]+(ds/2.0), nH2, igdp, ppGravity );

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
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = GetVALTemperature( Params.Lfull - ( s[iStep]+(ds/2.0) ), iVALTemperatureDP, ppVALTemperature );
        nT = P2 / ( BOLTZMANN_CONSTANT * T2 );
        nH2 = ( 1.0 / ( 1.0 + ( 1.0 - pHI->GetIonFrac( log10(T2) ) ) ) ) * nT;
        ne2 = nT - ( ( 1.0 - 1.44e-4 ) * nH2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( s[iStep]+(ds/2.0), nH2, igdp, ppGravity );

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
int AddChromospheres( int iTRplusCoronaplusTRSteps, double *s, double *P, double *T, double *nH, double *ne, PARAMETERS Params, int igdp, double **ppGravity )
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
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = T[iStep];
        nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, igdp, ppGravity );

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
        dPbyds = CalcdPbyds( s[iStep], nH[iStep], igdp, ppGravity );

        P2 = P[iStep] + ( dPbyds * (ds/2.0) );
        T2 = T[iStep];
        nH2 = P2 / ( 2.0 * BOLTZMANN_CONSTANT * T2 );

// *****************************************************************************
// *    STEP 2                                                                 *
// *****************************************************************************

        // Get the pressure gradient
        dPbyds = CalcdPbyds( (s[iStep]+(ds/2.0)), nH2, igdp, ppGravity );

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