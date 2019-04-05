// ****
// *
// * Miscellaneous routines used by the hydrostatic code
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
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/Utils/regPoly/regpoly.h"
#include "../../Resources/Utils/regPoly/nrutil.h"


void GetConfigurationParameters( PARAMETERS *pParams )
{
FILE *pFile;

pFile = fopen( "Initial_Conditions/config/initial_conditions.cfg", "r" );

// Output path and filename
fscanf( pFile, "%s", pParams->szOutputFilename );

// Loop geometry
ReadDouble( pFile, &(pParams->Lfull) );
ReadDouble( pFile, &(pParams->Inc) );
ReadDouble( pFile, &(pParams->s0) );

// Lower boundary conditions
ReadDouble( pFile, &(pParams->T0) );
ReadDouble( pFile, &(pParams->n0) );

// Solution parameter space
// Heating
ReadDouble( pFile, &(pParams->sH0) );
ReadDouble( pFile, &(pParams->sH) );
ReadDouble( pFile, &(pParams->Log_10H0[0]) );
ReadDouble( pFile, &(pParams->Log_10H0[1]) );
ReadDouble( pFile, &(pParams->dLog_10H0) );
ReadDouble( pFile, &(pParams->Hintervals) );

fclose( pFile );
}

#ifdef USE_POLY_FIT_TO_GRAVITY
#else // USE_POLY_FIT_TO_GRAVITY
void GenerateSemiCircularLoop( PARAMETERS Params )
{
double CalcSolarGravity( double s, double Lfull, double Inc );

FILE *pFile;
char szGravityFilename[512];
double s, ds;
int i;

// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
void BasisFuncs( double x, double *bfunc, int ma );
double *x, *y, *sig, *a, chisq;
double **u, **v, *w;
int ndat, ma;
// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****

ndat = MIN_CELLS + 1;

x = vector( 1, ndat );
y = vector( 1, ndat );
sig = vector( 1, ndat );

s = 0.0;
ds = Params.Lfull / ( ndat - 1 );
for( i=1; i<=ndat; i++ )
{
    x[i] = s / Params.Lfull;
    y[i] = CalcSolarGravity( s, Params.Lfull, Params.Inc );
    sig[i] = 1.0;
    s += ds;
}

ma = POLY_ORDER + 1;
a = vector( 1, ma );
u = matrix( 1, ndat, 1, ma );
v = matrix( 1, ma, 1, ma );
w = vector( 1, ma );

// Call the function that returns the coefficients of the best-fit polynomial (Single Value Decomposition)
svdfit( x, y, sig, ndat, a, ma, u, v, w, &chisq, BasisFuncs );

sprintf( szGravityFilename, "%s.gravity", Params.szOutputFilename );
pFile = fopen( szGravityFilename, "w" );
for( i=1; i<=ma; i++ )
    fprintf( pFile, "%.16e\t", a[i] );
fclose( pFile );

free_vector( w, 1, ma );
free_matrix( v, 1, ma, 1, ma );
free_matrix( u, 1, ndat, 1, ma );
free_vector( a, 1, ma );
free_vector( sig, 1, ndat );
free_vector( y, 1, ndat );
free_vector( x, 1, ndat );
}

double CalcSolarGravity( double s, double Lfull, double Inc )
{
double fTheta, fHeight, fPhi1, fPhi2, r, result;

fTheta = ( _PI_ * s ) / Lfull;
fHeight = ( Lfull / _PI_ ) * sin( fTheta );

// Now allow for the possibility of loop inclination away from the perpendicular
fPhi1 = ( _PI_ / 180.0 ) * Inc;
fPhi2 = atan( ( fHeight * tan( fPhi1 ) ) / ( SOLAR_RADIUS + fHeight ) );
r = ( SOLAR_RADIUS + fHeight ) / cos( fPhi2 );

result = - SOLAR_SURFACE_GRAVITY * ( SOLAR_RADIUS_SQUARED / ( r * r ) ) * cos( fTheta ) * cos( fPhi1 );

return result;
}

void BasisFuncs( double x, double *bfunc, int ma )
{
int i=1;

bfunc[i] = 1.0;
for( i=2; i<=ma; i++ )
	bfunc[i] = bfunc[i-1] * x;	
}
#endif // USE_POLY_FIT_TO_GRAVITY

void WriteAMRFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params )
{
FILE *pFile;

double ds, sAMR;
double x[5], y[5], error, fConservedQuantities[3];

long iMAX_CELLS;
#ifdef ADAPT
long **ppiID, iID;
long iRL, iBlockLength;
int j;
#endif // ADAPT
long i, k;

ds = Params.Lfull / MIN_CELLS;
#ifdef ADAPT
ds /= pow( 2.0, MAX_REFINEMENT_LEVEL );
#endif // ADAPT
iMAX_CELLS = (long)(Params.Lfull/ds);

// Ensure the maximum number of grid cells is an even number
if( iMAX_CELLS & 1 )
{
	iMAX_CELLS++;
	ds = Params.Lfull / (double)iMAX_CELLS;
}

#ifdef ADAPT
// Allocate sufficient memory to hold the unique ID numbers that pair grid cells following refinement
ppiID = (long**)malloc( sizeof(long) * iMAX_CELLS );
for( i=0; i<iMAX_CELLS; i++ )
    ppiID[i] = (long*)malloc( sizeof(long) * MAX_REFINEMENT_LEVEL );

// Create the blocks of unique ID numbers for each refinement level
iID = 0;
for( j=0; j<MAX_REFINEMENT_LEVEL; j++ )
{
    iRL = j + 1;
    iBlockLength = pow( 2, ( MAX_REFINEMENT_LEVEL - iRL + 1 ) );
    for( i=0; i<iMAX_CELLS; i+=iBlockLength )
    {
        for( k=0; k<iBlockLength; k++ )
            ppiID[i+k][j] = iID;

        iID++;
    }
}
#endif // ADAPT

sAMR = ds / 2.0;

#if defined (OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE)
char szAMRFilename[512];
sprintf( szAMRFilename, "%s.orig", Params.szOutputFilename );
pFile = fopen( szAMRFilename, "w" );
#else // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
pFile = fopen( Params.szOutputFilename, "w" );
#endif // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
// Write the header information into the .amr file
fprintf( pFile, "0.0\n0\n%.16e\n%ld\n", Params.Lfull, iMAX_CELLS );
i = 2;
#ifdef ADAPT
j = 0;
#endif // ADAPT
for( ;; )
{
    if( sAMR > Params.Lfull ) break;

    while( sAMR > s[i] )
    {
        i++;
        if( i > iTotalSteps - 2 )
        {
            i = iTotalSteps - 2;
            break;
        }
    }

    x[1] = s[i-2];
    x[2] = s[i-1];
    x[3] = s[i];
    x[4] = s[i+1];

    y[1] = nH[i-2];
    y[2] = nH[i-1];
    y[3] = nH[i];
    y[4] = nH[i+1];
    FitPolynomial4( x, y, sAMR, &(fConservedQuantities[0]), &error );

    y[1] = ne[i-2];
    y[2] = ne[i-1];
    y[3] = ne[i];
    y[4] = ne[i+1];
    FitPolynomial4( x, y, sAMR, &(fConservedQuantities[1]), &error );

    y[1] = T[i-2];
    y[2] = T[i-1];
    y[3] = T[i];
    y[4] = T[i+1];
    FitPolynomial4( x, y, sAMR, &(fConservedQuantities[2]), &error );

#ifdef ADAPT
    fprintf( pFile, "%.16e\t%.16e\t%.16e\t0.0\t%.16e\t%.16e\t%i", sAMR, ds, (AVERAGE_PARTICLE_MASS*fConservedQuantities[0]),
                                                                  (1.5*BOLTZMANN_CONSTANT*fConservedQuantities[1]*fConservedQuantities[2]),
                                                                  (1.5*BOLTZMANN_CONSTANT*fConservedQuantities[0]*fConservedQuantities[2]),
                                                                  MAX_REFINEMENT_LEVEL );

    for( k=0; k<MAX_REFINEMENT_LEVEL; k++ )
        fprintf( pFile, "\t%ld", ppiID[j][k] );
    fprintf( pFile, "\n" );
    j++;
#else // ADAPT
    fprintf( pFile, "%.16e\t%.16e\t%.16e\t0.0\t%.16e\t%.16e\t0", sAMR, ds, (AVERAGE_PARTICLE_MASS*fConservedQuantities[0]),
                                                                  (1.5*BOLTZMANN_CONSTANT*fConservedQuantities[1]*fConservedQuantities[2]),
                                                                  (1.5*BOLTZMANN_CONSTANT*fConservedQuantities[0]*fConservedQuantities[2]) );

    for( k=0; k<MAX_REFINEMENT_LEVEL; k++ )
        fprintf( pFile, "\t0" );
    fprintf( pFile, "\n" );
#endif // ADAPT

    sAMR += ds;
}
fclose( pFile );

#ifdef ADAPT
// Free the memory allocated to the unique cell ID numbers
for( i=0; i<iMAX_CELLS; i++ )
    free( ppiID[i] );
free( ppiID );
#endif // ADAPT
}

void WritePHYFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params )
{
FILE *pFile;
char szPHYFilename[512];
double v, rho, Pe, PH, P, Cs, Fce, FcH;
int i;

// In the hydrostatic case the velocity equals zero
v = 0.0;

#if defined (OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE)
sprintf( szPHYFilename, "%s.orig.phy", Params.szOutputFilename );
#else // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE
sprintf( szPHYFilename, "%s.phy", Params.szOutputFilename );
#endif // OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE

pFile = fopen( szPHYFilename, "w" );

for( i=0; i<iTotalSteps; i++ )
{
    rho = AVERAGE_PARTICLE_MASS * nH[i];
    Pe = BOLTZMANN_CONSTANT * ne[i] * T[i];
    PH = BOLTZMANN_CONSTANT * nH[i] * T[i];
    P = Pe + PH;
    Cs = sqrt( ( GAMMA * P ) / rho );

    if( i > 1 && i < iTotalSteps-2 )
    {
	Fce = - SPITZER_ELECTRON_CONDUCTIVITY * pow( T[i], 2.5 ) * ( ( T[i+1] - T[i-1] ) / ( s[i+1] - s[i-1] ) );
	FcH = - SPITZER_ION_CONDUCTIVITY * pow( T[i], 2.5 ) * ( ( T[i+1] - T[i-1] ) / ( s[i+1] - s[i-1] ) );
    }
    else
    {
	Fce = 0.0;
	FcH = 0.0;
    }

    fprintf( pFile, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", s[i], v, Cs, ne[i], nH[i], Pe, PH, T[i], T[i], Fce, FcH );
}

fclose( pFile );
}

void WriteSOLFile( double finalH0, PARAMETERS Params )
{
FILE *pFile;
char szSOLFilename[512];

sprintf( szSOLFilename, "%s.sol", Params.szOutputFilename );
pFile = fopen( szSOLFilename, "w" );
fprintf( pFile, "%.16e\n%% Peak heating rate for the solution", finalH0 );
fclose( pFile );
}
