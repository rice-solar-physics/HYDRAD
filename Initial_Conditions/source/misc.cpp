// ****
// *
// * Miscellaneous routines used by the hydrostatic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 28/10/2013
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

#ifdef USE_TABULATED_GRAVITY
#else // USE_TABULATED_GRAVITY
void GenerateSemiCircularLoop( PARAMETERS Params )
{
double CalcSolarGravity( double s, double Lfull, double Inc );

FILE *pFile;
char szGravityFilename[256];

double s, ds;
int i, igdp;

igdp = MIN_CELLS / 5; // Hard wired

s = 0.0;
ds = Params.Lfull / ((double)igdp);

sprintf( szGravityFilename, "%s.gravity", Params.szOutputFilename );
pFile = fopen( szGravityFilename, "w" );
fprintf( pFile, "%i\n", igdp+1 );
for( i=0; i<=igdp; i++ )
{
    fprintf( pFile, "%.16e\t%.16e\n", s, CalcSolarGravity( s, Params.Lfull, Params.Inc ) );
    s += ds;
}
fclose( pFile );
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
#endif // USE_TABULATED_GRAVITY

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

pFile = fopen( Params.szOutputFilename, "w" );
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
char szPHYFilename[256];
int i;

sprintf( szPHYFilename, "%s.phy", Params.szOutputFilename );
pFile = fopen( szPHYFilename, "w" );
for( i=0; i<iTotalSteps; i++ )
    fprintf( pFile, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", s[i], ne[i], nH[i], T[i], BOLTZMANN_CONSTANT * ne[i] * T[i], BOLTZMANN_CONSTANT * nH[i] * T[i] );
fclose( pFile );
}

void WriteSOLFile( double finalH0, PARAMETERS Params )
{
FILE *pFile;
char szSOLFilename[256];

sprintf( szSOLFilename, "%s.sol", Params.szOutputFilename );
pFile = fopen( szSOLFilename, "w" );
fprintf( pFile, "%.16e\n%% Peak heating rate for the solution", finalH0 );
fclose( pFile );
}
