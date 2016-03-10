// ****
// *
// * A heating code to simulate different forms of heat deposition
// *
// * Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *     
// * Date last modified: 07/12/2015
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "heat.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"


CHeat::CHeat( char *szFilename, double fL )
{
Initialise( szFilename, fL );
}

CHeat::~CHeat( void )
{
FreeAll();
}

void CHeat::Initialise( char *szFilename, double fL )
{
// Get the length of the loop
fLoopLength = fL;

GetHeatingData( szFilename );
GetVALHeatingData();
}

void CHeat::FreeAll( void )
{
int i;

for( i=0; i<2; i++ )
    free( ppVALHeating[i] );
free( ppVALHeating );

if( !NumActivatedEvents )
    return;

// Location, scale-length and maximum energy
free( s0episodic );
free( sHepisodic );
free( E0episodic );

// Start and end times of rise phase
free( tsRepisodic );
free( teRepisodic );

// Start and end times of decay phase
free( tsDepisodic );
free( teDepisodic );
}

void CHeat::GetHeatingData( char *szFilename )
{
FILE *pConfigFile;
int i;

// Open and read the configuration file
pConfigFile = fopen( szFilename, "r" );

// Get the duration of the heating simulation
ReadDouble( pConfigFile, &fDuration );

// Get the quiescent heating parameter values
ReadDouble( pConfigFile, &s0quiescent );
ReadDouble( pConfigFile, &sHquiescent );
ReadDouble( pConfigFile, &E0quiescent );

// Get the episodic heating event parameter values
fscanf( pConfigFile, "%i", &NumActivatedEvents );

if( !NumActivatedEvents )
{
    // Close the configuration file
    fclose( pConfigFile );
    return;
}

// Allocate sufficient memory to store the positioning and timing information for each event

// Location, scale-length and maximum energy
s0episodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
sHepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
E0episodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

// Start and end times of rise phase
tsRepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
teRepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

// Start and end times of decay phase
tsDepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
teDepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

for( i=0; i<NumActivatedEvents; i++ )
{
    // Location, scale-length and maximum energy
    ReadDouble( pConfigFile, &(s0episodic[i]) );
    ReadDouble( pConfigFile, &(sHepisodic[i]) );
    ReadDouble( pConfigFile, &(E0episodic[i]) );

    // Start and end times of rise phase
    ReadDouble( pConfigFile, &(tsRepisodic[i]) );
    ReadDouble( pConfigFile, &(teRepisodic[i]) );

    // Start and end times of decay phase
    ReadDouble( pConfigFile, &(tsDepisodic[i]) );
    ReadDouble( pConfigFile, &(teDepisodic[i]) );
}

// Close the configuration file
fclose( pConfigFile );
}

void CHeat::GetVALHeatingData( void )
{
FILE *pFile;
int i;

// Open and read the configuration file
pFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat", "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iVALHeatingDP );

// Allocate sufficient memory to hold the heating data
ppVALHeating = (double**)malloc( sizeof(double) * 2 );
for( i=0; i<2; i++ )
    ppVALHeating[i] = (double*)malloc( sizeof(double) * iVALHeatingDP );

for( i=0; i<iVALHeatingDP; i++ )
{
    // Array index [0][i] contain the column mass density and [1][i] contain the volumetric heating rate
    ReadDouble( pFile, &(ppVALHeating[0][i]) );
    ReadDouble( pFile, &(ppVALHeating[1][i]) );
}

fclose( pFile );
}

double CHeat::CalculateQuiescentHeating( double s )
{
double fQuiescent = 0.0, term1;

// Left footpoint
if( E0quiescent )
{
    term1 = ( s - s0quiescent );
    term1 *= term1;

    fQuiescent = E0quiescent * exp( - term1 / (2.0*sHquiescent*sHquiescent) );
}

return fQuiescent;
}

double CHeat::CalculateEpisodicHeating( double s, double t )
{
double fEpisodic = 0.0, term1, term2;
int i;

for( i=0; i<NumActivatedEvents; i++ )
{
    if( t >= tsRepisodic[i] && t <= teDepisodic[i] )
    {
        if( t >= tsRepisodic[i] && t <= teRepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
            term1 *= term1;

            term2 = ( t - tsRepisodic[i] ) / ( teRepisodic[i] - tsRepisodic[i] );

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) ) * term2;
        }
        else if( t > teRepisodic[i] && t < tsDepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
        	term1 *= term1;

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) );
        }
        else if( t >= tsDepisodic[i] && t <= teDepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
            term1 *= term1;

            term2 = 1.0 - ( ( t - tsDepisodic[i] ) / ( teDepisodic[i] - tsDepisodic[i] ) );

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) ) * term2;
        }
    }
}

return fEpisodic;
}

double CHeat::CalculateHeating( double s, double t )
{
double fQuiescentHeating = 0.0, fEpisodicHeating = 0.0;

if( E0quiescent )
    fQuiescentHeating = CalculateQuiescentHeating( s );

if( NumActivatedEvents )
    fEpisodicHeating = CalculateEpisodicHeating( s, t );

return( fQuiescentHeating + fEpisodicHeating );
}

double CHeat::CalculateVALHeating( double flog10_rho_c )
{
int i;
double x[3], y[3], fVALHeating;

if( !iVALHeatingDP ) return 0.0;

// Trap limits
if( flog10_rho_c < ppVALHeating[0][0] ) return pow( 10.0, ppVALHeating[1][0] );
else if( flog10_rho_c > ppVALHeating[0][iVALHeatingDP-1] ) return pow( 10.0, ppVALHeating[1][iVALHeatingDP-1] );

for( i=0; i<iVALHeatingDP; i++ )
{
    if( flog10_rho_c <= ppVALHeating[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iVALHeatingDP ) i = iVALHeatingDP - 1;

x[1] = ppVALHeating[0][i-1];
x[2] = ppVALHeating[0][i];
y[1] = ppVALHeating[1][i-1];
y[2] = ppVALHeating[1][i];

LinearFit( x, y, flog10_rho_c, &fVALHeating );

return pow( 10.0, fVALHeating );
}
