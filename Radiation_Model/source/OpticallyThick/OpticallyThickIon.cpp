// ****
// *
// * Optically-thick Ion Class Function Bodies for Radiative Emission Model
// *
// * Based on the formulation of Carlsson, M., & Leenaarts, J., 2012, A&A, 539, A39
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****

#include "OpticallyThickIon.h"

#ifdef OPTICALLY_THICK_RADIATION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../../Resources/source/file.h"
#include "../../../Resources/source/fitpoly.h"

COpticallyThickIon::COpticallyThickIon( int iZ, char *szIon, char *szAbundFilename )
{
Initialise( iZ, szIon, szAbundFilename );
}

COpticallyThickIon::~COpticallyThickIon( void )
{
FreeAll();
}

void COpticallyThickIon::Initialise( int iZ, char *szIon, char *szAbundFilename )
{
char szEmissFilename[256], szEscProbFilename[256], szIonFracFilename[256], szkappa_0Filename[256];

// Set the atomic number of the element
Z = iZ;

// Construct the filenames of the ion data files
sprintf( szIonFracFilename, "Radiation_Model/atomic_data/OpticallyThick/balances/%s.bal", szIon );
sprintf( szEmissFilename, "Radiation_Model/atomic_data/OpticallyThick/emissivities/%s.em", szIon );
sprintf( szEscProbFilename, "Radiation_Model/atomic_data/OpticallyThick/escape_probabilities/%s.esc", szIon );
sprintf( szkappa_0Filename, "Radiation_Model/atomic_data/OpticallyThick/thermal_conductivities/%s.tc", szIon );

// Get the ion data
GetAbundData( szAbundFilename );
GetIonFracData( szIonFracFilename );
GetEmissData( szEmissFilename );
GetEscProbData( szEscProbFilename );
Getkappa_0Data( szkappa_0Filename );
}

void COpticallyThickIon::FreeAll( void )
{
int i;

// Free the thermal conductivities
if( ikappa_0DP )
{
    for( i=0; i<2; i++ )
        free( ppkappa_0[i] );
    free( ppkappa_0 );
}

// Free the ion escape probabilities
for( i=0; i<2; i++ )
    free( ppEscProb[i] );
free( ppEscProb );

// Free the ion emissivities
for( i=0; i<2; i++ )
    free( ppEmiss[i] );
free( ppEmiss );

for( i=0; i<2; i++ )
    free( ppOriginalEmiss[i] );
free( ppOriginalEmiss );

// Free the ion population fractions
for( i=0; i<2; i++ )
    free( ppIonFrac[i] );
free( ppIonFrac );
}

void COpticallyThickIon::GetAbundData( char *szAbundFilename )
{
FILE *pFile;
double fAb[128], fm_g[128], fAbH, fSumAbund;
double fBuffer;
int iZ[128], iNumElements, iNumMasses;
int i = 0, j, iBuffer;

Abund = 0.0;

// Open the abundance file
pFile = fopen( szAbundFilename, "r" );
// Get the abundance data for hydrogen
fscanf( pFile, "%i", &(iZ[i]) );
ReadDouble( pFile, &fAbH );
fAbH = pow( 10.0, fAbH );
fAb[i] = fAbH;
// Get the element abundances and the abundance of the current element
for(;;)
{
    // Calculate the abundance of the element relative to hydrogen
    fAb[i] /= fAbH;

    // If the atomic number read is the atomic number of the current
    // element then store its abundance relative to hydrogen
    if( iZ[i] == Z )
        Abund = fAb[i];

    i++;

    // Get the atomic number of the element
    fscanf( pFile, "%i", &(iZ[i]) );
    // If the value read is -1 then break
    if( iZ[i] == -1 ) break;

    // Get the log_10 abundance of the element relative to H
    ReadDouble( pFile, &(fAb[i]) );
    fAb[i] = pow( 10.0, fAb[i] );
}
fclose( pFile );

// Calculate the fractional element abundances (e.g. normalised to 1.0)
iNumElements = i;
fSumAbund = 0.0;
for( i=0; i<iNumElements; i++ )
    fSumAbund += fAb[i];

for( i=0; i<iNumElements; i++ )
    fAb[i] /= fSumAbund;

// Open the atomic mass units file
pFile = fopen( AMU_FILENAME, "r" );

// Get the number of masses stored in the file
fscanf( pFile, "%i", &iNumMasses );
j = 0;
for( i=0; i<iNumMasses; i++ )
{
    fscanf( pFile, "%i", &iBuffer );
    // If the atomic number of this element matches the atomic number of the element
    // for which abundance data is available then store the atomic mass number
    if( iBuffer == iZ[j] )
    {
        ReadDouble( pFile, &(fm_g[j]) );
	// Convert from atomic mass units to grams
	fm_g[j] *= AMU;
	j++;
    } else {
        ReadDouble( pFile, &fBuffer );
    }
}
fclose( pFile );

// Calculate the number of hydrogen particles per gram of stellar material for 
// the current set of element abundances
fSumAbund = 0.0;
for( i=0; i<iNumElements; i++ )
    fSumAbund += fAb[i] * fm_g[i];
N_H = fAb[0] * ( 1.0 / fSumAbund );
}

void COpticallyThickIon::GetIonFracData( char *szIonFracFilename )
{
FILE *pFile;
int i;

// Open the ion data file
pFile = fopen( szIonFracFilename, "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iIonFracDP );

// Allocate sufficient memory to hold the ion data
ppIonFrac = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppIonFrac[i] = (double*)malloc( sizeof(double) * iIonFracDP );

for( i=0; i<iIonFracDP; i++ )
{
    // Array index [0][i] contain the log_10 temperatures and [1][i] contain the ion population fractions
    ReadDouble( pFile, &(ppIonFrac[0][i]) );
    ReadDouble( pFile, &(ppIonFrac[1][i]) );
}

fclose( pFile );
}

void COpticallyThickIon::GetEmissData( char *szEmissFilename )
{
FILE *pFile;
int i;

// Open the ion data file
pFile = fopen( szEmissFilename, "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iEmissDP );

// Allocate sufficient memory to hold the ion data
ppOriginalEmiss = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppOriginalEmiss[i] = (double*)malloc( sizeof(double) * iEmissDP );

// Allocate sufficient memory to hold the ion data
ppEmiss = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppEmiss[i] = (double*)malloc( sizeof(double) * iEmissDP );

for( i=0; i<iEmissDP; i++ )
{
    // Array index [0][i] contain the log_10 temperatures and [1][i] contain the log_10 emissivities
    ReadDouble( pFile, &(ppOriginalEmiss[0][i]) );
    ReadDouble( pFile, &(ppOriginalEmiss[1][i]) );

    ppEmiss[0][i] = ppOriginalEmiss[0][i];
    // Convert to the stored values
    ppOriginalEmiss[1][i] = pow( 10.0, ppOriginalEmiss[1][i] );
    ppOriginalEmiss[1][i] *= Abund;
    ppEmiss[1][i] = ppOriginalEmiss[1][i] * GetIonFrac( ppOriginalEmiss[0][i] );
}


fclose( pFile );
}

void COpticallyThickIon::GetEscProbData( char *szEscProbFilename )
{
FILE *pFile;
int i;

// Open the ion data file
pFile = fopen( szEscProbFilename, "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iEscProbDP );

// Allocate sufficient memory to hold the ion data
ppEscProb = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppEscProb[i] = (double*)malloc( sizeof(double) * iEscProbDP );

for( i=0; i<iEscProbDP; i++ )
{
    // Array index [0][i] contain the optical depths or log_10 column densities and [1][i] contain the escape probabilities
    ReadDouble( pFile, &(ppEscProb[0][i]) );
    ReadDouble( pFile, &(ppEscProb[1][i]) );
}

fclose( pFile );
}

void COpticallyThickIon::Getkappa_0Data( char *szkappa_0Filename )
{
FILE *pFile;
int i;

// Open the ion data file
pFile = fopen( szkappa_0Filename, "r" );
if( !pFile )
{
    ikappa_0DP = 0;
    return;
}

// Get the number of data points in the file
fscanf( pFile, "%i", &ikappa_0DP );

// Allocate sufficient memory to hold the ion data
ppkappa_0 = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppkappa_0[i] = (double*)malloc( sizeof(double) * ikappa_0DP );

for( i=0; i<ikappa_0DP; i++ )
{
    // Array index [0][i] contain the log_10 temperatures and [1][i] contain the log_10 thermal conductivities
    ReadDouble( pFile, &(ppkappa_0[0][i]) );
    ReadDouble( pFile, &(ppkappa_0[1][i]) );
}

fclose( pFile );
}

double COpticallyThickIon::GetIonFrac( double flog_10T )
{
int i;
double x[3], y[3], fIonFrac;

for( i=0; i<iIonFracDP; i++ )
{
    if( flog_10T <= ppIonFrac[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iIonFracDP ) i = iIonFracDP - 1;

x[1] = ppIonFrac[0][i-1];
x[2] = ppIonFrac[0][i];
y[1] = ppIonFrac[1][i-1];
y[2] = ppIonFrac[1][i];

LinearFit( x, y, flog_10T, &fIonFrac );

// Trap limits
if( fIonFrac < 0.0 ) fIonFrac = 0.0;
else if( fIonFrac > 1.0 ) fIonFrac = 1.0;

return fIonFrac;
}

double COpticallyThickIon::GetEmiss( double flog_10T )
{
int i;
double x[3], y[3], fEmiss;

// Trap limits
if( flog_10T < ppEmiss[0][0] ) return ppEmiss[1][0];
else if( flog_10T > ppEmiss[0][iEmissDP-1] ) return ppEmiss[1][iEmissDP-1];

for( i=0; i<iEmissDP; i++ )
{
    if( flog_10T <= ppEmiss[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iEmissDP ) i = iEmissDP - 1;

x[1] = ppEmiss[0][i-1];
x[2] = ppEmiss[0][i];
y[1] = ppEmiss[1][i-1];
y[2] = ppEmiss[1][i];

LinearFit( x, y, flog_10T, &fEmiss );

// Trap limits
if( fEmiss < 0.0 ) fEmiss = 0.0;

return fEmiss;
}

double COpticallyThickIon::GetEmiss( double flog_10T, double fIonFrac )
{
int i;
double x[3], y[3], fEmiss;

// Trap limits
if( flog_10T < ppOriginalEmiss[0][0] ) return ppOriginalEmiss[1][0] * fIonFrac;
else if( flog_10T > ppOriginalEmiss[0][iEmissDP-1] ) return ppOriginalEmiss[1][iEmissDP-1] * fIonFrac;

for( i=0; i<iEmissDP; i++ )
{
    if( flog_10T <= ppOriginalEmiss[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iEmissDP ) i = iEmissDP - 1;

x[1] = ppOriginalEmiss[0][i-1];
x[2] = ppOriginalEmiss[0][i];
y[1] = ppOriginalEmiss[1][i-1];
y[2] = ppOriginalEmiss[1][i];

LinearFit( x, y, flog_10T, &fEmiss );

// Trap limits
if( fEmiss < 0.0 ) fEmiss = 0.0;

return fEmiss * fIonFrac;
}

double COpticallyThickIon::GetEscProb( double fX )
{
int i;
double x[3], y[3], fEscProb;

for( i=0; i<iEscProbDP; i++ )
{
    if( fX <= ppEscProb[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iEscProbDP ) i = iEscProbDP - 1;

x[1] = ppEscProb[0][i-1];
x[2] = ppEscProb[0][i];
y[1] = ppEscProb[1][i-1];
y[2] = ppEscProb[1][i];

LinearFit( x, y, fX, &fEscProb );

// Trap limits
if( fEscProb < 0.0 ) fEscProb = 0.0;
else if( fEscProb > 1.0 ) fEscProb = 1.0;

return fEscProb;
}

double COpticallyThickIon::Getkappa_0( double flog_10T )
{
int i;
double x[3], y[3], fkappa_0;

if( !ikappa_0DP ) return 0.0;

// Trap limits
// Commented out to allow extrapolation to temperatures lower than 10^3.70 K
/*
if( flog_10T < ppkappa_0[0][0] ) return pow( 10.0, ppkappa_0[1][0] );
else if( flog_10T > ppkappa_0[0][ikappa_0DP-1] ) return pow( 10.0, ppkappa_0[1][ikappa_0DP-1] );
*/

for( i=0; i<ikappa_0DP; i++ )
{
    if( flog_10T <= ppkappa_0[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == ikappa_0DP ) i = ikappa_0DP - 1;

x[1] = ppkappa_0[0][i-1];
x[2] = ppkappa_0[0][i];
y[1] = ppkappa_0[1][i-1];
y[2] = ppkappa_0[1][i];

LinearFit( x, y, flog_10T, &fkappa_0 );

return pow( 10.0, fkappa_0 );
}

double COpticallyThickIon::GetVolumetricLossRate( double flog_10T, double fX, double n_e_rho )
{
return GetEmiss( flog_10T ) * GetEscProb( fX ) * N_H * n_e_rho;
}

double COpticallyThickIon::GetVolumetricLossRate( double flog_10T, double fIonFrac, double fX, double n_e_rho )
{
return GetEmiss( flog_10T, fIonFrac ) * GetEscProb( fX ) * N_H * n_e_rho;
}

#endif // OPTICALLY_THICK_RADIATION