// ****
// *
// * Element Class Function Bodies for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "element.h"
#include "config.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"


CElement::CElement( int iZ, char *szRangesFilename, char *szAbundFilename, char *szEmissFilename, char *szRatesFilename, char *szIonFracFilename )
{
Initialise( iZ, szRangesFilename, szAbundFilename, szEmissFilename, szRatesFilename, szIonFracFilename );
}

CElement::~CElement( void )
{
FreeAll();
}

void CElement::Initialise( int iZ, char *szRangesFilename, char *szAbundFilename, char *szEmissFilename, char *szRatesFilename, char *szIonFracFilename )
{
// Set the atomic number of the element
Z = iZ;

// Open the data files and initialise the element
OpenRangesFile( szRangesFilename );
OpenAbundanceFile( szAbundFilename );
OpenEmissivityFile( szEmissFilename );
OpenRatesFile( szRatesFilename );
OpenIonFracFile( szIonFracFilename );

// Calculate phi for each ion as a function of temperature and density
CalculatePhi();

// Calculate the total phi of all radiating elements as a function of temperature and density
CalculateTotalPhi();
}

void CElement::OpenRangesFile( char *szRangesFilename )
{
FILE *pFile;
int i;
char buffer;

// Open the ranges file
pFile = fopen( szRangesFilename, "r" );

// Read the comment line
buffer = 0;
while( buffer != '.' )
    fscanf( pFile, "%c", &buffer );

// Get the number of temperature and density values
fscanf( pFile, "%i", &NumTemp );
fscanf( pFile, "%i", &NumDen );

// Allocate arrays to hold the log_10 temperature and density values
pTemp = (double*)malloc( sizeof(double) * NumTemp );
pDen = (double*)malloc( sizeof(double) * NumDen );

// Read the comment line
buffer = 0;
while( buffer != '.' )
    fscanf( pFile, "%c", &buffer );

// Get the log_10 temperature values
for( i=0; i<NumTemp; i++ )
    ReadDouble( pFile, &(pTemp[i]) );

// Read the comment line
buffer = 0;
while( buffer != '.' )
    fscanf( pFile, "%c", &buffer );

// Get the log_10 density values
for( i=0; i<NumDen; i++ )
    ReadDouble( pFile, &(pDen[i]) );

fclose( pFile );
}

void CElement::OpenAbundanceFile( char *szAbundFilename )
{
FILE *pFile;
double fHAb;
double fTemp;
int buffer;

// Open the abundance file
pFile = fopen( szAbundFilename, "r" );

// Get the abundance data for hydrogen
fscanf( pFile, "%i", &buffer );
ReadDouble( pFile, &fTemp );
fHAb = pow( 10.0, fTemp );

// Search for the abundance of the required element

for(;;)
{
    // If the atomic number read is the atomic number of the element then calculate the abundance relative to hydrogen
    if( buffer == Z )
    {
        fAbund = pow( 10.0, fTemp ) / fHAb;
	break;
    }

    // Get the atomic number of the element
    fscanf( pFile, "%i", &buffer );

    // If the value read is -1 then break
    if( buffer == -1 ) break;

    // Get the log_10 abundance of the element relative to H
    ReadDouble( pFile, &fTemp );
}

// If the atomic number of the element was not found in the abundance file then set the abundance to zero
if( buffer == -1 )
    fAbund = 0.0;

fclose( pFile );
}

void CElement::OpenEmissivityFile( char *szEmissFilename )
{
FILE *pFile;
double fTemp;
int i, j, NumTempxNumDen, indexTemp, indexDen;
char buffer;

// Open the emissivity file
pFile = fopen( szEmissFilename, "r" );

// Get the ion data

// Read the comment line
buffer = 0;
while( buffer != '.' )
    fscanf( pFile, "%c", &buffer );
	
fscanf( pFile, "%i", &NumIons );

// Allocate an array to hold the ion list
pSpecNum = (int*)malloc( sizeof(int) * NumIons );

// Read the comment line
buffer = 0;
while( buffer != '.' )
    fscanf( pFile, "%c", &buffer );

// Get the spectroscopic numbers of each ion
for( i=0; i<NumIons; i++ )
    fscanf( pFile, "%i", &pSpecNum[i] );

// Get the emissivity values for each ion

// Allocate an array to hold the pointers to the emissivity for each ion
ppEmiss = (double**)malloc( sizeof(double*) * NumIons );

// Calculate the 2D array sizes
NumTempxNumDen = NumTemp * NumDen;

// Allocate arrays to hold the NumTemp * NumDen emissivity values for each ion and get the values from the file
for( i=0; i<NumIons; i++)
{
    ppEmiss[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );

    // Read the comment line
    buffer = 0;
    while( buffer != '.' )
        fscanf( pFile, "%c", &buffer );

    indexTemp = 0;
    indexDen = 0;
		
    for( j=0; j<NumTempxNumDen; j++ )
    {
        ReadDouble( pFile, &fTemp );
    
	// The value stored is the product of the Chianti calculated emissivity obtained using emiss_calc, the constant 0.83 and the abundance relative to hydrogen of the element, divided by the electron number density
	ppEmiss[i][j] = ( 0.83 * fAbund * fTemp ) / pow( 10.0, pDen[indexDen] );

	indexTemp++;
		
	if( indexTemp == NumTemp )
	{
            indexTemp = 0;
            indexDen++;
	}
    }
}

fclose( pFile );
}

void CElement::OpenRatesFile( char *szRatesFilename )
{
FILE *pFile;
int i, j, NumTempxNumDen;
char buffer[8];

// Allocate arrays to hold the pointers to the rates for each ion
ppIonRate = (double**)malloc( sizeof(double*) * Z );
ppRecRate = (double**)malloc( sizeof(double*) * Z );

#ifdef DENSITY_DEPENDENT_RATES
    // Calculate the 2D array sizes
    NumTempxNumDen = NumTemp * NumDen;
#else // DENSITY_DEPENDENT_RATES
    NumTempxNumDen = NumTemp; // NumDen = 1
#endif // DENSITY_DEPENDENT_RATES

// Open the rates file
pFile = fopen( szRatesFilename, "r" );

// Get the ionisation rates

for( i=0; i<Z; i++ )
{
    // Read the rates into memory
		
    ppIonRate[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );
			
    fscanf( pFile, "%s", buffer );
		
    for( j=0; j<NumTempxNumDen; j++ )
        ReadDouble( pFile, &(ppIonRate[i][j]) );
}

// Get the recombination rates

for( i=0; i<Z; i++ )
{
    // Read the rates into memory
	
    ppRecRate[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );

    fscanf( pFile, "%s", buffer );
		
    for( j=0; j<NumTempxNumDen; j++ )
        ReadDouble( pFile, &(ppRecRate[i][j]) );
}

fclose( pFile );
}

void CElement::OpenIonFracFile( char *szIonFracFilename )
{
FILE *pFile;
double fTemp;
int i, j, NumTempxNumDen;

// Allocate arrays to hold the pointers to the fractional populations for each ion
ppIonFrac = (double**)malloc( sizeof(double*) * ( Z + 1 ) );

#ifdef DENSITY_DEPENDENT_RATES
    // Calculate the 2D array sizes
    NumTempxNumDen = NumTemp * NumDen;
#else // DENSITY_DEPENDENT_RATES
    NumTempxNumDen = NumTemp; // NumDen = 1
#endif // DENSITY_DEPENDENT_RATES

// Allocate an array for each ion to contain the fractional population of that ion as a function of temperature
for( i=0; i<=Z; i++ )
    ppIonFrac[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );

pFile = fopen( szIonFracFilename, "r" );

// Get the fractional populations

for( j=0; j<NumTempxNumDen; j++ )
{
    // Read the log10 temperature value
    ReadDouble( pFile, &fTemp );

    for( i=0; i<=Z; i++ )
        ReadDouble( pFile, &(ppIonFrac[i][j]) );
}

fclose( pFile );
}

void CElement::CalculatePhi( void )
{
int i, j, NumTempxNumDen, indexTemp, indexDen;

// Allocate an array to hold the pointers to the values of phi( n, T ) for the element
ppPhi = (double**)malloc( sizeof(double*) * NumIons );

// Calculate the 2D array sizes
NumTempxNumDen = NumTemp * NumDen;

// Allocate arrays to hold the NumTemp * NumDen values of phi for each ion
for( i=0; i<NumIons; i++)
{
    ppPhi[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );

    indexTemp = 0;
    indexDen = 0;

    for( j=0; j<NumTempxNumDen; j++ )
    {
        // Calculate the factor phi( n, T ) for the current ion at the current temperature and density
#ifdef DENSITY_DEPENDENT_RATES
	ppPhi[i][j] = GetIonEmissivity( pSpecNum[i], pTemp[indexTemp], pDen[indexDen] ) * GetEquilIonFrac( pSpecNum[i], pTemp[indexTemp], pDen[indexDen] );
#else // DENSITY_DEPENDENT_RATES
	ppPhi[i][j] = GetIonEmissivity( pSpecNum[i], pTemp[indexTemp], pDen[indexDen] ) * GetEquilIonFrac( pSpecNum[i], pTemp[indexTemp] );
#endif // DENSITY_DEPENDENT_RATES
        indexTemp++;
                
        if( indexTemp == NumTemp )
        {
            indexTemp = 0;
           indexDen++;
        }
    }
}
}

void CElement::CalculateTotalPhi( void )
{
int i, j, NumTempxNumDen;

// Calculate the 2D array size
NumTempxNumDen = NumTemp * NumDen;

// Allocate an array to hold the total values of phi( n, T ) for the element
pTotalPhi = (double*)malloc( sizeof(double) * NumTempxNumDen );

for( j=0; j<NumTempxNumDen; j++ )
{
    pTotalPhi[j] = 0.0;

    // Calculate the total equilibrium radiation for the element at the current temperature and density
    for( i=0; i<NumIons; i++)
        pTotalPhi[j] += ppPhi[i][j];
}
}

void CElement::FreeAll( void )
{
int i;

free( pSpecNum );
free( pTemp );
free( pDen );
free( pTotalPhi );

for( i=0; i<NumIons; i++)
{
    free( ppEmiss[i] );
    free( ppPhi[i] );
}

for( i=0; i<Z; i++ )
{
    free( ppIonRate[i] );
    free( ppRecRate[i] );
}

for( i=0; i<=Z; i++ )
    free( ppIonFrac[i] );

free( ppEmiss );
free( ppIonRate );
free( ppRecRate );
free( ppIonFrac );
free( ppPhi );
}

double CElement::GetIonEmissivity( int iIon, double flog_10T, double flog_10n )
{
double x1[5], x2[5], **y, *pfTemp, result, error;
int i, j, k, l;

// Select the required ion

for( i=0; i<NumIons; i++ )
    if( iIon == pSpecNum[i] ) break;

// If the requested ion is not in the list of spectroscopic numbers then return 0.0
if( i == NumIons ) return 0.0;

// Select the four temperature values surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then set it to the appropriate limit
if( flog_10n < pDen[0] )
    flog_10n = pDen[0];
else if ( flog_10n > pDen[NumDen-1] )
    flog_10n = pDen[NumDen-1];

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog_10n ) break;

// Deal with the special cases where there aren't two values either side of the
// desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the i'th ion
// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen emissivity values corresponding to the grid

// Allocate an array of pointers to pointers for the emissivity values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the emissivity values corresponding to the set of temperatures at the current density
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the emissivity values for the required ion
    pfTemp = ppEmiss[i];

    // Skip to the emissivity set corresponding to the l'th density value and get the emissivities corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, &result, &error );

// Check emissivity is physically realistic
if( result < 0.0 ) result = 0.0;

// Return the emissivity for the required ion at the specified temperature and density
// The hydrocode must supply the electron density and the ionisation balance
return result;
}

double CElement::GetAbundance( void )
{
return fAbund;
}

void CElement::GetRates( int iIon, double flog_10T, double *pfIonRate, double *pfRecRate )
{
double x[5], Ion_y[5], Rec_y[5], error;
int i, j;

if( !iIon || iIon > Z )
{
    *pfIonRate = 0.0;
    *pfRecRate = 0.0;

    return;
}

// Select the required ion
i = iIon - 1;

// Select the four temperature values, ionisation and recombination rates surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x[1] = pTemp[j-2];
x[2] = pTemp[j-1];
x[3] = pTemp[j];
x[4] = pTemp[j+1];

Ion_y[1] = ppIonRate[i][j-2];
Ion_y[2] = ppIonRate[i][j-1];
Ion_y[3] = ppIonRate[i][j];
Ion_y[4] = ppIonRate[i][j+1];

Rec_y[1] = ppRecRate[i][j-2];
Rec_y[2] = ppRecRate[i][j-1];
Rec_y[3] = ppRecRate[i][j];
Rec_y[4] = ppRecRate[i][j+1];

// Perform the polynomial interpolation
FitPolynomial( x, Ion_y, 4, flog_10T, pfIonRate, &error );
FitPolynomial( x, Rec_y, 4, flog_10T, pfRecRate, &error );

// Check rates are physically realistic
if( *pfIonRate < 0.0 ) *pfIonRate = 0.0;
if( *pfRecRate < 0.0 ) *pfRecRate = 0.0;
}

void CElement::GetRates( int iIon, double flog_10T, double flog_10n, double *pfIonRate, double *pfRecRate )
{
double x1[5], x2[5], **y, *pfTemp, error;
int i, j, k, l;

if( !iIon || iIon > Z )
{
    *pfIonRate = 0.0;
    *pfRecRate = 0.0;

    return;
}

// Select the required ion
i = iIon - 1;

// Select the four temperature values, ionisation and recombination rates surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then set it to the appropriate limit
if( flog_10n < pDen[0] )
    flog_10n = pDen[0];
else if ( flog_10n > pDen[NumDen-1] )
    flog_10n = pDen[NumDen-1];

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog_10n ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the i'th ion
// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen ionisation and recombination rates corresponding to the grid

// Allocate an array of pointers to pointers for the ionisation and recombination rates
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the ionisation and recombination values corresponding to the set of temperatures at the current density
for( l=1; l<=4; l++ )
	y[l] = (double*)alloca( sizeof(double) * 5 );

// Ionisation
for( l=1; l<=4; l++ )
{
    // Point to the start of the ionisation rates for the required ion
    pfTemp = ppIonRate[i];

    // Skip to the ionisation rate set corresponding to the l'th density value and get the ionisation rates corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}
// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, pfIonRate, &error );

// Recombination
for( l=1; l<=4; l++ )
{
    // Point to the start of the recombination rates for the required ion
    pfTemp = ppRecRate[i];

    // Skip to the recombination rate set corresponding to the l'th density value and get the recombination rates corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}
// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, pfRecRate, &error );

// Check rates are physically realistic
if( *pfIonRate < 0.0 ) *pfIonRate = 0.0;
if( *pfRecRate < 0.0 ) *pfRecRate = 0.0;
}

double CElement::GetEquilIonFrac( int iIon, double flog_10T )
{
double x[5], y[5], IonFrac, error;
int i, j;

if( !iIon || iIon > Z+1 )
    return 0.0;

// Select the required ion
i = iIon - 1;

// Select the four temperature values and fractional populations surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x[1] = pTemp[j-2];
x[2] = pTemp[j-1];
x[3] = pTemp[j];
x[4] = pTemp[j+1];

y[1] = ppIonFrac[i][j-2];
y[2] = ppIonFrac[i][j-1];
y[3] = ppIonFrac[i][j];
y[4] = ppIonFrac[i][j+1];

// Perform the polynomial interpolation
FitPolynomial( x, y, 4, flog_10T, &IonFrac, &error );

// Ensure the minimum ion fraction remains above the cut-off and is physically realistic
if( IonFrac < CUTOFF_ION_FRACTION )
    IonFrac = CUTOFF_ION_FRACTION;

return IonFrac;
}

double CElement::GetEquilIonFrac( int iIon, double flog_10T, double flog_10n )
{
double x1[5], x2[5], **y, *pfTemp, IonFrac, error;
int i, j, k, l;

if( !iIon || iIon > Z+1 )
    return 0.0;

// Select the required ion
i = iIon - 1;

// Select the four temperature values and fractional populations surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then set it to the appropriate limit
if( flog_10n < pDen[0] )
    flog_10n = pDen[0];
else if ( flog_10n > pDen[NumDen-1] )
    flog_10n = pDen[NumDen-1];

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog_10n ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the i'th ion
// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen ion population fractions corresponding to the grid

// Allocate an array of pointers to pointers for the ion population fractions
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the ion population fractions corresponding to the set of temperatures at the current density
for( l=1; l<=4; l++ )
	y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the ion population fractions for the required ion
    pfTemp = ppIonFrac[i];

    // Skip to the ion population fraction set corresponding to the l'th density value and get the ion population fractions corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}
// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, &IonFrac, &error );

// Ensure the minimum ion fraction remains above the cut-off and is physically realistic
if( IonFrac < CUTOFF_ION_FRACTION )
    IonFrac = CUTOFF_ION_FRACTION;

return IonFrac;
}

double CElement::GetEmissivity( int iIon, double flog_10T, double flog_10n )
{
double x1[5], x2[5], **y, *pfTemp, result, error;
int i, j, k, l;

// Select the required ion

for( i=0; i<NumIons; i++ )
    if( iIon == pSpecNum[i] ) break;

// If the requested ion is not in the list of spectroscopic numbers then return 0.0
if( i == NumIons ) return 0.0;

// Select the four temperature values surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then set it to the appropriate limit
if( flog_10n < pDen[0] )
    flog_10n = pDen[0];
else if ( flog_10n > pDen[NumDen-1] )
    flog_10n = pDen[NumDen-1];

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog_10n ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the i'th ion
// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen phi( n, T ) values corresponding to the grid

// Allocate an array of pointers to pointers for the values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the values corresponding to the specified range of densities and temperatures
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the phi( n, T ) values for the required ion
    pfTemp = ppPhi[i];

    // Skip to the set corresponding to the l'th density value and get the
    // values of phi( n, T ) corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, &result, &error );

// Check the value of phi( n, T ) is physically realistic
if( result < 0.0 ) result = 0.0;

// Return the equilibrium emissivity for the required ion at the specified temperature and density
return result;
}

double CElement::GetEmissivity( double flog_10T, double flog_10n )
{
double x1[5], x2[5], **y, *pfTemp, result, error;
int j, k, l;

// Select the four temperature values surrounding the desired one

// If the temperature is out of range then set it to the appropriate limit
if( flog_10T < pTemp[0] )
    flog_10T = pTemp[0];
else if( flog_10T > pTemp[NumTemp-1] )
    flog_10T = pTemp[NumTemp-1];

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog_10T ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then set it to the appropriate limit
if( flog_10n < pDen[0] )
    flog_10n = pDen[0];
else if ( flog_10n > pDen[NumDen-1] )
    flog_10n = pDen[NumDen-1];

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog_10n ) break;

// Deal with the special cases where there aren't two values either side of the desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen phi( n, T ) values corresponding to the grid

// Allocate an array of pointers to pointers for the values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the values corresponding to the specified range of densities and temperatures
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the phi( n, T ) values for the element
    pfTemp = pTotalPhi;

    // Skip to the set corresponding to the l'th density value and get the values of total phi( n, T ) corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, &result, &error );

// Check the value of phi( n, T ) is physically realistic
if( result < 0.0 ) result = 0.0;

// Return the total equilibrium emissivity at the specified temperature and density
return result;
}

void CElement::Getdnibydt( double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double delta_s, double *pdnibydt, double *pTimeScale )
{
double ne, IonRate[2], RecRate[2], term1, term2, term3, term4, term5, delta_t1, delta_t2, TimeScale, SmallestTimeScale;
int iIndex, iSpecNum;

// Variables used for interpolation
double x[3], y[3];

// Variables used for transport (flux) calculation
double Q1, Q2, Q3, QT, ni0, ni2;

// Initialise the characteristic time-scales
TimeScale = SmallestTimeScale = LARGEST_DOUBLE;

// Calculate the electron number density
ne = pow( 10.0, flog_10n );

for( iIndex=0; iIndex<=Z; iIndex++ )
{
    // Reset the rates
    IonRate[0] = IonRate[1] = 0.0;
    RecRate[0] = RecRate[1] = 0.0;

    iSpecNum = iIndex + 1;

    // Calculate the fluxes to be used with Barton's Method

    if( pv[0] > 0.0 )
    {
        // Calculate the ion fraction
	    
        x[1] = s[0];
	x[2] = s[1];
	y[1] = pni0[iIndex];
	y[2] = pni1[iIndex];
	LinearFit( x, y, s_pos[0], &Q1 );

	x[1] = s[1];
	x[2] = s[2];
	y[1] = pni1[iIndex];
	y[2] = pni2[iIndex];
	LinearFit( x, y, s_pos[0], &Q2 );

	Q3 = pni1[iIndex];

	if( pni2[iIndex] <= pni1[iIndex] )
	{
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                ni0 = Q3;
            else
                ni0 = QT;
	}
	else
	{
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                ni0 = Q3;
            else
                ni0 = QT;
	}
    }
    else
    {
        // Calculate the ion fraction
	        
        x[1] = s[2];
        x[2] = s[3];
        y[1] = pni2[iIndex];
        y[2] = pni3[iIndex];
        LinearFit( x, y, s_pos[0], &Q1 );

        x[1] = s[1];
        x[2] = s[2];
        y[1] = pni1[iIndex];
        y[2] = pni2[iIndex];
        LinearFit( x, y, s_pos[0], &Q2 );

        Q3 = pni2[iIndex];

        if( pni2[iIndex] <= pni1[iIndex] )
        {
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                ni0 = Q3;
            else
                ni0 = QT;
        }
        else
        {
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                ni0 = Q3;
            else
                ni0 = QT;
        }
    }

    if( pv[2] > 0.0 )
    {
        // Calculate the ion fraction
	    
	x[1] = s[1];
	x[2] = s[2];
	y[1] = pni1[iIndex];
	y[2] = pni2[iIndex];
	LinearFit( x, y, s_pos[2], &Q1 );

	x[1] = s[2];
	x[2] = s[3];
	y[1] = pni2[iIndex];
	y[2] = pni3[iIndex];
	LinearFit( x, y, s_pos[2], &Q2 );

	Q3 = pni2[iIndex];

	if( pni3[iIndex] <= pni2[iIndex] )
	{
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                ni2 = Q3;
            else
                ni2 = QT;
	}
	else
	{
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                ni2 = Q3;
            else
                ni2 = QT;
	}
    }
    else
    {
        // Calculate the ion fraction
	        
	x[1] = s[3];
        x[2] = s[4];
        y[1] = pni3[iIndex];
        y[2] = pni4[iIndex];
        LinearFit( x, y, s_pos[2], &Q1 );

        x[1] = s[2];
        x[2] = s[3];
        y[1] = pni2[iIndex];
        y[2] = pni3[iIndex];
        LinearFit( x, y, s_pos[2], &Q2 );

        Q3 = pni3[iIndex];

        if( pni3[iIndex] <= pni2[iIndex] )
        {
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                ni2 = Q3;
            else
                ni2 = QT;
        }
        else
        {
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                ni2 = Q3;
            else
                ni2 = QT;
        }
    }

    term1 = - ( ( ni2 * pv[2] ) - ( ni0 * pv[0] ) ) / delta_s;

    if( iSpecNum > 1 )
    {
#ifdef DENSITY_DEPENDENT_RATES
        GetRates( iSpecNum-1, flog_10T, flog_10n, &IonRate[0], &RecRate[0] );
#else // DENSITY_DEPENDENT_RATES
        GetRates( iSpecNum-1, flog_10T, &IonRate[0], &RecRate[0] );
#endif // DENSITY_DEPENDENT_RATES
	term2 = pni2[iIndex-1] * IonRate[0];
    }
    else
        term2 = 0.0;
	
    if( iSpecNum < Z+1 )
    {
#ifdef DENSITY_DEPENDENT_RATES
        GetRates( iSpecNum, flog_10T, flog_10n, &IonRate[1], &RecRate[1] );
#else // DENSITY_DEPENDENT_RATES
        GetRates( iSpecNum, flog_10T, &IonRate[1], &RecRate[1] );
#endif // DENSITY_DEPENDENT_RATES
	term3 = pni2[iIndex+1] * RecRate[1];
    }
    else
        term3 = 0.0;

    term4 = - pni2[iIndex] * ( IonRate[1] + RecRate[0] );
    term5 = ne * ( term2 + term3 + term4 );

    if( term5 && pni2[iIndex] >= CUTOFF_ION_FRACTION )
    {
        delta_t1 = SAFETY_ATOMIC * ( EPSILON_D / fabs(term5) );
        delta_t2 = SAFETY_ATOMIC * pni2[iIndex] * ( EPSILON_R / fabs(term5) );

        TimeScale = min( delta_t1, delta_t2 );
/*
	// (1) Commenting in this block of code sets the ionization state to equilibrium in the lower chromosphere, but still allows it to be evolved by the flows
        if( TimeScale < MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE )
        {
#ifdef DENSITY_DEPENDENT_RATES
            pni2[iIndex] = GetEquilIonFrac( iIndex+1, flog_10T, flog_10n );
#else // DENSITY_DEPENDENT_RATES
            pni2[iIndex] = GetEquilIonFrac( iIndex+1, flog_10T );
#endif // DENSITY_DEPENDENT_RATES
            term5 = 0.0;
            TimeScale = LARGEST_DOUBLE;
        }
*/
	// (2) Commenting in this block of code evolves the ionization state on the minimum collision coupling timescale in the lower chromosphere
        if( TimeScale < MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE )
        {
            term5 *= ( TimeScale / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE );
	    TimeScale = MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
        }
    }
    else
        TimeScale = LARGEST_DOUBLE;

    if( TimeScale < SmallestTimeScale )
        SmallestTimeScale = TimeScale;

    pdnibydt[iIndex] = term1 + term5;
}

*pTimeScale = SmallestTimeScale;
}

double CElement::GetEmissivity( int iIon, double flog_10T, double flog_10n, double ni )
{
return GetIonEmissivity( iIon, flog_10T, flog_10n ) * ni;
}

double CElement::GetEmissivity( double flog_10T, double flog_10n, double *pni )
{
double Emiss = 0.0;
int i, iIndex;

for( i=0; i<NumIons; i++ )
{
    iIndex = pSpecNum[i] - 1;
    Emiss += GetIonEmissivity( pSpecNum[i], flog_10T, flog_10n ) * pni[iIndex];
}

return Emiss;
}
