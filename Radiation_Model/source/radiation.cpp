// ****
// *
// * Radiation Class Function Bodies for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/09/2020
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "radiation.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"


CRadiation::CRadiation( char *szFilename )
{
#if !defined (USE_POWER_LAW_RADIATIVE_LOSSES)
	// The power-law radiative losses ARE NOT being used
	// The atomic data are required to calculate the radiative losses
	Initialise( szFilename );
#else // !(USE_POWER_LAW_RADIATIVE_LOSSES)
	// The power-law radiative losses ARE being used
	#ifdef OPTICALLY_THICK_RADIATION
		#if defined (NLTE_CHROMOSPHERE) || defined(DECOUPLE_IONISATION_STATE_SOLVER)
			// If the NLTE chromosphere is being used or a non-equilibrium ionization calculation is being performed then the atomic data are required
			Initialise( szFilename );
		#endif // NLTE_CHROMOSPHERE || DECOUPLE_IONISATION_STATE_SOLVER
	#else // OPTICALLY_THICK_RADIATION
		#ifdef DECOUPLE_IONISATION_STATE_SOLVER
			// If a non-equilibrium ionisation calculation is being performed then the atomic data are required
			Initialise( szFilename );
		#endif // DECOUPLE_IONISATION_STATE_SOLVER
	#endif // OPTICALLY_THICK_RADIATION
#endif // !(USE_POWER_LAW_RADIATIVE_LOSSES)
}

CRadiation::~CRadiation( void )
{
#if !defined (USE_POWER_LAW_RADIATIVE_LOSSES)
	// The power-law radiative losses ARE NOT being used
	// The atomic data are required to calculate the radiative losses
	FreeAll();
#else // !(USE_POWER_LAW_RADIATIVE_LOSSES)
	// The power-law radiative losses ARE being used
	#ifdef OPTICALLY_THICK_RADIATION
		#if defined (NLTE_CHROMOSPHERE) || defined(DECOUPLE_IONISATION_STATE_SOLVER)
			// If the NLTE chromosphere is being used or a non-equilibrium ionization calculation is being performed then the atomic data are required
			FreeAll();
		#endif // NLTE_CHROMOSPHERE || DECOUPLE_IONISATION_STATE_SOLVER
	#else // OPTICALLY_THICK_RADIATION
		#ifdef DECOUPLE_IONISATION_STATE_SOLVER
			// If a non-equilibrium ionisation calculation is being performed then the atomic data are required
			FreeAll();
		#endif // DECOUPLE_IONISATION_STATE_SOLVER
	#endif // OPTICALLY_THICK_RADIATION
#endif // !(USE_POWER_LAW_RADIATIVE_LOSSES)
}

void CRadiation::Initialise( char *szFilename )
{
FILE *pFile;
char buffer1[16], buffer2[16], buffer3[16], buffer4[16];
char szRangesFilename[256], szAbundFilename[256], szEmissFilename[256], szRatesFilename[256], szIonFracFilename[256];
int i;

pFile = fopen( szFilename, "r" );

// Get the range definition filename
fscanf( pFile, "%s", buffer1 );
sprintf( szRangesFilename, "Radiation_Model/atomic_data/ranges/%s.rng", buffer1 );

// Get the ion emissivities
fscanf( pFile, "%s", buffer2 );

// Get the abundance set
fscanf( pFile, "%s", buffer3 );
sprintf( szAbundFilename, "Radiation_Model/atomic_data/abundances/%s.ab", buffer3 );

// Get the ionisation and recombination data
fscanf( pFile, "%s", buffer4 );

// Get the number of elements from the file
fscanf( pFile, "%i", &NumElements );

// Allocate sufficient memory to hold the pointers to each element object
ppElements = (PPELEMENT)malloc( sizeof( CElement* ) * NumElements );

// Allocate sufficient memory to hold the list of atomic numbers
pZ = (int*)malloc( sizeof(int) * NumElements );

for( i=0; i<NumElements; i++ )
{
    // Get the element symbol
    fscanf( pFile, "%s", buffer1 );
	
    // Get the atomic number
    fscanf( pFile, "%i", &(pZ[i]) );
	
    // Construct the filenames
    sprintf( szEmissFilename, "Radiation_Model/atomic_data/emissivities/%s/%s.em", buffer2, buffer1 );
    sprintf( szRatesFilename, "Radiation_Model/atomic_data/rates/%s/%s.rts", buffer4, buffer1 );
    sprintf( szIonFracFilename, "Radiation_Model/atomic_data/balances/%s/%s.bal", buffer4, buffer1 );

    // Instantiate each element object
    ppElements[i] = new CElement( pZ[i], szRangesFilename, szAbundFilename, szEmissFilename, szRatesFilename, szIonFracFilename );
}

fclose( pFile );

// Open the temperature and density ranges file and allocate memory to store these quantities
OpenRangesFile( szRangesFilename );

// Calculate the total phi of all radiating elements as a function of temperature and density
CalculateTotalPhi();

#ifdef USE_RADIATION_LOOKUP_TABLE
	OpenRadiationLookupTable();
#endif // USE_RADIATION_LOOKUP_TABLE
}

void CRadiation::OpenRangesFile( char *szRangesFilename )
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

void CRadiation::CalculateTotalPhi( void )
{
int i, j, k, l, NumTempxNumDen;

// Calculate the 2D array sizes
NumTempxNumDen = NumTemp * NumDen;

// Allocate an array to hold the total values of phi( n, T ) for the element
pTotalPhi = (double*)malloc( sizeof(double) * NumTempxNumDen );
for( l=0; l<NumTempxNumDen; l++ )
	pTotalPhi[l] = 0.0;

l=0;
for( k=0; k<NumDen; k++ )
    for( j=0; j<NumTemp; j++ )
    {
	for( i=0; i<NumElements; i++ )
	    pTotalPhi[l] += ppElements[i]->GetEmissivity( pTemp[j], pDen[k] );
	l++;
    }
}

#ifdef USE_RADIATION_LOOKUP_TABLE
void CRadiation::OpenRadiationLookupTable( void )
{
	FILE *pFile;
	int i;

	// Open the data file
	pFile = fopen( RADIATION_LOOKUP_TABLE, "r" );
		// Get the length (number of data points) in the lookup table
		fscanf( pFile, "%i", &iLookupTableLength);
		// Get the coefficients for the linear mapping from (log_10) temperature to array index for fast lookup
		ReadDouble( pFile, (&fMaptoIndex[0]) );
		ReadDouble( pFile, (&fMaptoIndex[1]) );
		// Allocate memory and read the data from the file
		ppfRadiationLookupTable = (double**)malloc( sizeof(double*) * iLookupTableLength );
		for( i=0; i<iLookupTableLength; i++ ) {
			ppfRadiationLookupTable[i] = (double*)malloc( sizeof(double) * 2 );
			ReadDouble( pFile, (&ppfRadiationLookupTable[i][0]) );
			ReadDouble( pFile, (&ppfRadiationLookupTable[i][1]) );
		}
	// Close the data file
	fclose( pFile );
}
#endif // USE_RADIATION_LOOKUP_TABLE

void CRadiation::FreeAll( void )
{
int i;

#ifdef USE_RADIATION_LOOKUP_TABLE
	for( i=0; i<iLookupTableLength; i++ ) {
		free( ppfRadiationLookupTable[i] );
	}
	free( ppfRadiationLookupTable );
#endif // USE_RADIATION_LOOKUP_TABLE

free( pTotalPhi );

free( pDen );
free( pTemp );

for( i=0; i<NumElements; i++ )
    delete ppElements[i];
	
free( ppElements );
free( pZ );
}

int* CRadiation::pGetAtomicNumbers( int *iNumElements )
{
*iNumElements = NumElements;
return pZ;
}

double CRadiation::GetAbundance( int iZ )
{
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return 0.0;

return ppElements[i]->GetAbundance();
}

void CRadiation::GetRates( int iZ, int iIon, double flog_10T, double *pfIonRate, double *pfRecRate )
{
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

ppElements[i]->GetRates( iIon, flog_10T, pfIonRate, pfRecRate );
}

void CRadiation::GetRates( int iZ, int iIon, double flog_10T, double flog_10n, double *pfIonRate, double *pfRecRate )
{
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

ppElements[i]->GetRates( iIon, flog_10T, flog_10n, pfIonRate, pfRecRate );
}

void CRadiation::GetEquilIonFrac( int iZ, double *pni, double flog_10T )
{
double fTotal = 0.0;
int i, j;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

// Get the set of equilibrium ion fractional populations for the specified element

for( j=0; j<=iZ; j++ )
{
    pni[j] = ppElements[i]->GetEquilIonFrac( j+1, flog_10T );
    fTotal += pni[j];
}

// Normalise the sum total of the ion fractional populations to 1
Normalise( iZ, pni, fTotal );
}

void CRadiation::GetEquilIonFrac( int iZ, double *pni, double flog_10T, double flog_10n )
{
double fTotal = 0.0;
int i, j;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

// Get the set of equilibrium ion fractional populations for the specified element

for( j=0; j<=iZ; j++ )
{
    pni[j] = ppElements[i]->GetEquilIonFrac( j+1, flog_10T, flog_10n );
    fTotal += pni[j];
}

// Normalise the sum total of the ion fractional populations to 1
Normalise( iZ, pni, fTotal );
}

void CRadiation::WriteEquilIonFracToFile( void *pFile, int iZ, double flog_10T )
{
int i, j;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

for( j=0; j<=iZ; j++ )
    fprintf( (FILE*)pFile, "\t%.8e", ppElements[i]->GetEquilIonFrac( j+1, flog_10T ) );

fprintf( (FILE*)pFile, "\n" );
}

void CRadiation::WriteEquilIonFracToFile( void *pFile, int iZ, double flog_10T, double flog_10n )
{
int i, j;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

for( j=0; j<=iZ; j++ )
    fprintf( (FILE*)pFile, "\t%.8e", ppElements[i]->GetEquilIonFrac( j+1, flog_10T, flog_10n ) );

fprintf( (FILE*)pFile, "\n" );
}

void CRadiation::WriteEquilIonFracToFile( void *pFile, double flog_10T )
{
int i, j;

for( i=0; i<NumElements; i++ )
{
    fprintf( (FILE*)pFile, "\n%i", pZ[i] );

    for( j=0; j<=pZ[i]; j++ )
        fprintf( (FILE*)pFile, "\t%.8e", ppElements[i]->GetEquilIonFrac( j+1, flog_10T ) );
}

fprintf( (FILE*)pFile, "\n" );
}

void CRadiation::WriteEquilIonFracToFile( void *pFile, double flog_10T, double flog_10n )
{
int i, j;

for( i=0; i<NumElements; i++ )
{
    fprintf( (FILE*)pFile, "\n%i", pZ[i] );

    for( j=0; j<=pZ[i]; j++ )
        fprintf( (FILE*)pFile, "\t%.8e", ppElements[i]->GetEquilIonFrac( j+1, flog_10T, flog_10n ) );
}

fprintf( (FILE*)pFile, "\n" );
}

void CRadiation::Normalise( int iZ, double *pni, double fTotal )
{
int i;

for( i=0; i<=iZ; i++ )
    pni[i] = pni[i] / fTotal;
}

double CRadiation::GetRadiation( int iZ, int iIon, double flog_10T, double fne, double fnH )
{
double fEmiss;
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return 0.0;

fEmiss = ppElements[i]->GetEmissivity( iIon, flog_10T, log10(fne) );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is NOT added here
}

double CRadiation::GetRadiation( int iZ, double flog_10T, double fne, double fnH )
{
double fEmiss;
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return 0.0;

fEmiss = ppElements[i]->GetEmissivity( flog_10T, log10(fne) );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is NOT added here
}

double CRadiation::GetRadiation( double flog_10T, double fne, double fnH )
{
double x1[5], x2[5], **y, *pfTemp, result, error, flog_10n = log10(fne);
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

// Select the sixteen phi( ne, T ) values corresponding to the grid

// Allocate an array of pointers to pointers for the values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the values corresponding to the specified range of densities and temperatures
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the phi( ne, T ) values for the element
    pfTemp = pTotalPhi;

    // Skip to the set corresponding to the l'th density value and get the values of total phi( ne, T ) corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog_10T, flog_10n, &result, &error );

// Check the value of phi( ne, T ) is physically realistic
if( result < 0.0 ) result = 0.0;

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * result;
// NOTE: free-free radiation is NOT added here
}

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiation::Getdnibydt( int iZ, double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double *pdnibydt, double *pTimeScale )
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiation::Getdnibydt( int iZ, double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double delta_s, double *pdnibydt, double *pTimeScale )
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
{
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return;

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
ppElements[i]->Getdnibydt( flog_10T, flog_10n, pni0, pni1, pni2, pni3, pni4, s, s_pos, pv, cross_section, cell_volume, pdnibydt, pTimeScale );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
ppElements[i]->Getdnibydt( flog_10T, flog_10n, pni0, pni1, pni2, pni3, pni4, s, s_pos, pv, delta_s, pdnibydt, pTimeScale );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
}

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiation::GetAlldnibydt( double flog_10T, double flog_10n, double **ppni0, double **ppni1, double **ppni2, double **ppni3, double **ppni4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double **ppdnibydt, double *pTimeScale )
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiation::GetAlldnibydt( double flog_10T, double flog_10n, double **ppni0, double **ppni1, double **ppni2, double **ppni3, double **ppni4, double *s, double *s_pos, double *pv, double delta_s, double **ppdnibydt, double *pTimeScale )
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
{
double TimeScale, SmallestTimeScale;
int i;

SmallestTimeScale = LARGEST_DOUBLE;

for( i=0; i<NumElements; i++ )
{
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    ppElements[i]->Getdnibydt( flog_10T, flog_10n, ppni0[i], ppni1[i], ppni2[i], ppni3[i], ppni4[i], s, s_pos, pv, cross_section, cell_volume, ppdnibydt[i], &TimeScale );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
    ppElements[i]->Getdnibydt( flog_10T, flog_10n, ppni0[i], ppni1[i], ppni2[i], ppni3[i], ppni4[i], s, s_pos, pv, delta_s, ppdnibydt[i], &TimeScale );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
    if( TimeScale < SmallestTimeScale )
        SmallestTimeScale = TimeScale;
}

*pTimeScale = SmallestTimeScale;
}

double CRadiation::GetRadiation( int iZ, int iIon, double flog_10T, double fne, double fnH, double ni )
{
double fEmiss;
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return 0.0;

fEmiss = ppElements[i]->GetEmissivity( iIon, flog_10T, log10(fne), ni );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is NOT added here
}

double CRadiation::GetRadiation( int iZ, double flog_10T, double fne, double fnH, double *pni )
{
double fEmiss;
int i;

// Find the required element
for( i=0; i<NumElements; i++ )
    if( iZ == pZ[i] ) break;

if( i == NumElements ) return 0.0;

fEmiss = ppElements[i]->GetEmissivity( flog_10T, log10(fne), pni );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is NOT added here
}

double CRadiation::GetRadiation( double flog_10T, double fne, double fnH, double **ppni )
{
double fEmiss = 0.0;
int i;

for( i=0; i<NumElements; i++ )
    fEmiss += ppElements[i]->GetEmissivity( flog_10T, log10(fne), ppni[i] );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is NOT added here
}

double CRadiation::GetPowerLawRad( double flog_10T, double fne, double fnH )
{
#ifdef USE_RADIATION_LOOKUP_TABLE
	// If a lookup table is being used then hijack this function to avoid further complicating the calculation of the 
	// radiation term in CEquations::EvaluateTerms. Either power laws or a lookup table can be used to find the radiative
	// losses, but not both at the same time, and so this should avoid introducing any bugs / conditions not accounted for 
	return GetRadiationFromLookupTable( flog_10T, fne, fnH );
#else // USE_RADIATION_LOOKUP_TABLE
double chi, alpha, fEmiss;

// The formulation used here is based on the calculations of John Raymond (1994, private communication) and twice the coronal abundances of Meyer (1985)

if( flog_10T <= 4.97 )
{
    chi = 1.09e-31;
    alpha = 2.0;
}
else if( flog_10T <= 5.67 )
{
    chi = 8.87e-17;
    alpha = -1.0;
}
else if( flog_10T <= 6.18 )
{
    chi = 1.90e-22;
    alpha = 0.0;
}
else if( flog_10T <= 6.55 )
{
    chi = 3.53e-13;
    // alpha = -3.0/2.0;
    alpha = -1.50;
}
else if( flog_10T <= 6.90 )
{
    chi = 3.46e-25;
    // alpha = 1.0/3.0;
    alpha = 0.33333333333333333333333333333333;
}
else if( flog_10T <= 7.63 )
{
    chi = 5.49e-16;
    alpha = -1.0;
}
else
{
    chi = 1.96e-27;
    // alpha = 1.0/2.0;
    alpha = 0.50;
}

fEmiss = chi * pow( 10.0, (alpha*flog_10T) );

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * fEmiss;
// NOTE: free-free radiation is included in the parameter values for log_10 T > 7.63
#endif // USE_RADIATION_LOOKUP_TABLE
}

double CRadiation::GetFreeFreeRad( double flog_10T, double fne, double fnH )
{
// Calculate the free-free emission due to thermal bremsstralung for a low-density and fully ionised plasma
// The previous formulation used here was taken from Mason & Monsignori Fossi (1994, Astron. Astrophys. Rev., 6, 123) where (1.96e-27) is replaced with (2.40e-27)
// The current formulation is taken from the power-law fit for log_10 T > 7.63

if( fne > MAX_OPTICALLY_THIN_DENSITY )
{
    fne = MAX_OPTICALLY_THIN_DENSITY;
    if( fnH > fne )
	fnH = fne;
}

return ( fne * fnH ) * (1.96e-27) * pow( 10.0, (0.5*flog_10T) );
}

#ifdef USE_RADIATION_LOOKUP_TABLE
double CRadiation::GetRadiationFromLookupTable( double flog_10T, double fne, double fnH )
{
	double x[5], y[5], result;
	int i;
	
	// Calculate the array index from (log_10) temperature using the linear mapping
	i = (int)(( fMaptoIndex[0] * flog_10T ) + fMaptoIndex[1]);
	
	// Check for array indices that are out-of-bounds
	if( i < 0 ) i = 0;
	else if ( i > iLookupTableLength - 2 ) i = iLookupTableLength - 2;
	
	x[1] = ppfRadiationLookupTable[i][0];
	x[2] = ppfRadiationLookupTable[i+1][0];

	y[1] = ppfRadiationLookupTable[i][1];
	y[2] = ppfRadiationLookupTable[i+1][1];

	LinearFit( x, y, flog_10T, &result );

	if( fne > MAX_OPTICALLY_THIN_DENSITY )
	{
    	fne = MAX_OPTICALLY_THIN_DENSITY;
    	if( fnH > fne )
		fnH = fne;
	}
	
	return ( fne * fnH ) * pow( 10.0, result );
}
#endif // USE_RADIATION_LOOKUP_TABLE