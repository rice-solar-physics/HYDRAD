// ****
// *
// * Ion Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 25/09/2013
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "ion.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"


CIon::CIon( char *pszLabel, int iZ, int iSpecNum, char *szEmissFilename, double *pRespFunc, int iNumDataPoints )
{
Initialise( pszLabel, iZ, iSpecNum, szEmissFilename, pRespFunc, iNumDataPoints );
}

CIon::~CIon( void )
{
FreeAll();
}

void CIon::Initialise( char *pszLabel, int iZ, int iSpecNum, char *szEmissFilename, double *pRespFunc, int iNumDataPoints )
{
FILE *pFile;
char fBuffer[256];
int iNumLines, i, iProtons, iNeutrons;

// Set the ion properties
sprintf( szLabel, "%s", pszLabel );
Z = iZ;
SpecNum = iSpecNum;

// Calculate the ion mass
pFile = fopen( "Radiation_Model/atomic_data/masses/masses.ma", "r" );
fscanf( pFile, "%i", &iNumLines );
fscanf( pFile, "%s", fBuffer ); fscanf( pFile, "%s", fBuffer );
for( i=0; i<iNumLines; i++ )
{
    // Get the number of protons and neutrons
    fscanf( pFile, "%i", &iProtons ); fscanf( pFile, "%i", &iNeutrons );
    if( iProtons == iZ ) break;
}
mi = ( ((double)iProtons) * PROTON_MASS ) + ( ((double)iNeutrons) * NEUTRON_MASS );
fclose( pFile );

// Open the ion emissivity data files and initialise the ion object
OpenEmissivityFile( szEmissFilename, pRespFunc, iNumDataPoints );

// Calculate the total line emission as a function of temperature and density (summed over the wavelength range)
CalculateTotalIonEmission();
}

void CIon::OpenEmissivityFile( char *szEmissFilename, double *pRespFunc, int iNumDataPoints )
{
FILE *pFile;
char szBuffer[256];
double fLambda_min, fLambda_max, *pfLambda, x[3], y[3], fIRC, fBuffer;
double term1, term2;
int iNumLines, iIndexDen, iIndexTemp, i, j;

// Get the wavelength sensitivity range of the instrument
fLambda_min = pRespFunc[0];
fLambda_max = pRespFunc[iNumDataPoints-2];

pFile = fopen( szEmissFilename, "r" );

// Get the ion label
fscanf( pFile, "%s", szBuffer );

// Get the total number of lines
fscanf( pFile, "%i", &iNumLines );

// Allocate sufficient temporary memory for the wavelengths
pfLambda = (double*)malloc( sizeof(double) * iNumLines );

// Get the wavelength range of the lines stored in the file
ReadDouble( pFile, &fBuffer );
ReadDouble( pFile, &fBuffer );

// Get the wavelengths of the lines

NumLines = 0;

for( i=0; i<iNumLines; i++ )
{
    ReadDouble( pFile, &(pfLambda[i]) );

    if( pfLambda[i] >= fLambda_min && pfLambda[i] <= fLambda_max )
        NumLines++;
}

// Allocate sufficient memory for the wavelengths in the sensitivity range 
// of the instrument 
pLambda = (double*)malloc( sizeof(double) * NumLines );

j = 0;

for( i=0; i<iNumLines; i++ )
{
    if( pfLambda[i] >= fLambda_min && pfLambda[i] <= fLambda_max )
    {
        pLambda[j] = pfLambda[i];
        j++;
    }
}

free( pfLambda );

// Get the number of density values
fscanf( pFile, "%i", &NumDen );

// Allocate sufficient memory for the density values
pDen = (double*)malloc( sizeof(double) * NumDen );

// Get the density values
for( i=0; i<NumDen; i++ )
    ReadDouble( pFile, &(pDen[i]) );

// Get the number of temperature values
fscanf( pFile, "%i", &NumTemp );

// Allocate sufficient memory for the temperature values
pTemp = (double*)malloc( sizeof(double) * NumTemp );

// Get the temperature values
for( i=0; i<NumTemp; i++ )
    ReadDouble( pFile, &(pTemp[i]) );

// Get the emissivity values

NumTempxNumDen = NumDen * NumTemp;

// Allocate sufficient memory for the emissivity values

ppEmiss = (double**)malloc( sizeof(double) * NumLines );

i = 0;

for(;;)
{
    iIndexDen = 0;
    iIndexTemp = 0;

    // Get the wavelength of the line
    ReadDouble( pFile, &fBuffer );

    // Necessary if the data for the final emission line in the file is read
    if( feof(pFile) ) break;

    // If the wavelength of the line is greater than the maximum wavelength
    // of the instrument sensitivity range then there is no point reading further
    if( fBuffer > fLambda_max ) break;

    if( fBuffer >= fLambda_min && fBuffer <= fLambda_max )
    {
        fIRC = 0.0;

        // Get the instrument response coefficient
        for( j=2; j<iNumDataPoints; j+=2 )
        {
            if( pRespFunc[j] >= pLambda[i] )
	    {
	        x[1] = pRespFunc[j-2];
		x[2] = pRespFunc[j];
		y[1] = pRespFunc[j-1];
		y[2] = pRespFunc[j+1];
		      
		// Perform a linear fit to the instrument response function because
		// the curve can vary both smoothly and extremely rapidly in different
		// wavelength regions
		LinearFit( x, y, pLambda[i], &fIRC );

		break;
            }
        }

        // Allocate sufficient memory for the emissivity values corresponding
	// to each wavelength
	ppEmiss[i] = (double*)malloc( sizeof(double) * NumTempxNumDen );

        // Get the individual emissivity values
        for( j=0; j<NumTempxNumDen; j++ )
        {
            // The values in the .em files have units of [erg s^-1] and were
            // calculated using the Chianti function EMISS_CALC
            ReadDouble( pFile, &fBuffer );
#ifdef UNITS_DN
            // Calculate the value to store, assuming:
            // the instrument response function in units of [DN pixel^-1 photon^-1 sr cm^2]
            // 0.83 = proton:electron ratio [dimensionless]
            // pLambda[i] = wavelength [A]
            // pow(10.0, pDen[iIndexDen]) = electron number density [cm^-3]
            // 2.495681203e-7 = hc * 4 _PI_ [erg A sr]
            term1 = fIRC * 0.83 * pLambda[i];
            term2 = pow(10.0, pDen[iIndexDen]) * (2.495681203e-7);
            // The stored value has units of [DN pixel^-1 photon^-1 sr cm^2] * [photon cm^3 s^-1 sr^-1]
            // --> [DN pixel^-1 s^-1 cm^5]
#endif // UNITS_DN
#ifdef UNITS_PHOTONS
            // Calculate the value to store, assuming:
            // 0.83 = proton:electron ratio [dimensionless]
            // pLambda[i] = wavelength [A]
            // pow(10.0, pDen[iIndexDen]) = electron number density [cm^-3]
            // 2.495681203e-7 = hc * 4 _PI_ [erg A sr]
            term1 = 0.83 * pLambda[i];
            term2 = pow(10.0, pDen[iIndexDen]) * (2.495681203e-7);
            // The stored value has units of [photon cm^3 s^-1 sr^-1]
#endif // UNITS_PHOTONS
#ifdef UNITS_ERG
            // Calculate the value to store, assuming:
            // 0.83 = proton:electron ratio [dimensionless]
            // pow(10.0, pDen[iIndexDen]) = electron number density [cm^-3]
            // 12.566370614359172953850573533118 = 4 _PI_ [sr]
            term1 = 0.83;
            term2 = pow(10.0, pDen[iIndexDen]) * (12.566370614359172953850573533118);
            // The stored value has units of [erg cm^3 s^-1 sr^-1]
#endif // UNITS_ERG
            ppEmiss[i][j] = fBuffer * ( term1 / term2 );
            // The stored value must then be multiplied by:
            // the element abundance relative to hydrogen [dimensionless]
            // the population fraction [dimensionless]
            // the volume emission measure = n^2 dV [cm^-3]
            // the inverse of the detector pixel area = 1 / A_pixel [cm^-2]
            // (OR, instead of the previous two values the column emission measure = n^2 dR [cm^-5] - in special case of perfect resolution)

            iIndexTemp++;

            if( iIndexTemp == NumTemp )
            {
                iIndexTemp = 0;
                iIndexDen++;
            }
        }

        i++;
    }
    else
    {
        // Skip the individual emissivity values
        for( j=0; j<NumTempxNumDen; j++ )
            ReadDouble( pFile, &fBuffer );
    }
}

fclose( pFile );
}

void CIon::CalculateTotalIonEmission( void )
{
double fTotal;
int iOffset;
int i, j, k;

pIonEmiss = (double*)malloc( sizeof(double) * NumTempxNumDen );

for( k=0; k<NumDen; k++ )
{
    iOffset = k * NumTemp;
    for( j=0; j<NumTemp; j++ )
    {
        fTotal = 0.0;
        // Sum the line emission over the wavelength range
        for( i=0; i<NumLines; i++ )
            fTotal += GetLineEmission( pLambda[i], pTemp[j], pDen[k] );

        pIonEmiss[iOffset++] = fTotal;
    }
}
}

void CIon::FreeAll( void )
{
int i;

// Free the allocated memory

free( pIonEmiss );

for( i=0; i<NumLines; i++ )
    free( ppEmiss[i] );

free( ppEmiss );

free( pTemp );
free( pDen );
free( pLambda );
}

int CIon::GetAtomicNumber( void )
{
return Z;
}

int CIon::GetSpecNumber( void )
{
return SpecNum;
}

double CIon::GetMass( void )
{
return mi;
}

int CIon::GetNumLines( void )
{
return NumLines;
}

void CIon::GetLineList( double *pfLineList )
{
int i;

for( i=0; i<NumLines; i++ )
    pfLineList[i] = pLambda[i];
}

void CIon::WriteIonLineListToFile( FILE *pFile )
{
int i;

fprintf( pFile, "\n\t%s: ", szLabel );

if( NumLines > 0 )
{
    fprintf( pFile, "%i lines (", NumLines );
    for( i=0; i<NumLines; i++ )
	fprintf( pFile, " %.3f", pLambda[i] );
    fprintf( pFile, " )\n" );
}
else
{
    fprintf( pFile, "%i lines\n", NumLines );
}
}

double CIon::GetLineEmission( double fLambda, double flog10T, double flog10n )
{
double x1[5], x2[5], **y, *pfTemp, result, error;
int i, j, k, l;

// There is no emission below the cut-off temperature
if( flog10T <= LOG_LINE_EMISSION_CUTOFF_TEMP ) return 0.0;

// Select the required line

for( i=0; i<NumLines; i++ )
{
    if( fLambda == pLambda[i] ) break;

#ifdef CHOOSE_NEAREST_WAVELENGTH
    // If the specified wavelength is not exactly equal to the stored
    // wavelength then select the first one that is greater than the
    // specified wavelength
    if( fLambda < pLambda[i] ) break;
#endif // CHOOSE_NEAREST_WAVELENGTH
}

// If the requested line is not available then return 0.0
if( i == NumLines ) return 0.0;

// Select the four temperature values surrounding the desired one

// If the temperature is out of range then return 0.0
if( flog10T < pTemp[0] || flog10T > pTemp[NumTemp-1] ) return 0.0;

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog10T ) break;

// Deal with the special cases where there aren't two values either side of the
// desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then return 0.0
if( flog10n < pDen[0] || flog10n > pDen[NumDen-1] ) return 0.0;

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog10n ) break;

// Deal with the special cases where there aren't two values either side of the
// desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the i'th line
// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen line emission values corresponding to the grid

// Allocate an array of pointers to pointers for the line emission values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the line emission values corresponding to the set of
// temperatures at the current density
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the line emission values for the required line
    pfTemp = ppEmiss[i];

    // Skip to the line emission set corresponding to the l'th density value and get the
    // line emissions corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog10T, flog10n, &result, &error );

// Check the result is physically realistic
if( result < 0.0 ) result = 0.0;

// Return the line emission
return result;
}

double CIon::GetIonEmission( double flog10T, double flog10n )
{
double x1[5], x2[5], **y, *pfTemp, result, error;
int j, k, l;

// There is no emission below the cut-off temperature
if( flog10T <= LOG_LINE_EMISSION_CUTOFF_TEMP ) return 0.0;

// Select the four temperature values surrounding the desired one

// If the temperature is out of range then return 0.0
if( flog10T < pTemp[0] || flog10T > pTemp[NumTemp-1] ) return 0.0;

for( j=0; j<NumTemp; j++ )
    if( pTemp[j] >= flog10T ) break;

// Deal with the special cases where there aren't two values either side of the
// desired one
if( j < 2 ) j = 2;
else if( j == NumTemp-1 ) j = NumTemp-2;

x1[1] = pTemp[j-2];
x1[2] = pTemp[j-1];
x1[3] = pTemp[j];
x1[4] = pTemp[j+1];

// Select the four density values surrounding the desired one

// If the density is out of range then return 0.0
if( flog10n < pDen[0] || flog10n > pDen[NumDen-1] ) return 0.0;

for( k=0; k<NumDen; k++ )
    if( pDen[k] >= flog10n ) break;

// Deal with the special cases where there aren't two values either side of the
// desired one
if( k < 2 ) k = 2;
else if( k == NumDen-1 ) k = NumDen-2;

x2[1] = pDen[k-2];
x2[2] = pDen[k-1];
x2[3] = pDen[k];
x2[4] = pDen[k+1];

// We are using the j-2, j-1, j and j+1 'th temperature values
// We are using the k-2, k-1, k and k+1 'th density values

// Select the sixteen ion emission values corresponding to the grid

// Allocate an array of pointers to pointers for the ion emission values
y = (double**)alloca( sizeof(double) * 5 );

// Allocate an array to hold the ion emission values corresponding to the set of
// temperatures at the current density
for( l=1; l<=4; l++ )
    y[l] = (double*)alloca( sizeof(double) * 5 );

for( l=1; l<=4; l++ )
{
    // Point to the start of the ion emission values
    pfTemp = pIonEmiss;

    // Skip to the ion emission set corresponding to the l'th density value and get the
    // ion emissions corresponding to the four temperature values
    pfTemp += ( k + l - 3 ) * NumTemp;

    y[1][l] = *( pfTemp + ( j - 2 ) );
    y[2][l] = *( pfTemp + ( j - 1 ) );
    y[3][l] = *( pfTemp + j );
    y[4][l] = *( pfTemp + ( j + 1 ) );
}

// Perform the 2D polynomial interpolation
FitPolynomial2D( x1, x2, y, 4, 4, flog10T, flog10n, &result, &error );

// Check the result is physically realistic
if( result < 0.0 ) result = 0.0;

// Return the ion emission
return result;
}