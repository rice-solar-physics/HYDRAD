// ****
// *
// * Instrument Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/19/2018
// *
// ****


#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "instrument.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"


CInstrument::CInstrument( char *pszName )
{
Initialise( pszName );
}

CInstrument::~CInstrument( void )
{
FreeAll();
}

void CInstrument::Initialise( char *pszName )
{
FILE *pFile;
char szFilename[256];

// Set the name of the instrument
sprintf( szName, "%s", pszName );

printf( "Initialising: %s\n", szName );

// Create the virtual instrument configuration file path
sprintf( szFilename, "Forward_Model/instruments/%s.ins", szName );

// Open the configuration file
pFile = fopen( szFilename, "r" );

// Read the instrument response function part of the file into memory
GetRespFunc( pFile );

// Read the ion data part of the file into memory
GetIonData( pFile );

// Close the configuration file
fclose( pFile );

// Initialise the virtual detectors
ppfSolarXY = NULL;
pfEQDetector = NULL;
pfNEQDetector = NULL;
pfLambda = NULL;
ppfEQDetector_SPEC = NULL;
ppfNEQDetector_SPEC = NULL;

// Initialise the Radiation objects
pEQRadiation = new CRadiation( "Radiation_Model/config/elements_eq.cfg" );
pNEQRadiation = new CRadiation( "Radiation_Model/config/elements_neq.cfg" );

printf( "\n\n" );
}

void CInstrument::GetRespFunc( FILE *pFile )
{
int i;

// Get the number of data points in the file
fscanf( pFile, "%i", &iNumRespFuncDataPoints );
iNumRespFuncDataPoints *= 2;

// Allocate sufficient memory for the response function
pRespFunc = (double*)malloc( sizeof(double) * iNumRespFuncDataPoints );

// Read the response function values as a function of wavelength from the data file
// The response function has units of: Effective area [cm^2 electron photon^-1] * platescale [sr pixel^-1] / gain [electron DN^-1] --> [DN pixel^-1 photon^-1 sr cm^2]
for( i=0; i<iNumRespFuncDataPoints; i++ )
    ReadDouble( pFile, &(pRespFunc[i]) );

// Get the wavelength sensitivity range of the instrument
fLambda_min = pRespFunc[0];
fLambda_max = pRespFunc[iNumRespFuncDataPoints-2];

// Get the spectral resolution of the instrument
ReadDouble( pFile, &fSpecRes );

// Get the instrumental line width
ReadDouble( pFile, &fInstrumentalLineWidth );

// Get the spatial and temporal resolution of the instrument
for( i=0; i<2; i++ )
    ReadDouble( pFile, &(fPixel_arcsec[i]) );

fPixel_Area_arcsec = 1.0;
fPixel_Area_cm = 1.0;
for( i=0; i<2; i++ )
{
    fPixel_cm[i] = fPixel_arcsec[i] * CM_PER_ARCSEC;
    fPixel_Area_arcsec *= fPixel_arcsec[i];
    fPixel_Area_cm *= fPixel_cm[i];
}

ReadDouble( pFile, &fCadence );
}

void CInstrument::GetIonData( FILE *pFile )
{
FILE *pEmissFilePath;
char szLabel[16], szEmissFilePath[256], szEmissFilename[256];
int iZ, iSpecNum, i;

// Get the path of the wavelength resolved ion emissivity files
pEmissFilePath = fopen( "Forward_Model/ion_emiss_tables/ion_emiss_tables.path", "r" );
fscanf( pEmissFilePath, "%s", szEmissFilePath );
fclose( pEmissFilePath );

// **** IMPORTANT ****
// The ions should be stored in order of increasing atomic number and spectroscopic number
// ie. HI, HII, HeI, HeII, HeIII, .... , FeI, FeII, FeIII, FeIV, etc.

// Get the number of emitting ions
fscanf( pFile, "%i", &iNumIons );

// Allocate sufficient memory for the ion objects
ppIon = (PPION)malloc( sizeof(CIon*) * iNumIons );

printf( "Adding ion: " );

for( i=0; i<iNumIons; i++ )
{
    // Get the ion label. Format: fe_X  (where X is the spectroscopic number)
    //                            fe_Xd (where d signifies 'relative to the satellite lines')
    fscanf( pFile, "%s", szLabel );
    printf( " %s", szLabel );

    // Get the atomic number
    fscanf( pFile, "%i", &iZ );

    // Get the spectroscopic number
    fscanf( pFile, "%i", &iSpecNum );

    // Construct the emissivity filename
    sprintf( szEmissFilename, "%s/%s.em", szEmissFilePath, szLabel );

    // Instantiate each ion object
    ppIon[i] = new CIon( szLabel, iZ, iSpecNum, szEmissFilename, pRespFunc, iNumRespFuncDataPoints );
}
}

void CInstrument::FreeAll( void )
{
int i;

delete pNEQRadiation;
delete pEQRadiation;

DeleteVirtualDetector();

for( i=0; i<iNumIons; i++ )
    delete ppIon[i];

free( ppIon );

free( pRespFunc );
}

double CInstrument::GetMass( int iZ )
{
int i;
double fmi = 0.0;

// Select the required ion

for( i=0; i<iNumIons; i++ )
{
    if( iZ == ppIon[i]->GetAtomicNumber() )
        fmi = ppIon[i]->GetMass();

    // No point in continuing if the search has gone past the specified ion
    if( iZ <= ppIon[i]->GetAtomicNumber() ) break;
}

return fmi;
}

int CInstrument::GetNumLines( int iZ, int iSpecNum )
{
int i, iNumLines = 0;

// Select the required ion

for( i=0; i<iNumIons; i++ )
{
    if( ( iZ == ppIon[i]->GetAtomicNumber() ) && ( iSpecNum == ppIon[i]->GetSpecNumber() ) )
       iNumLines = ppIon[i]->GetNumLines();

    // No point in continuing if the search has gone past the specified ion
    if( ( iZ <= ppIon[i]->GetAtomicNumber() ) && ( iSpecNum < ppIon[i]->GetSpecNumber() ) ) break;
}

return iNumLines;
}

void CInstrument::GetLineList( int iZ, int iSpecNum, double *pfLineList )
{
int i;

// Select the required ion

for( i=0; i<iNumIons; i++ )
{
    if( ( iZ == ppIon[i]->GetAtomicNumber() ) && ( iSpecNum == ppIon[i]->GetSpecNumber() ) )
       ppIon[i]->GetLineList( pfLineList );

    // No point in continuing if the search has gone past the specified ion
    if( ( iZ <= ppIon[i]->GetAtomicNumber() ) && ( iSpecNum < ppIon[i]->GetSpecNumber() ) ) break;
}
}

double CInstrument::GetLineEmission( int iZ, int iSpecNum, double fLambda, double flog10T, double flog10n )
{
double fLineEmission;
int i;

fLineEmission = 0.0;

// Select the required ion

for( i=0; i<iNumIons; i++ )
{
    if( ( iZ == ppIon[i]->GetAtomicNumber() ) && ( iSpecNum == ppIon[i]->GetSpecNumber() ) )
        fLineEmission += ppIon[i]->GetLineEmission( fLambda, flog10T, flog10n );

    // No point in continuing if the search has gone past the specified ion
    if( ( iZ <= ppIon[i]->GetAtomicNumber() ) && ( iSpecNum < ppIon[i]->GetSpecNumber() ) ) break;
}

return fLineEmission;
}

double CInstrument::GetIonEmission( int iZ, int iSpecNum, double flog10T, double flog10n )
{
double fIonEmission;
int i;

fIonEmission = 0.0;

// Select the required ion

for( i=0; i<iNumIons; i++ )
{
    if( ( iZ == ppIon[i]->GetAtomicNumber() ) && ( iSpecNum == ppIon[i]->GetSpecNumber() ) )
        fIonEmission += ppIon[i]->GetIonEmission( flog10T, flog10n );

    // No point in continuing if the search has gone past the specified ion
    if( ( iZ <= ppIon[i]->GetAtomicNumber() ) && ( iSpecNum < ppIon[i]->GetSpecNumber() ) ) break;
}

return fIonEmission;
}

void CInstrument::Detect( double fs, double fds, double fSD, double fLength, double fEQCounts, double fNEQCounts )
{
double fX_pos_arcsec[3], fPixel_pos_arcsec[3];
double fTheta, fThetaL, fThetaR, fLOS;
double fsL, fsR, fdsL, fdsR, fdsM;
double fTemp;
int iStartPixel, iEndPixel, iRemainingPixels, i;

// Don't try to use the detectors if they have not been created
if( !pfEQDetector || !pfNEQDetector )
    return;

// Convert the field-aligned grid cell left boundary, centre and right boundary coordinates to solar-X displacement in arcsec
fX_pos_arcsec[1] = ( fR * ( 1.0 - cos(fs/fR) ) ) / CM_PER_ARCSEC;

fsL = fs - ( fds / 2.0 );
fX_pos_arcsec[0] = ( fR * ( 1.0 - cos(fsL/fR) ) ) / CM_PER_ARCSEC;

fsR = fsL + fds;
fX_pos_arcsec[2] = ( fR * ( 1.0 - cos(fsR/fR) ) ) / CM_PER_ARCSEC;

iStartPixel = 0;
for( i=0; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[0] < fPixel_pos_arcsec[2] )
        break;

    iStartPixel++;
}

iEndPixel = iStartPixel;
for( i=iStartPixel; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[2] <= fPixel_pos_arcsec[2] )
        break;

    iEndPixel++;
}

if( iStartPixel == iEndPixel )
{
    // The grid cell is entirely contained within the pixel
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fTheta = ( _PI_ * fs ) / fLength;
    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fds * cos(fTheta) * cos(fTheta) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );
}
else
{
    // Calculate the fraction of emission for the left-most pixel
    // Find the field-aligned coordinate 's' of the right-hand boundary of the start pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iStartPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[2]*CM_PER_ARCSEC) / fR ) );
    fdsL = fTemp - fsL;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaL = ( _PI_ * (0.5*(fsL+fTemp)) ) / fLength;
    fLOS = ( fdsL * cos(fThetaL) * cos(fThetaL) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );

    // Calculate the fraction of emission for the right-most pixel
    // Find the field-aligned coordinate 's' of the left-hand boundary of the end pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iEndPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[0]*CM_PER_ARCSEC) / fR ) );
    fdsR = fsR - fTemp;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaR = ( _PI_ * (0.5*(fTemp+fsR)) ) / fLength;
    fLOS = ( fdsR * cos(fThetaR) * cos(fThetaR) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iEndPixel] += ( fLOS * fNEQCounts );

    fTemp = fSD * ( fdsL / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaL) * sin(fThetaL);
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );

    fTemp = fSD * ( fdsR / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaR) * sin(fThetaR);
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iEndPixel] += ( fLOS * fNEQCounts );

    // Calculate the fraction of emission for the remaining pixels
    iRemainingPixels = ( iEndPixel - iStartPixel - 1 );
    if( iRemainingPixels > 0 )
    {
        fdsM = ( fds - fdsL - fdsR ) / ((double)iRemainingPixels);
	fTemp = fsL + fdsL + ( fdsM / 2.0 );
	for( i=iStartPixel+1; i<iEndPixel; i++ )
	{
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
	    fTheta = ( _PI_ * fTemp ) / fLength;
	    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fdsM * cos(fTheta) * cos(fTheta) );
	    fTemp += fdsM;
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
            pfEQDetector[i] += ( fLOS * fEQCounts );
            pfNEQDetector[i] += ( fLOS * fNEQCounts );
	}
    }
}
}

void CInstrument::Detect( double fs, double fds, double fSD, double fLength, double fEQCounts )
{
double fX_pos_arcsec[3], fPixel_pos_arcsec[3];
double fTheta, fThetaL, fThetaR, fLOS;
double fsL, fsR, fdsL, fdsR, fdsM;
double fTemp;
int iStartPixel, iEndPixel, iRemainingPixels, i;

// Don't try to use the detector if it has not been created
if( !pfEQDetector )
    return;

// Convert the field-aligned grid cell left boundary, centre and right boundary coordinates to solar-X displacement in arcsec
fX_pos_arcsec[1] = ( fR * ( 1.0 - cos(fs/fR) ) ) / CM_PER_ARCSEC;

fsL = fs - ( fds / 2.0 );
fX_pos_arcsec[0] = ( fR * ( 1.0 - cos(fsL/fR) ) ) / CM_PER_ARCSEC;

fsR = fsL + fds;
fX_pos_arcsec[2] = ( fR * ( 1.0 - cos(fsR/fR) ) ) / CM_PER_ARCSEC;

iStartPixel = 0;
for( i=0; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[0] < fPixel_pos_arcsec[2] )
        break;

    iStartPixel++;
}

iEndPixel = iStartPixel;
for( i=iStartPixel; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[2] <= fPixel_pos_arcsec[2] )
        break;

    iEndPixel++;
}

if( iStartPixel == iEndPixel )
{
    // The grid cell is entirely contained within the pixel
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fTheta = ( _PI_ * fs ) / fLength;
    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fds * cos(fTheta) * cos(fTheta) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
}
else
{
    // Calculate the fraction of emission for the left-most pixel
    // Find the field-aligned coordinate 's' of the right-hand boundary of the start pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iStartPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[2]*CM_PER_ARCSEC) / fR ) );
    fdsL = fTemp - fsL;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaL = ( _PI_ * (0.5*(fsL+fTemp)) ) / fLength;
    fLOS = ( fdsL * cos(fThetaL) * cos(fThetaL) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );

    // Calculate the fraction of emission for the right-most pixel
    // Find the field-aligned coordinate 's' of the left-hand boundary of the end pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iEndPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[0]*CM_PER_ARCSEC) / fR ) );
    fdsR = fsR - fTemp;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaR = ( _PI_ * (0.5*(fTemp+fsR)) ) / fLength;
    fLOS = ( fdsR * cos(fThetaR) * cos(fThetaR) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );

    fTemp = fSD * ( fdsL / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaL) * sin(fThetaL);
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );

    fTemp = fSD * ( fdsR / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaR) * sin(fThetaR);
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );

    // Calculate the fraction of emission for the remaining pixels
    iRemainingPixels = ( iEndPixel - iStartPixel - 1 );
    if( iRemainingPixels > 0 )
    {
        fdsM = ( fds - fdsL - fdsR ) / ((double)iRemainingPixels);
	fTemp = fsL + fdsL + ( fdsM / 2.0 );
	for( i=iStartPixel+1; i<iEndPixel; i++ )
	{
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
	    fTheta = ( _PI_ * fTemp ) / fLength;
	    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fdsM * cos(fTheta) * cos(fTheta) );
	    fTemp += fdsM;
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
	    pfEQDetector[i] += ( fLOS * fEQCounts );
	}
    }
}
}

void CInstrument::Detect( PHYDATA PHYData, double fSD, double fLength, double fEQCounts, double fNEQCounts, double fRestLambda, double *fSpectralProperties )
{
double fX_pos_arcsec[3], fPixel_pos_arcsec[3];
double fTheta, fThetaL, fThetaR, fLOS;
double fs, fds, fsL, fsR, fdsL, fdsR, fdsM;
double fIndex, fTemp;
int iStartPixel, iEndPixel, iRemainingPixels, i, j;

// Don't try to use the detectors if they have not been created
if( !ppfEQDetector_SPEC || !ppfNEQDetector_SPEC )
    return;

fs = PHYData.fs;
fds = PHYData.fds;

// Convert the field-aligned grid cell left boundary, centre and right boundary coordinates to solar-X displacement in arcsec
fX_pos_arcsec[1] = ( fR * ( 1.0 - cos(fs/fR) ) ) / CM_PER_ARCSEC;

fsL = fs - ( fds / 2.0 );
fX_pos_arcsec[0] = ( fR * ( 1.0 - cos(fsL/fR) ) ) / CM_PER_ARCSEC;

fsR = fsL + fds;
fX_pos_arcsec[2] = ( fR * ( 1.0 - cos(fsR/fR) ) ) / CM_PER_ARCSEC;

// Find the detector pixel within which the projected grid cell begins
iStartPixel = 0;
for( i=0; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[0] < fPixel_pos_arcsec[2] )
        break;

    iStartPixel++;
}

// Find the detector pixel within which the projected grid cell ends
iEndPixel = iStartPixel;
for( i=iStartPixel; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[2] <= fPixel_pos_arcsec[2] )
        break;

    iEndPixel++;
}

if( iStartPixel == iEndPixel )
{
    // The grid cell is entirely contained within the pixel
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fTheta = ( _PI_ * fs ) / fLength;
    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fds * cos(fTheta) * cos(fTheta) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
        ppfNEQDetector_SPEC[iStartPixel][j] += ( fTemp * fNEQCounts ) * exp( fIndex );
    }
}
else
{
    // Calculate the fraction of emission for the left-most pixel
    // Find the field-aligned coordinate 's' of the right-hand boundary of the start pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iStartPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[2]*CM_PER_ARCSEC) / fR ) );
    fdsL = fTemp - fsL;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaL = ( _PI_ * (0.5*(fsL+fTemp)) ) / fLength;
    fLOS = ( fdsL * cos(fThetaL) * cos(fThetaL) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
        ppfNEQDetector_SPEC[iStartPixel][j] += ( fTemp * fNEQCounts ) * exp( fIndex );
    }

    // Calculate the fraction of emission for the right-most pixel
    // Find the field-aligned coordinate 's' of the left-hand boundary of the end pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iEndPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[0]*CM_PER_ARCSEC) / fR ) );
    fdsR = fsR - fTemp;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaR = ( _PI_ * (0.5*(fTemp+fsR)) ) / fLength;
    fLOS = ( fdsR * cos(fThetaR) * cos(fThetaR) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iEndPixel] += ( fLOS * fNEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iEndPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
        ppfNEQDetector_SPEC[iEndPixel][j] += ( fTemp * fNEQCounts ) * exp( fIndex );
    }

    fTemp = fSD * ( fdsL / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaL) * sin(fThetaL);
    fTemp = fLOS / fSpectralProperties[2];
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iStartPixel] += ( fLOS * fNEQCounts );
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
        ppfNEQDetector_SPEC[iStartPixel][j] += ( fTemp * fNEQCounts ) * exp( fIndex );
    }

    fTemp = fSD * ( fdsR / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaR) * sin(fThetaR);
    fTemp = fLOS / fSpectralProperties[2];
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    pfNEQDetector[iEndPixel] += ( fLOS * fNEQCounts );
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iEndPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
        ppfNEQDetector_SPEC[iEndPixel][j] += ( fTemp * fNEQCounts ) * exp( fIndex );
    }

    // Calculate the fraction of emission for the remaining pixels
    iRemainingPixels = ( iEndPixel - iStartPixel - 1 );
    if( iRemainingPixels > 0 )
    {
        fdsM = ( fds - fdsL - fdsR ) / ((double)iRemainingPixels);
	fTemp = fsL + fdsL + ( fdsM / 2.0 );
        for( i=iStartPixel+1; i<iEndPixel; i++ )
	{
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
	    fTheta = ( _PI_ * fTemp ) / fLength;
	    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fdsM * cos(fTheta) * cos(fTheta) );
	    fTemp += fdsM;
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
            pfEQDetector[i] += ( fLOS * fEQCounts );
            pfNEQDetector[i] += ( fLOS * fNEQCounts );
            for( j=0; j<iNumSpectralDataPoints; j++ )
            {
                // Calculate the index of the exponential
		fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
		fIndex *= fIndex;
		fIndex /= 2.0 * fSpectralProperties[1];
		fIndex = -fIndex;
                // Calculate the spectral line intensity
                ppfEQDetector_SPEC[i][j] += ( (fLOS / fSpectralProperties[2]) * fEQCounts ) * exp( fIndex );
                ppfNEQDetector_SPEC[i][j] += ( (fLOS / fSpectralProperties[2]) * fNEQCounts ) * exp( fIndex );
            }
	}
    }
}
}

void CInstrument::Detect( PHYDATA PHYData, double fSD, double fLength, double fEQCounts, double fRestLambda, double *fSpectralProperties )
{
double fX_pos_arcsec[3], fPixel_pos_arcsec[3];
double fTheta, fThetaL, fThetaR, fLOS;
double fs, fds, fsL, fsR, fdsL, fdsR, fdsM;
double fIndex, fTemp;
int iStartPixel, iEndPixel, iRemainingPixels, i, j;

// Don't try to use the detectors if they have not been created
if( !ppfEQDetector_SPEC )
    return;

fs = PHYData.fs;
fds = PHYData.fds;

// Convert the field-aligned grid cell left boundary, centre and right boundary coordinates to solar-X displacement in arcsec
fX_pos_arcsec[1] = ( fR * ( 1.0 - cos(fs/fR) ) ) / CM_PER_ARCSEC;

fsL = fs - ( fds / 2.0 );
fX_pos_arcsec[0] = ( fR * ( 1.0 - cos(fsL/fR) ) ) / CM_PER_ARCSEC;

fsR = fsL + fds;
fX_pos_arcsec[2] = ( fR * ( 1.0 - cos(fsR/fR) ) ) / CM_PER_ARCSEC;

// Find the detector pixel within which the projected grid cell begins
iStartPixel = 0;
for( i=0; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[0] < fPixel_pos_arcsec[2] )
        break;

    iStartPixel++;
}

// Find the detector pixel within which the projected grid cell ends
iEndPixel = iStartPixel;
for( i=iStartPixel; i<iXY[0]; i++ )
{
    fPixel_pos_arcsec[1] = ppfSolarXY[0][i];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];

    if( fX_pos_arcsec[2] <= fPixel_pos_arcsec[2] )
        break;

    iEndPixel++;
}

if( iStartPixel == iEndPixel )
{
    // The grid cell is entirely contained within the pixel
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fTheta = ( _PI_ * fs ) / fLength;
    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fds * cos(fTheta) * cos(fTheta) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
    }
}
else
{
    // Calculate the fraction of emission for the left-most pixel
    // Find the field-aligned coordinate 's' of the right-hand boundary of the start pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iStartPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[2]*CM_PER_ARCSEC) / fR ) );
    fdsL = fTemp - fsL;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaL = ( _PI_ * (0.5*(fsL+fTemp)) ) / fLength;
    fLOS = ( fdsL * cos(fThetaL) * cos(fThetaL) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
    }

    // Calculate the fraction of emission for the right-most pixel
    // Find the field-aligned coordinate 's' of the left-hand boundary of the end pixel
    fPixel_pos_arcsec[1] = ppfSolarXY[0][iEndPixel];
    fPixel_pos_arcsec[0] = fPixel_pos_arcsec[1] - ( fPixel_arcsec[0] / 2.0 );
    fPixel_pos_arcsec[2] = fPixel_pos_arcsec[0] + fPixel_arcsec[0];
    fTemp = fR * acos( 1.0 - ( (fPixel_pos_arcsec[0]*CM_PER_ARCSEC) / fR ) );
    fdsR = fsR - fTemp;
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    fThetaR = ( _PI_ * (0.5*(fTemp+fsR)) ) / fLength;
    fLOS = ( fdsR * cos(fThetaR) * cos(fThetaR) );
    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    fTemp = fLOS / fSpectralProperties[2];
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iEndPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
    }

    fTemp = fSD * ( fdsL / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaL) * sin(fThetaL);
    fTemp = fLOS / fSpectralProperties[2];
    pfEQDetector[iStartPixel] += ( fLOS * fEQCounts );
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iStartPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
    }

    fTemp = fSD * ( fdsR / ( fdsL + fdsR ) );
    fLOS = fTemp * sin(fThetaR) * sin(fThetaR);
    fTemp = fLOS / fSpectralProperties[2];
    pfEQDetector[iEndPixel] += ( fLOS * fEQCounts );
    for( j=0; j<iNumSpectralDataPoints; j++ )
    {
        // Calculate the index of the exponential
        fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
	fIndex *= fIndex;
	fIndex /= 2.0 * fSpectralProperties[1];
        fIndex = -fIndex;
        // Calculate the spectral line intensity
        ppfEQDetector_SPEC[iEndPixel][j] += ( fTemp * fEQCounts ) * exp( fIndex );
    }

    // Calculate the fraction of emission for the remaining pixels
    iRemainingPixels = ( iEndPixel - iStartPixel - 1 );
    if( iRemainingPixels > 0 )
    {
        fdsM = ( fds - fdsL - fdsR ) / ((double)iRemainingPixels);
	fTemp = fsL + fdsL + ( fdsM / 2.0 );
        for( i=iStartPixel+1; i<iEndPixel; i++ )
	{
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
	    fTheta = ( _PI_ * fTemp ) / fLength;
	    fLOS = ( fSD * sin(fTheta) * sin(fTheta) ) + ( fdsM * cos(fTheta) * cos(fTheta) );
	    fTemp += fdsM;
	    // **** CALCULATE THE LINE-OF-SIGHT DEPTH ****
            pfEQDetector[i] += ( fLOS * fEQCounts );
            for( j=0; j<iNumSpectralDataPoints; j++ )
            {
                // Calculate the index of the exponential
		fIndex = ( pfLambda[j] - fRestLambda - fSpectralProperties[3] );
		fIndex *= fIndex;
		fIndex /= 2.0 * fSpectralProperties[1];
		fIndex = -fIndex;
                // Calculate the spectral line intensity
                ppfEQDetector_SPEC[i][j] += ( (fLOS / fSpectralProperties[2]) * fEQCounts ) * exp( fIndex );
            }
	}
    }
}
}

void CInstrument::ForwardModel_IMAGING( PLOOP pLoop )
{
FILE *pStrandsFile;
char szStrandsFilename[256];
double fSD, fLength, fEQCounts, fNEQCounts, fAb, *pfEQ, fNEQ, fIonEmiss, fEQElementEmiss, fNEQElementEmiss;
int iNumStrands, iStrandRange[3], iNumCells, iNumEQElements, iNumNEQElements, iZ, *piZ;
int i, j, k, m;
PHYDATA PHYData;

iNumStrands = pLoop->GetNumStrands();
pLoop->GetStrandRange( iStrandRange );

// Open a data file to store the emission as a function of position along each strand, for all strands
sprintf( szStrandsFilename, "Results/%s[%ito%i].strands", szName, iStrandRange[0], iStrandRange[1] );
pStrandsFile = fopen( szStrandsFilename, "w" );
// Write the length of the loop to the data file
fprintf( pStrandsFile, "%.8e\n", pLoop->GetLength(0) );
// Write the instrument pixel size (arcsec) to the data file
fprintf( pStrandsFile, "%g\n", sqrt( fPixel_Area_arcsec ) );
// Write the instrument cadence to the data file
fprintf( pStrandsFile, "%g\n", fCadence );
// Write the number of strands to the data file
fprintf( pStrandsFile, "%i\n", iNumStrands );

printf( "Calculating %s:\n\n", szName );
for( i=0; i<iNumStrands; i++ )
{
    // Get the number of grid cells comprising the current strand
    iNumCells = pLoop->GetNumCells( i );

    // Get the diameter of the current strand
    fSD = pLoop->GetDiameter( i );

    // Get the length of the current strand
    fLength = pLoop->GetLength( i );

    // Write the timestamp and the number of grid cells comprising the current strand to the data file
    fprintf( pStrandsFile, "%g\n%i\n", pLoop->GetTimeStamp(i), iNumCells );

    printf( "Strand %i", iStrandRange[0] + ( i * iStrandRange[2] ) );

    // Determine whether any non-equilibrium ion population data exists for the current strand
    iNumNEQElements = pLoop->GetNumNEQElements( i );

    if( iNumNEQElements > 0 )
    {
        for( j=0; j<iNumCells; j++ )
    		{
            pLoop->GetPHYData( i, j, &PHYData );

            // Reset the total counts detected from the current grid cell
            fEQCounts = 0.0;
            fNEQCounts = 0.0;

            for( k=0; k<iNumNEQElements; k++ )
            {
                // Reset the total emission due to this element (summed over all ionisation states)
				fEQElementEmiss = 0.0;
				fNEQElementEmiss = 0.0;

				// Get the atomic number of the element
				iZ = pLoop->GetAtomicNumber( i, k );

				// Get the element abundance
				fAb = pNEQRadiation->GetAbundance( iZ );

				// Allocate sufficient space to store the equilibrium ion population fractions
				pfEQ = (double*)malloc( sizeof(double) * (iZ+1) );
#ifdef DENSITY_DEPENDENT_RATES
				pNEQRadiation->GetEquilIonFrac( iZ, pfEQ, log10(PHYData.fTe), log10(PHYData.fne) );
#else // DENSITY_DEPENDENT_RATES
				pNEQRadiation->GetEquilIonFrac( iZ, pfEQ, log10(PHYData.fTe) );
#endif // DENSITY_DEPENDENT_RATES

				for( m=0; m<=iZ; m++ )
				{
                    // Get the ion population fraction
                    fNEQ = pLoop->GetNEQ( i, j, k, m );
                    // Get the ion emissivity
                    fIonEmiss = GetIonEmission( iZ, m+1, log10(PHYData.fTe), log10(PHYData.fne) );
                    // Multiply the stored ion emissivity by the ion population fraction and sum
                    fEQElementEmiss += ( pfEQ[m] * fIonEmiss );
                    fNEQElementEmiss += ( fNEQ * fIonEmiss );
				}

				// Multiply by the element abundance relative to hydrogen and the column emission measure, and divide by the number of strands (to avoid summing over time)
				fEQElementEmiss *= ( fAb * PHYData.fne * PHYData.fnH ) / ((double)iNumStrands);
				fNEQElementEmiss *= ( fAb * PHYData.fne * PHYData.fnH ) / ((double)iNumStrands);

                // Keep track of the total counts detected from the current grid cell
                fEQCounts += fEQElementEmiss;
                fNEQCounts += fNEQElementEmiss;

				free( pfEQ );
            }
			// Record the specified counts [DN pixel^-1 s^-1] at the appropriate detector pixel, taking into account the line-of-sight depth
			Detect( PHYData.fs, PHYData.fds, fSD, fLength, fEQCounts, fNEQCounts );
            // Write the equilibrium and non-equilibrium counts to the strands data file [DN pixel^-1 s^-1]
            fprintf( pStrandsFile, "%.8e\t%.8e\t%.8e\t%.8e\n", PHYData.fs, PHYData.fds, ( PHYData.fds * fEQCounts ), ( PHYData.fds * fNEQCounts ) );
		}
		printf( " - DONE!\n" );
    }
    else
    {
        // Assume an equilibrium ionisaton state
		piZ = pEQRadiation->pGetAtomicNumbers( &iNumEQElements );

		for( j=0; j<iNumCells; j++ )
		{
            pLoop->GetPHYData( i, j, &PHYData );

            // Reset the total counts detected from the current grid cell
            fEQCounts = 0.0;

            for( k=0; k<iNumEQElements; k++ )
            {
                // Reset the total emission due to this element (summed over all ionisation states)
				fEQElementEmiss = 0.0;

				// Get the element abundance
				fAb = pEQRadiation->GetAbundance( piZ[k] );

				// Allocate sufficient space to store the equilibrium ion population fractions
				pfEQ = (double*)malloc( sizeof(double) * (piZ[k]+1) );
#ifdef DENSITY_DEPENDENT_RATES
				pEQRadiation->GetEquilIonFrac( piZ[k], pfEQ, log10(PHYData.fTe), log10(PHYData.fne) );
#else // DENSITY_DEPENDENT_RATES
				pEQRadiation->GetEquilIonFrac( piZ[k], pfEQ, log10(PHYData.fTe) );
#endif // DENSITY_DEPENDENT_RATES

				for( m=0; m<=piZ[k]; m++ )
				{
                    // Multiply the stored ion emissivity by the ion population fraction and sum
                    fEQElementEmiss += ( pfEQ[m] * GetIonEmission( piZ[k], m+1, log10(PHYData.fTe), log10(PHYData.fne) ) );
				}

				// Multiply by the element abundance relative to hydrogen and the column emission measure
				fEQElementEmiss *= ( fAb * PHYData.fne * PHYData.fnH ) / ((double)iNumStrands);

	          	// Keep track of the total counts detected from the current grid cell
                fEQCounts += fEQElementEmiss;

				free( pfEQ );
            }
			// Record the specified counts [DN pixel^-1 s^-1] at the appropriate detector pixel, taking into account the line-of-sight depth
			Detect( PHYData.fs, PHYData.fds, fSD, fLength, fEQCounts );
            // Write the equilibrium counts to the strands data file [DN pixel^-1 s^-1]
            fprintf( pStrandsFile, "%.8e\t%.8e\t%.8e\n", PHYData.fs, PHYData.fds, ( PHYData.fds * fEQCounts ) );
		}
		printf( " - DONE!\n" );
	}
}

// Close the strands data file
fclose( pStrandsFile );

printf( "\n" );
}

void CInstrument::ForwardModel_SPECTROSCOPIC( PLOOP pLoop )
{
FILE *pStrandsFile;
char szStrandsFilename[256];
double fSD, fLength, fEQCounts, fNEQCounts, fAb, fConst, fvi_th, *pfEQ, fNEQ, *pfLineList, fSpectralProperties[4], fLineEmiss;
int iNumStrands, iStrandRange[3], iNumCells, iNumEQElements, iNumNEQElements, iZ, *piZ, iNumLines;
// i = strand no., j = cell no., k = element no., m = ion no., n = emission line no.
int i, j, k, m, n;
PHYDATA PHYData;

iNumStrands = pLoop->GetNumStrands();
pLoop->GetStrandRange( iStrandRange );

// Open a data file to store the emission as a function of position along each strand, for all strands
sprintf( szStrandsFilename, "Results/%s[%ito%i].strands", szName, iStrandRange[0], iStrandRange[1] );
pStrandsFile = fopen( szStrandsFilename, "w" );
// Write the length of the loop to the data file
fprintf( pStrandsFile, "%.8e\n", pLoop->GetLength(0) );
// Write the instrument pixel size (arcsec) to the data file
fprintf( pStrandsFile, "%g\n", sqrt( fPixel_Area_arcsec ) );
// Write the instrument cadence to the data file
fprintf( pStrandsFile, "%g\n", fCadence );
// Write the number of strands to the data file
fprintf( pStrandsFile, "%i\n", iNumStrands );

printf( "Calculating %s:\n\n", szName );
for( i=0; i<iNumStrands; i++ )
{
    // Get the number of grid cells comprising the current strand
    iNumCells = pLoop->GetNumCells( i );

    // Get the diameter of the current strand
    fSD = pLoop->GetDiameter( i );

    // Get the length of the current strand
    fLength = pLoop->GetLength( i );

    // Write the timestamp and the number of grid cells comprising the current strand to the data file
    fprintf( pStrandsFile, "%g\n%i\n", pLoop->GetTimeStamp(i), iNumCells );

    printf( "Strand %i", iStrandRange[0] + ( i * iStrandRange[2] ) );

    // Determine whether any non-equilibrium ion population data exists for the current strand
    iNumNEQElements = pLoop->GetNumNEQElements( i );

    if( iNumNEQElements > 0 )
    {
        for( j=0; j<iNumCells; j++ )
		{
            pLoop->GetPHYData( i, j, &PHYData );

            // Reset the total counts detected from the current grid cell
            fEQCounts = 0.0;
            fNEQCounts = 0.0;

            for( k=0; k<iNumNEQElements; k++ )
            {
                // Get the atomic number of the element
				iZ = pLoop->GetAtomicNumber( i, k );

				// Get the element abundance
				fAb = pNEQRadiation->GetAbundance( iZ );

				// Allocate sufficient space to store the equilibrium ion population fractions
				pfEQ = (double*)malloc( sizeof(double) * (iZ+1) );
#ifdef DENSITY_DEPENDENT_RATES
				pNEQRadiation->GetEquilIonFrac( iZ, pfEQ, log10(PHYData.fTe), log10(PHYData.fne) );
#else // DENSITY_DEPENDENT_RATES
				pNEQRadiation->GetEquilIonFrac( iZ, pfEQ, log10(PHYData.fTe) );
#endif // DENSITY_DEPENDENT_RATES

                // Calculate the multiplicative factor used to convert to [DN pixel^-1 s^-1]
                fConst = fAb * PHYData.fne * PHYData.fnH;

                // Calculate the most probable ion speed (page 221, Dere, K. P., & Mason, H. E. 1993, Sol. Phys., 144, 217)
                fvi_th = sqrt( (2.0*BOLTZMANN_CONSTANT*PHYData.fTi) / GetMass( iZ ) );

				for( m=0; m<=iZ; m++ )
				{
                    // Get the ion population fraction
                    fNEQ = pLoop->GetNEQ( i, j, k, m );

                    // Get the number of emission lines for the current ion
                    iNumLines = GetNumLines( iZ, m+1 );

                    // Allocate sufficient space to store the wavelength list for the current ion
                    pfLineList = (double*)malloc( sizeof(double) * iNumLines );

                    // Get the wavelength list for the current ion
                    GetLineList( iZ, m+1, pfLineList );

                    for( n=0; n<iNumLines; n++ )
                    {
                       	// Calculate the spectral properties of the line
						// (Eq. (2), Dere, K. P., & Mason, H. E. 1993, Sol. Phys., 144, 217)
                        fSpectralProperties[0] = ( pfLineList[n] * pfLineList[n] ) / ( 2.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT );
						// Calculate the Doppler width plus instrumental width (the instrumental width is the FWHM so needs to be converted into Gaussian width)
						// FWHM = 2*SQRT(2*LN(2)) * Gaussian Width -> 1.0/8*LN(2) = 0.18033688
						fSpectralProperties[1] = ( fSpectralProperties[0] * fvi_th * fvi_th ) + ( 0.18033688 * fInstrumentalLineWidth * fInstrumentalLineWidth );
						// SQRT( 2*PI ) = 2.50662827
						fSpectralProperties[2] = 2.50662827 * sqrt( fSpectralProperties[1] );
						fSpectralProperties[3] = sqrt( 2.0 * fSpectralProperties[0] ) * PHYData.fvp;
                        
						// Get the line emissivity
                        fLineEmiss = GetLineEmission( iZ, m+1, pfLineList[n], log10(PHYData.fTe), log10(PHYData.fne) );
                        // Multiply the stored line emissivity by the equilibrium and non-equilibrium ion population fractions, the element abundance relative to hydrogen, the density squared,
						// and divide by the number of strands (to avoid summing over time)
                        fLineEmiss *= fConst / ((double)iNumStrands);
                        Detect( PHYData, fSD, fLength, fLineEmiss*pfEQ[m], fLineEmiss*fNEQ, pfLineList[n], fSpectralProperties );

                        // Keep track of the total counts detected from the current grid cell
                        fEQCounts += fLineEmiss*pfEQ[m];
                        fNEQCounts += fLineEmiss*fNEQ;
                    }
                    free( pfLineList );
				}
				free( pfEQ );
            }
            // Write the equilibrium and non-equilibrium counts to the strands data file [DN pixel^-1 s^-1]
            fprintf( pStrandsFile, "%.8e\t%.8e\t%.8e\t%.8e\n", PHYData.fs, PHYData.fds, ( PHYData.fds * fEQCounts ), ( PHYData.fds * fNEQCounts ) );
		}
        printf( " - DONE!\n" );
    }
    else
    {
        // Assume an equilibrium ionisation state
		piZ = pEQRadiation->pGetAtomicNumbers( &iNumEQElements );

		for( j=0; j<iNumCells; j++ )
		{
            pLoop->GetPHYData( i, j, &PHYData );

            // Reset the total counts detected from the current grid cell
            fEQCounts = 0.0;

            for( k=0; k<iNumEQElements; k++ )
            {
                // Get the element abundance
				fAb = pEQRadiation->GetAbundance( piZ[k] );

				// Allocate sufficient space to store the equilibrium ion population fractions
				pfEQ = (double*)malloc( sizeof(double) * (piZ[k]+1) );
#ifdef DENSITY_DEPENDENT_RATES
				pEQRadiation->GetEquilIonFrac( piZ[k], pfEQ, log10(PHYData.fTe), log10(PHYData.fne) );
#else // DENSITY_DEPENDENT_RATES
				pEQRadiation->GetEquilIonFrac( piZ[k], pfEQ, log10(PHYData.fTe) );
#endif // DENSITY_DEPENDENT_RATES

                // Calculate the multiplicative factor used to convert to [DN pixel^-1 s^-1]
                fConst = fAb * PHYData.fne * PHYData.fnH;

                // Calculate the most probable ion speed (page 221, Dere, K. P., & Mason, H. E. 1993, Sol. Phys., 144, 217)
                fvi_th = sqrt( (2.0*BOLTZMANN_CONSTANT*PHYData.fTi) / GetMass( piZ[k] ) );

				for( m=0; m<=piZ[k]; m++ )
				{
                    // Get the number of emission lines for the current ion
                    iNumLines = GetNumLines( piZ[k], m+1 );

                    // Allocate sufficient space to store the wavelength list for the current ion
                    pfLineList = (double*)malloc( sizeof(double) * iNumLines );

                    // Get the wavelength list for the current ion
                    GetLineList( piZ[k], m+1, pfLineList );

                    for( n=0; n<iNumLines; n++ )
                    {
                        // Calculate the spectral properties of the line
						// (Eq. (2), Dere, K. P., & Mason, H. E. 1993, Sol. Phys., 144, 217)
                        fSpectralProperties[0] = ( pfLineList[n] * pfLineList[n] ) / ( 2.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT );
						// Calculate the Doppler width plus instrumental width (the instrumental width is the FWHM so needs to be converted into Gaussian width)
						// FWHM = 2*SQRT(2*LN(2)) * Gaussian Width -> 1.0/8*LN(2) = 0.18033688
						fSpectralProperties[1] = ( fSpectralProperties[0] * fvi_th * fvi_th ) + ( 0.18033688 * fInstrumentalLineWidth * fInstrumentalLineWidth );
						// SQRT( 2*PI ) = 2.50662827
						fSpectralProperties[2] = 2.50662827 * sqrt( fSpectralProperties[1] );
						fSpectralProperties[3] = sqrt( 2.0 * fSpectralProperties[0] ) * PHYData.fvp;

                        // Get the line emissivity
                        fLineEmiss = GetLineEmission( piZ[k], m+1, pfLineList[n], log10(PHYData.fTe), log10(PHYData.fne) );
                        // Multiply the stored line emissivity by the equilibrium and non-equilibrium ion population fractions, the element abundance relative to hydrogen, the density squared,
                        // and divide by the number of strands (to avoid summing over time)
						fLineEmiss *= fConst / ((double)iNumStrands);
                        Detect( PHYData, fSD, fLength, fLineEmiss*pfEQ[m], pfLineList[n], fSpectralProperties );

                        // Keep track of the total counts detected from the current grid cell
                        fEQCounts += fLineEmiss*pfEQ[m];
                    }
                    free( pfLineList );
				}
				free( pfEQ );
            }
            // Write the equilibrium and non-equilibrium counts to the strands data file [DN pixel^-1 s^-1]
            fprintf( pStrandsFile, "%.8e\t%.8e\t%.8e\n", PHYData.fs, PHYData.fds, ( PHYData.fds * fEQCounts ) );
		}
		printf( " - DONE!\n" );
    }
}

// Close the strands data file
fclose( pStrandsFile );

printf( "\n" );
}

void CInstrument::GetInstrumentName( char *pszName )
{
strcpy( pszName, szName );
}

void CInstrument::WriteIonLineListToFile( FILE *pFile )
{
int i;

fprintf( pFile, "Instrument: %s (%.3f to %.3f)\n", szName, fLambda_min, fLambda_max );

for( i=0; i<iNumIons; i++ )
    ppIon[i]->WriteIonLineListToFile( pFile );

fprintf( pFile, "\n" );
}

void CInstrument::CreateVirtualDetector( PLOOP pLoop )
{
int i, j, iNumStrands;
bool iNEQ = false;

// If a virtual detector has already been created then delete it
if( ppfSolarXY )
    DeleteVirtualDetector();

// Determine whether any non-equilibrium ion population data exists for the current loop
iNumStrands = pLoop->GetNumStrands();
for( i=0; i<iNumStrands; i++ )
    if( pLoop->GetNumNEQElements( i ) > 0 )
    {
        iNEQ = true;
        break;
    }

// Assume the loop length is equal to the length of strand 0
fR = (pLoop->GetLength( 0 )) / _PI_;

// Calculate the number of detector pixels in the solar-X and Y directions
for( i=0; i<2; i++ )
	iXY[i] = (int)( ( ( 2.0 * fR ) / fPixel_cm[i] ) + 1.0 );

// Calculate the solar-X and Y coordinates in arcsec at centre of each pixel
ppfSolarXY = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
{
    ppfSolarXY[i] = (double*)malloc( sizeof(double) * iXY[i] );
    for( j=0; j<iXY[i]; j++ )
        ppfSolarXY[i][j] = ( (double)j + 0.5 ) * fPixel_arcsec[i];
}

// **** IMAGING DETECTOR ****
// Create a 1D detector initially, with pixels in the solar-X direction
pfEQDetector = (double*)malloc( sizeof(double) * iXY[0] );

if( iNEQ )
    pfNEQDetector = (double*)malloc( sizeof(double) * iXY[0] );
// **** IMAGING DETECTOR ****

// **** SPECTROSCOPIC DETECTOR ****
if( fSpecRes )
{
    // Calculate the number of data points needed to plot the spectrum for each pixel
    iNumSpectralDataPoints = (int)( ( fLambda_max - fLambda_min ) / fSpecRes );

    // Calculate the wavelength scale in angstroms
    pfLambda = (double*)malloc( sizeof(double) * iNumSpectralDataPoints );
    for( i=0; i<iNumSpectralDataPoints; i++ )
        pfLambda[i] = fLambda_min + ( ((double)i) * fSpecRes );

    // Create a 1D detector initially, with pixels in the solar-X direction, and a spectrum associated with each pixel
    ppfEQDetector_SPEC = (double**)malloc( sizeof(double*) * iXY[0] );
    for( i=0; i<iXY[0]; i++ )
        ppfEQDetector_SPEC[i] = (double*)malloc( sizeof(double) * iNumSpectralDataPoints );

    if( iNEQ )
    {
        // Create a 1D detector initially, with pixels in the solar-X direction, and a spectrum associated with each pixel
        ppfNEQDetector_SPEC = (double**)malloc( sizeof(double*) * iXY[0] );
        for( i=0; i<iXY[0]; i++ )
            ppfNEQDetector_SPEC[i] = (double*)malloc( sizeof(double) * iNumSpectralDataPoints );
    }
}
// **** SPECTROSCOPIC DETECTOR ****

// Reset the virtual detector
ResetVirtualDetector();
}

bool CInstrument::IsEquilibriumDetector( void )
{
if( pfEQDetector ) return true;
    
return false;
}

bool CInstrument::IsNonEquilibriumDetector( void )
{
if( pfNEQDetector ) return true;

return false;
}

void CInstrument::DeleteVirtualDetector( void )
{
int i;

if( ppfSolarXY )
{
    if( fSpecRes )
    {
        if( ppfNEQDetector_SPEC )
        {
            for( i=0; i<iXY[0]; i++ )
                free( ppfNEQDetector_SPEC[i] );

            free( ppfNEQDetector_SPEC );
        }

        for( i=0; i<iXY[0]; i++ )
            free( ppfEQDetector_SPEC[i] );

        free( ppfEQDetector_SPEC );

        free( pfLambda );
    }

    if( pfNEQDetector )
        free( pfNEQDetector );

    free( pfEQDetector );

    for( i=0; i<2; i++ )
        free( ppfSolarXY[i] );

    free( ppfSolarXY );
}

ppfSolarXY = NULL;
pfEQDetector = NULL;
pfNEQDetector = NULL;
pfLambda = NULL;
ppfEQDetector_SPEC = NULL;
ppfNEQDetector_SPEC = NULL;
}

void CInstrument::ResetVirtualDetector( void )
{
int i;

if( pfEQDetector )
    memset( pfEQDetector, 0, sizeof(double) * iXY[0] );

if( pfNEQDetector )
    memset( pfNEQDetector, 0, sizeof(double) * iXY[0] );

if( fSpecRes )
{
    for( i=0; i<iXY[0]; i++ )
        memset( ppfEQDetector_SPEC[i], 0, sizeof(double) * iNumSpectralDataPoints );

    if( ppfNEQDetector_SPEC )
    {
        for( i=0; i<iXY[0]; i++ )
            memset( ppfNEQDetector_SPEC[i], 0, sizeof(double) * iNumSpectralDataPoints );
    }
}
}

void CInstrument::WriteDetectorToFile( char *pszWorkingDirectory )
{
FILE *pFile;
char szFilename[256];
int i, j;

if( !ppfSolarXY )
    return;

sprintf( szFilename, "%s/%s.det", pszWorkingDirectory, szName );
printf( "%s", szFilename );

pFile = fopen( szFilename, "w" );

for( i=0; i<iXY[0]; i++ )
{
    // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
    // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
    fprintf( pFile, "%.8e\t%.8e", ppfSolarXY[0][i], pfEQDetector[i] );
    if( pfNEQDetector )
        fprintf( pFile, "\t%.8e",  pfNEQDetector[i] );
    fprintf( pFile, "\n" );
}

fclose( pFile );

if( fSpecRes )
{
    sprintf( szFilename, "%s/%s.spec.det", pszWorkingDirectory, szName );
    printf( ";\t%s", szFilename );

    // FILE STRUCTURE:
    // Number of spectral data points
    // Wavelength values
    // Solar-X
    // Equilibrium ionisation spectrum
    // Non-equilibrium ionisation spectrum
    // Solar-X
    // ..
    // ..

    pFile = fopen( szFilename, "w" );

    // Write the wavelength scale into the data file
    fprintf( pFile, "%i\n", iNumSpectralDataPoints );
    fprintf( pFile, "%.8e", pfLambda[0] );
    for( i=1; i<iNumSpectralDataPoints; i++ )
        fprintf( pFile, "\t%.8e", pfLambda[i] );
    fprintf( pFile, "\n" );

    for( i=0; i<iXY[0]; i++ )
    {
        // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
        // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
	fprintf( pFile, "%.8e\n%.8e", ppfSolarXY[0][i], ppfEQDetector_SPEC[i][0] );
        for( j=1; j<iNumSpectralDataPoints; j++ )
            fprintf( pFile, "\t%.8e", ppfEQDetector_SPEC[i][j] );
        fprintf( pFile, "\n" );

        if( ppfNEQDetector_SPEC )
        {
            fprintf( pFile, "%.8e",  ppfNEQDetector_SPEC[i][0] );
            for( j=1; j<iNumSpectralDataPoints; j++ )
                fprintf( pFile, "\t%.8e",  ppfNEQDetector_SPEC[i][j] );
            fprintf( pFile, "\n" );
        }
    }

    fclose( pFile );
}

printf( " - DONE!\n" );
}

void CInstrument::WriteDetectorToFile( char *pszWorkingDirectory, int iFileNumber )
{
FILE *pFile;
char szFilename[256];
int i, j;

if( !ppfSolarXY )
    return;

sprintf( szFilename, "%s/%s.%i.det", pszWorkingDirectory, szName, iFileNumber );
printf( "%s", szFilename );

pFile = fopen( szFilename, "w" );

for( i=0; i<iXY[0]; i++ )
{
    // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
    // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
    fprintf( pFile, "%.8e\t%.8e", ppfSolarXY[0][i], pfEQDetector[i] );
    if( pfNEQDetector )
        fprintf( pFile, "\t%.8e",  pfNEQDetector[i] );
    fprintf( pFile, "\n" );
}

fclose( pFile );

if( fSpecRes )
{
    sprintf( szFilename, "%s/%s.%i.spec.det", pszWorkingDirectory, szName, iFileNumber );
    printf( ";\t%s", szFilename );

    // FILE STRUCTURE:
    // Number of spectral data points
    // Wavelength values
    // Solar-X
    // Equilibrium ionisation spectrum
    // Non-equilibrium ionisation spectrum
    // Solar-X
    // ..
    // ..

    pFile = fopen( szFilename, "w" );

    // Write the wavelength scale into the data file
    fprintf( pFile, "%i\n", iNumSpectralDataPoints );
    fprintf( pFile, "%.8e", pfLambda[0] );
    for( i=1; i<iNumSpectralDataPoints; i++ )
        fprintf( pFile, "\t%.8e", pfLambda[i] );
    fprintf( pFile, "\n" );

    for( i=0; i<iXY[0]; i++ )
    {
        // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
        // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
	fprintf( pFile, "%.8e\n%.8e", ppfSolarXY[0][i], ppfEQDetector_SPEC[i][0] );
        for( j=1; j<iNumSpectralDataPoints; j++ )
            fprintf( pFile, "\t%.8e", ppfEQDetector_SPEC[i][j] );
        fprintf( pFile, "\n" );

        if( ppfNEQDetector_SPEC )
        {
            fprintf( pFile, "%.8e",  ppfNEQDetector_SPEC[i][0] );
            for( j=1; j<iNumSpectralDataPoints; j++ )
                fprintf( pFile, "\t%.8e",  ppfNEQDetector_SPEC[i][j] );
            fprintf( pFile, "\n" );
        }
    }

    fclose( pFile );
}

printf( " - DONE!\n" );
}

void CInstrument::ForwardModel( PLOOP pLoop )
{
    if( fSpecRes )
        ForwardModel_SPECTROSCOPIC( pLoop );
    else
        ForwardModel_IMAGING( pLoop );
}

void CInstrument::EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory )
{
FILE *pFile;
PRADIATION pRadiation;
char szFilename[256];
double fAb, flog10_T, fLineEmiss, *pfI_per_EMc;
double *pfLineList, *pfEQ;
int iNumLines = 0, iNumElements, iNumVals;
int *piZ;
int i, k, m, j, n;

// Check to see whether the virtual detectors have been created
if( !ppfSolarXY )
    return;

sprintf( szFilename, "%s/%s.EM_Loci", pszWorkingDirectory, szName );
printf( "%s", szFilename );

pFile = fopen( szFilename, "w" );

// Write the temperature range
iNumVals = (int)( ( ( flog10_Tmax - flog10_Tmin ) / flog10_step ) + 1.0 );
fprintf( pFile, "%i\n", iNumVals );
for( i=0; i<iNumVals; i++ )
{
    flog10_T = flog10_Tmin+(((double)i)*flog10_step);
    fprintf( pFile, "\t%g", flog10_T );
}
fprintf( pFile, "\n" );

// Allocate sufficient space to store the equilibrium intensity per column emission measure for each temperature value
pfI_per_EMc = (double*)malloc( sizeof(double) * iNumVals );
memset( pfI_per_EMc, 0, sizeof(double) * iNumVals );

// Find out whether non-equilibrium ion populations exist
if( pfNEQDetector )
    pRadiation = pNEQRadiation;
else
    pRadiation = pEQRadiation;

piZ = pRadiation->pGetAtomicNumbers( &iNumElements );

for( k=0; k<iNumElements; k++ )
{
    // Get the element abundance
    fAb = pRadiation->GetAbundance( piZ[k] );

    // Allocate sufficient space to store the equilibrium ion population fractions
    pfEQ = (double*)malloc( sizeof(double) * (piZ[k]+1) );

    for( m=0; m<=piZ[k]; m++ )
    {
        // Get the number of emission lines for the current ion
        iNumLines = GetNumLines( piZ[k], m+1 );

        // The spectroscopic number of the emitting ion is m+1 and its array index is m
        // Allocate sufficient space to store the wavelength list for the emitting ion
        pfLineList = (double*)malloc( sizeof(double) * iNumLines );

        // Get the wavelength list for the emitting ion
        GetLineList( piZ[k], m+1, pfLineList );

        // Write the equilibrium ionization emission measure loci
        // Loop through the temperature range
        for( j=0; j<iNumVals; j++ )
        {
            flog10_T = flog10_Tmin+(((double)j)*flog10_step);
            // Get the ionization fraction for the emitting ion
#ifdef DENSITY_DEPENDENT_RATES
	    pRadiation->GetEquilIonFrac( piZ[k], pfEQ, flog10_T, flog10_n );
#else // DENSITY_DEPENDENT_RATES
	    pRadiation->GetEquilIonFrac( piZ[k], pfEQ, flog10_T );
#endif // DENSITY_DEPENDENT_RATES
	    fLineEmiss = 0.0;
            for( n=0; n<iNumLines; n++ )
            {
                // Get the total line emissivity
                fLineEmiss += GetLineEmission( piZ[k], m+1, pfLineList[n], flog10_T, flog10_n );
            }
            pfI_per_EMc[j] += ( fAb * pfEQ[m] * fLineEmiss );
        }
        free( pfLineList );
    }
    free( pfEQ );
}

// Write the rest of the data file
for( i=0; i<iXY[0]; i++ )
{
    // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
    // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
    fprintf( pFile, "%.8e\n", ppfSolarXY[0][i] );

    // These quantities are the true (model) emission measure
    for( j=0; j<iNumVals; j++ )
    {
            if( !pfI_per_EMc[j] )
                fprintf( pFile, "\t0.0" );
            else
                fprintf( pFile, "\t%.8e", pfEQDetector[i] / pfI_per_EMc[j] );
    }
    fprintf( pFile, "\n" );

    if( pfNEQDetector )
    {
        // These quantities are the emission measure shifted according to the non-equilibrium ionization state
        for( j=0; j<iNumVals; j++ )
        {
            if( !pfI_per_EMc[j] )
                fprintf( pFile, "\t0.0" );
            else
                fprintf( pFile, "\t%.8e", pfNEQDetector[i] / pfI_per_EMc[j] );
        }
        fprintf( pFile, "\n" );
    }
}

free( pfI_per_EMc );

fclose( pFile );

printf( " - DONE!\n" );
}

void CInstrument::EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory, int iFileNumber )
{
FILE *pFile;
PRADIATION pRadiation;
char szFilename[256];
double fAb, flog10_T, fLineEmiss, *pfI_per_EMc;
double *pfLineList, *pfEQ;
int iNumLines = 0, iNumElements, iNumVals;
int *piZ;
int i, k, m, j, n;

// Check to see whether the virtual detectors have been created
if( !ppfSolarXY )
    return;

sprintf( szFilename, "%s/%s.%i.EM_Loci", pszWorkingDirectory, szName, iFileNumber );
printf( "%s", szFilename );

pFile = fopen( szFilename, "w" );

// Write the temperature range
iNumVals = (int)( ( ( flog10_Tmax - flog10_Tmin ) / flog10_step ) + 1.0 );
fprintf( pFile, "%i\n", iNumVals );
for( i=0; i<iNumVals; i++ )
{
    flog10_T = flog10_Tmin+(((double)i)*flog10_step);
    fprintf( pFile, "\t%g", flog10_T );
}
fprintf( pFile, "\n" );

// Allocate sufficient space to store the equilibrium intensity per column emission measure for each temperature value
pfI_per_EMc = (double*)malloc( sizeof(double) * iNumVals );
memset( pfI_per_EMc, 0, sizeof(double) * iNumVals );

// Find out whether non-equilibrium ion populations exist
if( pfNEQDetector )
    pRadiation = pNEQRadiation;
else
    pRadiation = pEQRadiation;

piZ = pRadiation->pGetAtomicNumbers( &iNumElements );

for( k=0; k<iNumElements; k++ )
{
    // Get the element abundance
    fAb = pRadiation->GetAbundance( piZ[k] );

    // Allocate sufficient space to store the equilibrium ion population fractions
    pfEQ = (double*)malloc( sizeof(double) * (piZ[k]+1) );

    for( m=0; m<=piZ[k]; m++ )
    {
        // Get the number of emission lines for the current ion
        iNumLines = GetNumLines( piZ[k], m+1 );

        // The spectroscopic number of the emitting ion is m+1 and its array index is m
        // Allocate sufficient space to store the wavelength list for the emitting ion
        pfLineList = (double*)malloc( sizeof(double) * iNumLines );

        // Get the wavelength list for the emitting ion
        GetLineList( piZ[k], m+1, pfLineList );

        // Write the equilibrium ionization emission measure loci
        // Loop through the temperature range
        for( j=0; j<iNumVals; j++ )
        {
            flog10_T = flog10_Tmin+(((double)j)*flog10_step);
            // Get the ionization fraction for the emitting ion
#ifdef DENSITY_DEPENDENT_RATES
	    pRadiation->GetEquilIonFrac( piZ[k], pfEQ, flog10_T, flog10_n );
#else // DENSITY_DEPENDENT_RATES
	    pRadiation->GetEquilIonFrac( piZ[k], pfEQ, flog10_T );
#endif // DENSITY_DEPENDENT_RATES
	    fLineEmiss = 0.0;
            for( n=0; n<iNumLines; n++ )
            {
                // Get the total line emissivity
                fLineEmiss += GetLineEmission( piZ[k], m+1, pfLineList[n], flog10_T, flog10_n );
            }
            pfI_per_EMc[j] += ( fAb * pfEQ[m] * fLineEmiss );
        }
        free( pfLineList );
    }
    free( pfEQ );
}

// Write the rest of the data file
for( i=0; i<iXY[0]; i++ )
{
    // ppfSolarXY[0][i] = Solar-X coordinates (arcsec)
    // ppfSolarXY[1][i] = Solar-Y coordinates (arcsec)
    fprintf( pFile, "%.8e\n", ppfSolarXY[0][i] );

    // These quantities are the true (model) emission measure
    for( j=0; j<iNumVals; j++ )
    {
            if( !pfI_per_EMc[j] )
                fprintf( pFile, "\t0.0" );
            else
                fprintf( pFile, "\t%.8e", pfEQDetector[i] / pfI_per_EMc[j] );
    }
    fprintf( pFile, "\n" );

    if( pfNEQDetector )
    {
        // These quantities are the emission measure shifted according to the non-equilibrium ionization state
        for( j=0; j<iNumVals; j++ )
        {
            if( !pfI_per_EMc[j] )
                fprintf( pFile, "\t0.0" );
            else
                fprintf( pFile, "\t%.8e", pfNEQDetector[i] / pfI_per_EMc[j] );
        }
        fprintf( pFile, "\n" );
    }
}

free( pfI_per_EMc );

fclose( pFile );

printf( " - DONE!\n" );
}