// ****
// *
// * A heating code to simulate different forms of heat deposition
// *
// * Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *     
// * Date last modified: 09/12/2018
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "heat.h"
#include "../../Resources/source/gammabeta.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"


CHeat::CHeat( void )
{
Initialise();
}

CHeat::~CHeat( void )
{
FreeAll();
}

void CHeat::Initialise( void )
{
GetHeatingData();
#ifdef BEAM_HEATING
	GetBeamHeatingData();
#endif // BEAM_HEATING
#ifdef OPTICALLY_THICK_RADIATION
	GetVALHeatingData();
#endif // OPTICALLY_THICK_RADIATION
}

void CHeat::FreeAll( void )
{
#ifdef OPTICALLY_THICK_RADIATION
int i;

for( i=0; i<2; i++ )
    free( ppVALHeating[i] );
free( ppVALHeating );
#endif // OPTICALLY_THICK_RADIATION

#ifdef BEAM_HEATING
if( iBeamHeatingDP )
{
	free( pfBeamSpectralIndex );
	free( pfBeamCutOff );
	free( pfBeamEnergyFlux );
	free( pfBeamTime );
}
#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		if( pfQbeam )
			delete[] pfQbeam;
	#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // BEAM_HEATING

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

void CHeat::GetHeatingData( void )
{
FILE *pConfigFile;
int i;

// Open and read the configuration file
pConfigFile = fopen( "Heating_Model/config/heating_model.cfg", "r" );

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

#ifdef BEAM_HEATING
void CHeat::GetBeamHeatingData( void )
{
FILE *pConfigFile;
char szBuffer[256];
int i;

// Open and read the configuration file
pConfigFile = fopen( "Heating_Model/config/beam_heating_model.cfg", "r" );

// Get the number of tabulated values for the beam heating parameters
fscanf( pConfigFile, "%i", &iBeamHeatingDP );
if( !iBeamHeatingDP ) return;	// There is no beam heating in the current run

// Get the header information
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );

pfBeamTime = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamEnergyFlux = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamCutOff = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamSpectralIndex = (double*)malloc( sizeof(double) * iBeamHeatingDP );

if( iBeamHeatingDP == 1 )
{
	// The beam heating parameters are time-independent
	// The first time value in the beam configuration file will be the beam duration in this case
	ReadDouble( pConfigFile, &(pfBeamTime[0]) );
	// Get the remaining beam parameter values
	ReadDouble( pConfigFile, &(pfBeamEnergyFlux[0]) );
	ReadDouble( pConfigFile, &(pfBeamCutOff[0]) );
	ReadDouble( pConfigFile, &(pfBeamSpectralIndex[0]) );
} else {
	// The beam heating parameters are time-dependent
	// Get the tabulated beam parameter values
	for( i=0; i<iBeamHeatingDP; i++ )
	{
		ReadDouble( pConfigFile, &(pfBeamTime[i]) );
		// Get the remaining beam parameter values
		ReadDouble( pConfigFile, &(pfBeamEnergyFlux[i]) );
		ReadDouble( pConfigFile, &(pfBeamCutOff[i]) );
		ReadDouble( pConfigFile, &(pfBeamSpectralIndex[i]) );
	}
}

fclose( pConfigFile );
}
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
void CHeat::GetVALHeatingData( void )
{
FILE *pFile;
int i;

// Open and read the configuration file
pFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat", "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iVALHeatingDP );

// Allocate sufficient memory to hold the heating data
ppVALHeating = (double**)malloc( sizeof(double*) * 2 );
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
#endif // OPTICALLY_THICK_RADIATION

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

#ifdef BEAM_HEATING
void CHeat::CalculateBeamParameters( double t, double *pBeamParams )
{
double x[3], y[3];
int i;

if( !iBeamHeatingDP ) return;
else if( iBeamHeatingDP == 1 ) {
	if( t > pfBeamTime[0] ) return;
	pBeamParams[0] = pfBeamEnergyFlux[0]; 
	pBeamParams[1] = pfBeamCutOff[0];
	pBeamParams[2] = pfBeamSpectralIndex[0];
} else {
	if( t < pfBeamTime[0] || t > pfBeamTime[iBeamHeatingDP-1] ) return;

	// Find the time interval for the interpolation of the parameter values
	for( i=0; i<iBeamHeatingDP-1; i++ )
		if( t < pfBeamTime[i+1] ) break;
	x[1] = pfBeamTime[i];
	x[2] = pfBeamTime[i+1];

	// Energy flux
	y[1] = pfBeamEnergyFlux[i];
	y[2] = pfBeamEnergyFlux[i+1];
	LinearFit( x, y, t, &(pBeamParams[0]) );

	// Cut-off
	y[1] = pfBeamCutOff[i];
	y[2] = pfBeamCutOff[i+1];
	LinearFit( x, y, t, &(pBeamParams[1]) );

	// Spectral index
	y[1] = pfBeamSpectralIndex[i];
	y[2] = pfBeamSpectralIndex[i+1];
	LinearFit( x, y, t, &(pBeamParams[2]) );
}

return;
}

double CHeat::CalculateBeamHeating( double t, double *pBeamParams, double nds, double Nstar, double n_e, double n_H, double x )
{
double energy_flux, cutoff_energy, delta;
double mu0 = 1.0;	// The cosine of the pitch angle is hard-wired for the moment
double fBeamHeating;
double term1, term2;

	// Trap the special cases
	if( ( !iBeamHeatingDP ) || ( iBeamHeatingDP == 1 && t > pfBeamTime[0] ) || ( iBeamHeatingDP > 1 && ( t < pfBeamTime[0] || t > pfBeamTime[iBeamHeatingDP-1] ) ) ) return 0.0;

	// The total energy flux in the beam at the injection site (erg cm^-2 s^-1)	
	energy_flux = pBeamParams[0];

	// Cut-off energy of the beam (erg)
	// 1.602e-9 erg / keV
	cutoff_energy = (1.602e-9) * pBeamParams[1];

	// The power-law index of the beam's energy distribution
	// delta must strictly be greater than 2
	delta = pBeamParams[2];

	//Average electron energy -> Derived from weighted average
	double avg_energy = ( ( 1.0 - delta ) / ( 2.0 - delta ) ) * cutoff_energy;

	// The three Coulomb logarithms used:
	double Lambda1 = 66.0 + ( 1.5 * log(avg_energy) ) - ( 0.5 * log(n_e) );
	double Lambda2 = 25.1 + log(avg_energy);

	// The special values of gamma and beta as defined in H&F 1994 or Emslie 1978 naming them g and b to avoid confusion with the mathematical Gamma and Beta functions
	double g = ( x * Lambda1 ) + ( ( 1.0 - x ) * Lambda2 );
	// double b = 2.0; // As prescribed by H&F
	
	// Nc and Ncstar, as defined in H&F 1994
	// double Nc = ( mu0 * cutoff_energy * cutoff_energy ) / ( g * (2.0 + (b/2.0)) * 2.0 * _PI_ * pow( ELECTRON_CHARGE, 4.0 ) );
	// Simplified to:
	double Nc = ( mu0 * cutoff_energy * cutoff_energy ) / ( g * (1.0006128424679109518248146922857e-36) );

	// double Ncstar = ( mu0 * cutoff_energy * cutoff_energy ) / ( Lambda1 * (2.0 + (b/2.0)) * 2.0 * _PI_ * pow( ELECTRON_CHARGE, 4.0 ) );
	// Simplified to:	
	double Ncstar = Nc * ( g / Lambda1 );

	// The lower bound energy of the heating integral
	// double integral_lowenergy = sqrt( (2.0+b/2.0) * g * 2.0 * _PI_ * pow(ELECTRON_CHARGE,4.0) * nds / mu0 );
	// Simplified to:
	double integral_lowenergy = sqrt( ( g * nds * (1.0006128424679109518248146922857e-36) ) / mu0 );
	if( isnan( integral_lowenergy ) ) return 0.0;
	
	term1 = (3.335376141559703172749382307619e-37) * n_H * g * ( delta - 2.0 ) * energy_flux * pow( (Nstar/Ncstar), (-delta/2.0) ) * Beta( (delta/2.0), 0.33333333 );
	term2 = 2.0 * cutoff_energy * cutoff_energy;

	if ( cutoff_energy > integral_lowenergy )
		term1 *= incompleteBeta( (delta/2.0), 0.33333333, (nds/Nc) );

	fBeamHeating = term1 / term2;
	if( isnan(fBeamHeating) || fBeamHeating < 0.0 ) fBeamHeating = 0.0;

	return( fBeamHeating );
}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void CHeat::InitQbeam( double fAvgEE, int iNumCells )
{
	// Set the average electron energy
	fAverageElectronEnergy = fAvgEE;

	// If memory has already been allocated to a beam at a previous time-step then free that memory ready for the beam at the current time-step
	if( pfQbeam )
		delete[] pfQbeam;

	// Reset the index to the beam quantities
	iQbeamIndex = 0;
	 
	 iQbeamIndex_max = iNumCells * 2;
	pfQbeam = new double[iQbeamIndex_max];
}

void CHeat::SetQbeam( double s, double Qbeam )
{
	pfQbeam[iQbeamIndex] = s;
	pfQbeam[iQbeamIndex+1] = Qbeam;

	iQbeamIndex += 2;
	if( iQbeamIndex == iQbeamIndex_max )
		iQbeamIndex = 0;
}

double CHeat::GetQbeam( double s )
{
	if( !pfQbeam ) return 0.0;

	double x[3], y[3], Qbeam;
	int i;

#ifdef OPENMP
	for( i=0; i<=iQbeamIndex_max-2; i+=2 )
#else // OPENMP
	for( i=iQbeamIndex; i<=iQbeamIndex_max-2; i+=2 )
#endif // OPENMP
		if( s <= pfQbeam[i] ) break;

	if( !i ) i+=2;
	else if( i > iQbeamIndex_max-2) i = iQbeamIndex_max-2;

	x[1] = pfQbeam[i-2];
	x[2] = pfQbeam[i];
	y[1] = pfQbeam[i-1];
	y[2] = pfQbeam[i+1];
	LinearFit( x, y, s, &Qbeam );
	if( Qbeam < 0.0 ) Qbeam = 0.0;

	iQbeamIndex = i;

	return Qbeam;
}

double CHeat::GetAvgEE( void )
{
	return fAverageElectronEnergy;	
}

void CHeat::WriteQbeam( void )
{
	FILE *pFile;
	int i;
	
	pFile = fopen( "Heating_Model/config/Qbeam.dat", "w" );
		fprintf( pFile, "%.16e\n", fAverageElectronEnergy );
		fprintf( pFile, "%i\n", iQbeamIndex_max>>1 );
		for( i=0; i<iQbeamIndex_max; i+=2 )
			fprintf( pFile, "%.16e\t%.16e\n", pfQbeam[i], pfQbeam[i+1] );
	fclose( pFile );
}

void CHeat::ReadQbeam( void )
{
	FILE *pFile;
	double fAvgEE, fs, fQbeam;
	int iNumCells, i;
	
	pFile = fopen( "Heating_Model/config/Qbeam.dat", "r" );
		ReadDouble( pFile, &fAvgEE );
		fscanf( pFile, "%i", &iNumCells );
		InitQbeam( fAvgEE, iNumCells );
		for( i=0; i<iNumCells; i++ )
		{
			ReadDouble( pFile, &fs );
			ReadDouble( pFile, &fQbeam );
			SetQbeam( fs, fQbeam );
		}
	fclose( pFile );
}
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
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
#endif // OPTICALLY_THICK_RADIATION