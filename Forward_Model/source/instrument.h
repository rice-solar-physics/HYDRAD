// ****
// *
// * Instrument Class definition
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


#include "ion.h"
#include "loop.h"
#include "../../Radiation_Model/source/radiation.h"
#include "../../Radiation_Model/source/config.h"


#define CM_PER_ARCSEC	7.6E7


class CInstrument {

    private:
        
    // The name of the instrument
    char szName[256];

    // The instrument response function. The wavelengths and coefficients are paired in the array:
    // pRespFunc[0] = Lambda1, pRespFunc[1] = IRF1, pRespFunc[2] = Lambda2, pRespFunc[3] = IRF2, etc.
    int iNumRespFuncDataPoints;
    double *pRespFunc;

    // The wavelength sensitivity range of the instrument obtained from the response function
    double fLambda_min, fLambda_max;

    // The spectral resolution and instrumental line width for spectroscopic instruments
    double fSpecRes, fInstrumentalLineWidth;

    // The instrument pixel size in arcsec and cm ([0] = solar-X, [1] = solar-Y)
    double fPixel_arcsec[2], fPixel_cm[2];
    double fPixel_Area_arcsec, fPixel_Area_cm;

    // The instrument temporal cadence in seconds
    double fCadence;

    // The virtual detectors (should eventually be 2D)
    int iXY[2];                                             // Number of detector pixels in the solar-X and Y directions
    int iNumSpectralDataPoints;                             // Number of data points used to plot the spectrum for each pixel
    double fR;                                              // Loop radius
    double **ppfSolarXY;                                    // Solar-X and Y coordinates in arcsec at centre of each pixel
    double *pfEQDetector, *pfNEQDetector;                   // EQ and NEQ DN pixel^-1 s^-1
    double *pfLambda;                                       // Wavelength scale (A)
    double **ppfEQDetector_SPEC, **ppfNEQDetector_SPEC;     // Array of X pixels, each with a corresponding spectrum
                                                            // Note: should eventually be 2D (X x Y array of pixels)

    // The number of ions that emit lines within the wavelength range of interest
    int iNumIons;
        
    // Pointer to an array of pointers of type PION
    PPION ppIon;
        
    // Pointers to the Radiation objects
    PRADIATION pEQRadiation, pNEQRadiation;

    // Function to initialise the instrument object with a set of ions and response functions
    void Initialise( char *pszName );

    // Function to read the instrument response function part of the virtual instrument configuration file
    void GetRespFunc( FILE *pFile );

    // Function to read the ion data part of the virtual instrument configuration file
    void GetIonData( FILE *pFile );

    // Function to free all allocated memory
    void FreeAll( void );

    // Function to return the mass (g) of the ion
    double GetMass ( int iZ );

    // Function to return the number of emission lines for a particular ion
    int GetNumLines( int iZ, int iSpecNum );

    // Function to write the emission line wavelength list for a particular ion into a pre-allocated array
    void GetLineList( int iZ, int iSpecNum, double *pfLineList );

    // Function to return the emission of a particular line from a particular ion as a function of temperature and density
    double GetLineEmission( int iZ, int iSpecNum, double fLambda, double flog10T, double flog10n );

    // Function to return the total line emission from a particular ion as a function of temperature and density (summed over the wavelength range)
    double GetIonEmission( int iZ, int iSpecNum, double flog10T, double flog10n );
 
    // Functions to record the specified counts (DN pixel^-1 s^-1) at the appropriate detector pixel
    void Detect( double fs, double fds, double fSD, double fLength, double fEQCount, double fNEQCounts );
    void Detect( double fs, double fds, double fSD, double fLength, double fEQCount );

    // Functions to add the total counts for an emission line (DN pixel^-1 s^-1) to the spectrum at the appropriate detector pixel (DN pixel^-1 s^-1 A^-1)
    void Detect( PHYDATA PHYData, double fSD, double fLength, double fEQCounts, double fNEQCounts, double fRestLambda, double *fSpectralProperties );
    void Detect( PHYDATA PHYData, double fSD, double fLength, double fEQCounts, double fRestLambda, double *fSpectralProperties );

    // Functions to forward model the emission along the loop for spectral and imaging instruments
    // Units: DN pixel^-1 s^-1 (imaging); DN pixel^-1 s^-1 A^-1 (spectroscopic)
    void ForwardModel_IMAGING( PLOOP pLoop );
    void ForwardModel_SPECTROSCOPIC( PLOOP pLoop );

    public:

    // Constructor
    CInstrument( char *pszName );
        
    // Destructor
    ~CInstrument( void );

    // Function to return the name of the instrument
    void GetInstrumentName( char *pszName );

    // Function to write a file containing a list of the ions and emission lines within the wavelength response region of the instrument
    void WriteIonLineListToFile( FILE *pFile );

    // Function to create a virtual detector
    void CreateVirtualDetector( PLOOP pLoop );

    // Functions to report the presence of detectors to store equilibrium and non-equilibrium emission
    bool IsEquilibriumDetector( void );
    bool IsNonEquilibriumDetector( void );

    // Function to delete a previously created virtual detector
    void DeleteVirtualDetector( void );

    // Function to reset the virtual detector by setting measured intensities to zero
    void ResetVirtualDetector( void );

    // Functions to write the detector data to a file
    void WriteDetectorToFile( char *pszWorkingDirectory );
    void WriteDetectorToFile( char *pszWorkingDirectory, int iFileNumber );

    // Function to forward model the emission along the loop
    void ForwardModel( PLOOP pLoop );

    // Functions to generate EM loci plots
    void EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory );
    void EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory, int iFileNumber );
	
};

typedef CInstrument* PINSTRUMENT;
typedef CInstrument** PPINSTRUMENT;
