// ****
// *
// * Forward Model class definition
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/03/2015
// *
// ****


#include "instrument.h"


class CForward {

    private:

    // The number of virtual instruments
    int iNumInstruments;

    // Array of pointers to virtual instruments
    PPINSTRUMENT ppInstrument;

    // Function to initialise the Forward Model object with a set of instruments
    void Initialise( char *szFilename );

    // Function to free all allocated memory
    void FreeAll( void );

    public:

    // Constructor
    CForward( char *szFilename );
        
    // Destructor
    ~CForward( void );

    // Function to return the number of virtual instruments
    int GetNumInstruments( void );

    // Function to write a file containing a list of the ions and emission lines
    // within the wavelength response region of the instruments
    void WriteIonLineListToFile( char *szFilename );

    // Function to create a virtual detector
    void CreateVirtualDetector( PLOOP pLoop );

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

typedef CForward* PFORWARD;
