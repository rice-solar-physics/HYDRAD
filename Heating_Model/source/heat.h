// ****
// *
// * A heating code to simulate different forms of heat deposition
// *
// * Class function definitions
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/06/2017
// *
// ****

#include "config.h"

class CHeat {

    private:

    // Quiescent heating parameters
	// Location, scale-length and maximum energy
    	double s0quiescent, sHquiescent, E0quiescent;
	
    // Episodic heating parameters
	// The number of episodic heating events to be activated
    	int NumActivatedEvents;

        // Location, scale-length and maximum energy
        double *s0episodic, *sHepisodic, *E0episodic;

        // Start and end times of rise phase
        double *tsRepisodic, *teRepisodic;
        // Start and end times of decay phase
        double *tsDepisodic, *teDepisodic;

#ifdef BEAM_HEATING
    // Beam heating parameters
        int iBeamHeatingDP;
        double *pfBeamTime, *pfBeamEnergyFlux, *pfBeamCutOff, *pfBeamSpectralIndex;
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Background volumetric heating required to maintain the VAL atmosphere in equilibrium
        int iVALHeatingDP;
        double **ppVALHeating;
#endif // OPTICALLY_THICK_RADIATION

    void Initialise( void );
    void FreeAll( void );

    // Time independent / dependent (thermal) heating functionality
    void GetHeatingData ( void );

#ifdef BEAM_HEATING
    // Beam heating functionality
    void GetBeamHeatingData( void );
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Heating to maintain the lower atmosphere
    void GetVALHeatingData( void );
#endif // OPTICALLY_THICK_RADIATION

    public:

    CHeat( void );
    ~CHeat( void );

    // Time independent / dependent (thermal) heating functionality
    double CalculateQuiescentHeating( double s );
    double CalculateEpisodicHeating( double s, double t );
    double CalculateHeating( double s, double t );

#ifdef BEAM_HEATING
    // Beam heating functionality
    void CalculateBeamParameters( double t, double *pBeamParams );
    double CalculateBeamHeating( double t, double *pBeamParams, double nds, double Nstar, double n_e, double n_H, double x );
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Heating to maintain the lower atmosphere
    double CalculateVALHeating( double flog10_rho_c );
#endif // OPTICALLY_THICK_RADIATION

};

typedef CHeat* PHEAT;
