// ****
// *
// * Code to calculate the electron distribution function, using the Spitzer-Harm method for
// * the bulk of the distribution function and the BGK approximation for the tail,
// * given a spatial temperature and density distribution in one dimension
// *
// * Class function definitions
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 05/11/2012
// *
// ****


// Definitions for the properties of the distribution function

#define V_RANGE_UPPER			12.0	// Upper velocity bound relative to the greatest thermal speed
#define DISTRIBUTION_DATA_POINTS	1000	// Number of discrete data points

#define F_MIN				1e-300	// The smallest value of the distribution function (to avoid underflow errors)

#define KNUDSEN_NUMBER			2e-3	// Ratio of the electron mean-free-path to the temperature scale-length at
						// which the Spitzer-Harm heat flux becomes invalid

// #define USE_APPROXIMATE_NU_EI		// Use to speed up the calculation of nu_ei (else uses more detailed equation in NRL Plasma Formulary)


double fLogLambda_ee( double Te, double n );
double fLogLambda_ei( double Te, double Ti, double n );
double fSHMeanFreePath( double u_th, double Te, double Ti, double n );


class CKinetic {

	private:

	// The velocities, collision frequencies, Maxwellian and non-Maxwellian distribution functions
	double u_th, *pupsilon, *pnu_ee, *pnu_ei, *plambda_ei, *pMaxDFN_ee, *pMaxDFN_ei, *pNonMaxDFN;

	void Initialise( void );
	void FreeAll( void );

	double CalculateThermalSpeed( double T );

	public:

	CKinetic( void );
	~CKinetic( void );

	void CalculateVelocityRange( double T, double Tmax );
	void CalculateCollisionalProperties( double Te, double Ti, double n );
	void CalculateMaxDFN( double Te, double Ti, double n );

	double Get_u_th( void );
	
	double* Get_pupsilon( void );
	double* Get_pnu_ee( void );
	double* Get_pnu_ei( void );
	double* Get_plambda_ei( void );
	double* Get_pMaxDFN_ee( void );
	double* Get_pMaxDFN_ei( void );
	double* Get_pNonMaxDFN( void );

	double Get_rho( void );
	double Get_TE_KE( void );
	double Get_Fc( void );

};

typedef CKinetic* PKINETIC;
