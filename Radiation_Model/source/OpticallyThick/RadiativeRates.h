// ****
// *
// * Class definition of the radiative transition rates
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 06/06/2018
// *
// ****

#include "../config.h"

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE

#define MAX_ITERATIONS			300
// #define CONVERGENCE_EPSILON		0.10		// Original
#define CONVERGENCE_EPSILON		0.25
#define CONVERGENCE_CONDITION		1E-6


// **** RADIATIVE TRANSITION RATES CLASS ****

// Define the radiative transition rates class
class CRadiativeRates {
	
    private:

    // Radiative bound-bound rates:
	int iBB_TRvals, iBB_Tvals;
	double *pfBB_logT;
	// Photoexcitation:
	double **ppfBB_lu;
	// Radiative decay:
	double **ppfBB_ul;

    // Radiative bound-free rates:
	int iBF_TRvals, iBF_Tvals;
	double *pfBF_logT;
	// Photoionization
	double **ppfBF;

    // Radiative free-bound rates:
	// The radiation temperatures are in columns and the electron temperatures are in rows
	int iFB_TRvals, iFB_radTvals, iFB_eTvals;
	double *pfFB_logradT, *pfFB_logeT;
	// Radiative de-excitation
	double ***pppfFB;
    
    // Collisional bound-bound and bound-free rates
	int iColl_Tvals;
	double *pfColl_logT;
	// Collisional excitation
	int iColl_exTRvals;
	double **ppfColl_ex_lu, **ppfColl_ex_ul;
	// Collisional ionization
	int iColl_ionTRvals;
	double **ppfColl_ion, **ppfColl_rec;
    
    int iNBBT, iNBFT;
    double *pTrt, *pnu0, *pTeZ_c, *pZ_c_LEFT, *pZ_c_RIGHT, *pMcZ_c_LEFT, *pMcZ_c_RIGHT;
    double *pterm1, *pterm2, *pH;

    // Functions to calculate the singular value decomposition of a matrix A and back-substitute to solve a matrix equation A x = b
    // (Used by SolveHIIFraction)
        // Singular value decomposition of matrix A of size m x n
	int svdcmp(double **a, int m, int n, double *w, double **v);
	// Back substitution of the singular value decomposition to solve for vector x
	void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);
	// Pythagorean distance
	double pythag(double a, double b);

    void GetBBRates( char *pszBBRatesFile );
    void GetBFRates( char *pszBFRatesFile );
    void GetFBRates( char *pszFBRatesFile );
    void GetCollRates( char *pszCollRatesFile );

    void Initialise( char *pszRatesFiles );
    void FreeAll( void );

    public:

    // Constructor
    CRadiativeRates( char *pszRatesFiles );
	
    // Destructor
    ~CRadiativeRates( void );

    // Radiative transition rates class functions

    void GetBoundBoundRates( double *pfBB_lu, double *pfBB_ul, double *pflog10T );
    void GetBoundFreeRates( double *pfBF, double *pflog10T );
    void GetFreeBoundRates( double *pfFB, double *pflog10radT, double flog10eT, double fne );
    void GetFreeBoundRatesRH( double *pfFB, double *pflog10radT, double flog10eT, double fne, double *pfLevel_Ratio );
    void GetCollisionalRates( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne );
    void GetCollisionalRatesRH( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne );
    
    void SolveHIIFraction( double *pfHstate, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB );
        
    void GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double delta_s, double *pDel_Hstate_dot_v );
    void Normalize( double *pHstate );

    int GetNumberOfTransitions( void );
    double GetZ_c_LEFT( int i );
    double GetZ_c_RIGHT( int i );
    double GetMcZ_c_LEFT( int i );
    double GetMcZ_c_RIGHT( int i );
    double GetTeZ_c( int i );
    double GetTerm1( int i );
    double GetTerm2( int i );
    double GetNu0( int i );
    double GetH( int i );
    void SetZ_c_LEFT( int i, double Z_c );
    void SetZ_c_RIGHT( int i, double Z_c );
    void SetMcZ_c_LEFT( int i, double McZ_c );
    void SetMcZ_c_RIGHT( int i, double McZ_c );
    
};

// Define a type for the radiative transition rates class
typedef CRadiativeRates* PRADIATIVERATES;

#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
