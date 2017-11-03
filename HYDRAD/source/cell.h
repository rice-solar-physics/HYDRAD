// ****
// *
// * Class definition of the adaptive mesh cells
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/01/2017
// *
// ****


#include "params.h"


#define RIGHT			1		// The array index to reference the grid cell to the right of the current cell
#define LEFT			2		// The array index to reference the grid cell to the left of the current cell
#define BOTTOM			3		// The array index to reference the grid cell below (i.e. one time-step behind) the current cell

#define FALSE			0		// Logical FALSE
#define TRUE			1		// Logical TRUE

#define SPECIES			2		// The number of particle species comprising the plasma

#define ELECTRON		0		// The array index to reference the electron species
#define HYDROGEN		1		// The array index to reference the hydrogen species

#define MASS_TERMS		1		// The number of terms in the mass conservation equation
#define MOMENTUM_TERMS          5		// The number of terms in the momentum conservation equation
#define ENERGY_TERMS            9		// The number of terms in the energy conservation equations


// **** ADAPTIVE MESH CELL PROPERTIES ****

// Define the adaptive mesh cell properties structure
struct AdaptiveMeshCellProperties {

    // The prolongation condition and refinement level of the current cell, and
    // the ID number of the only cell it can be merged with at the prolongation step
    int iRefinementLevel, iUniqueID[MAX_REFINEMENT_LEVEL+1];

    // Information about the location of the cell on the grid
    double s[3], cell_width;

    // Conserved quantities
    double rho[3], rho_v[3], TE_KE[3][SPECIES];

    // Physical quantities
    double n[SPECIES], v[3], T[SPECIES], P[3][SPECIES], TE_KE_P[3][SPECIES], nu_ie, Cs, M;

    // Thermal flux
    double Fc[3][SPECIES];

    // Viscosity
    double eta, Feta[3];

    // Numerical viscosity
    double Fnumerical[3];

#ifdef NON_EQUILIBRIUM_RADIATION
    // Fractional population of the ions in the current cell
    PIONFRAC pIonFrac;
	
    // The smallest ionisation / recombination time-scale in the current cell
    double atomic_delta_t;
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
    // The neutral hydrogen population fraction, and the column number and mass densities
    double HI, HI_c, rho_c;
	#ifdef BEAM_HEATING
	    // The hydrogen and modified hydrogen column densities
	    double nH_c, nH_star_c;
	#endif // BEAM_HEATING
	#ifdef NLTE_CHROMOSPHERE
	    double rho_e, Mc, Trad[10], Hstate[6], Del_Hstate_dot_v[6];
	#endif // NLTE_CHROMOSPHERE
#else // // OPTICALLY_THICK_RADIATION
	#ifdef BEAM_HEATING
	    // The neutral hydrogen population fraction, and the hydrogen and modified hydrogen column densities
	    double HI, nH_c, nH_star_c;
	#endif // BEAM_HEATING
#endif // OPTICALLY_THICK_RADIATION

#ifdef USE_KINETIC_MODEL
    // Kinetic model: contains Maxwellian and non-Maxwellian distributions; and functions to calculate moments
    PKINETIC pKinetic;
#endif // USE_KINETIC_MODEL

    // Terms of the equations
	
    // Mass conservation equation
    double rho_term[MASS_TERMS], drhobydt;

    // Momentum equation
    double rho_v_term[MOMENTUM_TERMS], drho_vbydt;

    // Energy equation
    double TE_KE_term[ENERGY_TERMS][SPECIES], dTE_KEbydt[SPECIES];

    // Characteristic time-scales in the current cell
    double advection_delta_t, conduction_delta_t[SPECIES], collision_delta_t, radiation_delta_t, viscosity_delta_t;

};

// Define types for the adaptive mesh cell properties structure
typedef AdaptiveMeshCellProperties CELLPROPERTIES;
typedef AdaptiveMeshCellProperties* PCELLPROPERTIES;

// **** ADAPTIVE MESH CELL CLASS ****

// Define the adaptive mesh cell class
class CAdaptiveMeshCell : private AdaptiveMeshCellProperties {
	
    private:
	
    // Pointers to the adjacent cells
    CAdaptiveMeshCell *pNearestCell[4];

    public:
	
    // Constructor
    CAdaptiveMeshCell( PCELLPROPERTIES pInitCellProperties );
	
    // Destructor
    ~CAdaptiveMeshCell( void );

    // Adaptive mesh cell functions

    void UpdateCellProperties( PCELLPROPERTIES pCellProperties );
    void GetCellProperties( PCELLPROPERTIES pCellProperties );

    void SetPointer( int iPointer, CAdaptiveMeshCell *pPointer );
    CAdaptiveMeshCell *pGetPointer( int iPointer );

};

// Define a type for the adaptive mesh cell class
typedef CAdaptiveMeshCell* PCELL;
