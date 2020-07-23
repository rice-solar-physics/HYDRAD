// ****
// *
// * Defines for particular properties of the simulation and definition of the parameters structure
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****


#include "config.h"

#include "../../Resources/Utils/generatePieceWiseFit/source/piecewisefit.h"

#include "../../Heating_Model/source/heat.h"

#ifdef NON_EQUILIBRIUM_RADIATION
#include "../../Radiation_Model/source/ionfrac.h"
#else // NON_EQUILIBRIUM_RADIATION
#include "../../Radiation_Model/source/radiation.h"
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#include "../../Radiation_Model/source/OpticallyThick/OpticallyThickIon.h"
#ifdef NLTE_CHROMOSPHERE
#include "../../Radiation_Model/source/OpticallyThick/RadiativeRates.h"
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#ifdef USE_KINETIC_MODEL
#include "../../Kinetic_Model/source/kinetic.h"
#endif // USE_KINETIC_MODEL

// **** PARAMETERS STRUCTURE ****

// Define the parameters structure
struct Parameters {

    // Mass density, momentum density and energy density profiles
    char Profiles[256];

    // Loop length
    double L;

	// Number of grid cells
	int iNumberOfCells;

    // Gravitational acceleration polynomial-fit coefficients filename
    char GravityFilename[256];

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
	// Magnetic field polynomial-fit coefficients filename
    char MagneticFieldFilename[256];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

    // Duration and profile output period
    double Duration, OutputPeriod;

};

// Define a type for the parameters structure
typedef Parameters PARAMETERS;
