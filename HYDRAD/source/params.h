// ****
// *
// * Defines for particular properties of the simulation and definition of the parameters structure
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/31/2017
// *
// ****


#include "config.h"

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

    // Initial (t=0) loop density, momentum density and energy density profiles
    char Profiles[256];

    // Loop length
    double L;

    // Gravity look-up table filename
    char GravityFilename[256];

#ifdef USE_TABULATED_CROSS_SECTION
    char CrossSectionFilename[256];
#endif // USE_TABULATED_CROSS_SECTION

    // Duration and profile output period
    double Duration, OutputPeriod;

};

// Define a type for the parameters structure
typedef Parameters PARAMETERS;
