// ****
// *
// * #defines for configuring the hydrostatic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by HYDRAD_GUI on 16-11-2017 16:04:50
// *
// ****

// **** Output ****
// **** End of Output ****

// **** Physics ****
// Radiation //
#include "../../Radiation_Model/source/config.h"
// End of Radiation //
// Gravity //
// End of Gravity //
// **** End of Physics ****

// **** Solver ****
#define EPSILON 0.01
// **** End of Solver ****

// **** Grid ****
#define ADAPT
#define MIN_CELLS 150
#define MAX_CELLS 30000
#define MAX_REFINEMENT_LEVEL 12
#define MIN_DS 1e0
#define MAX_VARIATION 1.10
// **** End of Grid ****
