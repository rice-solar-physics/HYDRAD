// ****
// *
// * Definition of the user parameters structure
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 05/11/2012
// *
// ****


// Define the user parameters structure
struct Parameters {

    // Output path and filename
    char szOutputFilename[256];

    // Loop geometry
    double Lfull, s0, Inc;

    // Lower boundary conditions
    double T0, n0;

    // Solution parameter space
    double sH0, sH, Log_10H0[2], dLog_10H0, Hintervals;

};

// Define a type for the user parameters structure
typedef Parameters PARAMETERS;
