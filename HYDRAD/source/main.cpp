// ****
// *
// * HYDrodynamics and RADiation code (HYDRAD)
// * An astrophysical fluid dynamics code for 
// * solar and stellar atmospheres
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/02/2007
// *
// ****


#include "mesh.h"


int main( void )
{
PMESH pMesh;

// Set up the problem, create the mesh and solve the equations
pMesh = new CAdaptiveMesh;

// Delete the mesh
delete pMesh;

return 0;
}
