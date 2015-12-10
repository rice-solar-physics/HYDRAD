// ****
// *
// * Function bodies for the class definition of
// * the adaptive mesh cells
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 05/11/2012
// *
// ****


#include <string.h>

#include "cell.h"


// **** ADAPTIVE MESH CELL CLASS ****

// Constructor
CAdaptiveMeshCell::CAdaptiveMeshCell( PCELLPROPERTIES pInitCellProperties )
{
UpdateCellProperties( pInitCellProperties );
}

// Destructor
CAdaptiveMeshCell::~CAdaptiveMeshCell( void )
{
}

void CAdaptiveMeshCell::SetPointer( int iPointer, CAdaptiveMeshCell *pPointer )
{
pNearestCell[iPointer] = pPointer;
}

CAdaptiveMeshCell *CAdaptiveMeshCell::pGetPointer( int iPointer )
{
return pNearestCell[iPointer];
}

void CAdaptiveMeshCell::UpdateCellProperties( PCELLPROPERTIES pCellProperties )
{
memcpy( this, pCellProperties, sizeof(CELLPROPERTIES) );
}

void CAdaptiveMeshCell::GetCellProperties( PCELLPROPERTIES pCellProperties )
{
memcpy( pCellProperties, this, sizeof(CELLPROPERTIES) );
}
