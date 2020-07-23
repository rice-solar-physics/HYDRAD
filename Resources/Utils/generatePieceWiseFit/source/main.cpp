
#include <stdio.h>

#include "piecewisefit.h"


int main( void )
{
	PPIECEWISEFIT pMagneticField;

	pMagneticField = new CPieceWiseFit( (char*)"fits/B(closed).txt", (char*)"fits/initial.amr.B" );
	// pMagneticField = new CPieceWiseFit( (char*)"fits/B(open).txt", (char*)"fits/initial.amr.B" );
	// pMagneticField = new CPieceWiseFit( (char*)"fits/initial.amr.B" );
#ifdef VERBOSE
	pMagneticField->ShowPieceWiseFit();
#endif // VERBOSE
	pMagneticField->GetPieceWiseFit( 0.25 );
	pMagneticField->GetPieceWiseFit( 0.5 );
	pMagneticField->GetPieceWiseFit( 0.75 );

	pMagneticField->GetDerivative( 0.25, 1 );
	pMagneticField->GetDerivative( 0.5, 1 );
	pMagneticField->GetDerivative( 0.75, 1 );

	delete pMagneticField;

	return 0;
}
