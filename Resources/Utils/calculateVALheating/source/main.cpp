//
//  updateVAL.cpp
//
//  Updates the VAL heating profile in the case that the initial
//  density profile has been altered.
//
//  Created by Jeffrey Reep on 7/27/15.
//
//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../source/constants.h"
#include "../../../../Radiation_Model/source/OpticallyThick/OpticallyThickIon.h"

using namespace std;

int main()
{
    
    /*ifstream fin("Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat");
    int heat_num_lines,i;
    double **heat_rows;
    fin >> heat_num_lines;
    heat_rows = (double**)malloc( sizeof(double) * 2 );
    for( i=0; i<2; ++i )
        heat_rows[i] = (double*)malloc( sizeof(double) * heat_num_lines );
    for( i=0; i<heat_num_lines; ++i )
    {
        // Array index [0][i] contain the spatial coordinate and [1][i] contain the volumetric heating rate
        fin >> heat_rows[0][i] >> heat_rows[1][i];
    }
    
    fin.close();
    */
    
    
    // Read in the loop length
    char dummy[256];
    double fLength;
    ifstream fin2("Initial_Conditions/config/initial_conditions.cfg");
    fin2 >> dummy;
    fin2 >> fLength;
    fin2.close();
    
    unsigned int number_of_lines = 0;
    FILE *infile = fopen("Initial_Conditions/profiles/initial.amr.phy", "r");
    int ch;
    
    while (EOF != (ch=getc(infile)))
        if ('\n' == ch)
            ++number_of_lines;
    
    fclose( infile );
    
    ifstream fin3("Initial_Conditions/profiles/initial.amr.phy");

    double frho_c = 0.;
    double fHI_c =0.;
    double previous_s = fLength/2.;
    double cell_width_cos_theta = 0.;
    double fDensityDifference,fRadiation=0.;
    double s,n_e,n_H,T,P_e,P_H;
    
    COpticallyThickIon *pHI, *pMgII,*pCaII;
    
    
    
    
    pHI = new COpticallyThickIon( 1, (char *)"h_1", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    pMgII = new COpticallyThickIon( 12, (char *)"mg_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    pCaII = new COpticallyThickIon( 20, (char *)"ca_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
    
    // Create vectors to store the log of radiation and column density
    vector<double> fRad;
    vector<double> fRho;
    
    // Reserve space for the maximum size of the arrays
    fRad.reserve(number_of_lines/2);
    fRho.reserve(number_of_lines/2);
    int num_elements = 0;
    int i;
    
    for( i=0;i<number_of_lines;++i )
    {
        // Read in the line and values on the line
        fin3 >> s >> n_e >> n_H >> T >> P_e >> P_H;
        
        // For positions in the second half of the loop, calculate the column densities
        if( i > number_of_lines/2 )
        {
            
            if( T < OPTICALLY_THICK_TEMPERATURE )
            {
                cell_width_cos_theta = fabs(s - previous_s) * fabs( cos( ( _PI_ * s ) / fLength ) );
                
                frho_c += n_H * AVERAGE_PARTICLE_MASS * cell_width_cos_theta;
                
                // Calculate the neutral hydrogen column number density
                fDensityDifference = n_H - ( n_e / 1.000144 );
                if ( fabs(fDensityDifference)  < 1.0 )
                    fHI_c += 0.0;
                else
                    fHI_c += ( n_H - ( n_e / 1.000144 ) ) * cell_width_cos_theta;
                
                fRho.push_back(log10(frho_c));
                
                fRadiation = pHI->GetVolumetricLossRate( log10(T), log10((4e-14)*fHI_c), n_e * n_H * AVERAGE_PARTICLE_MASS) + pMgII->GetVolumetricLossRate( log10(T), log10(frho_c), n_e * n_H * AVERAGE_PARTICLE_MASS) + pCaII->GetVolumetricLossRate( log10(T), log10(frho_c), n_e * n_H * AVERAGE_PARTICLE_MASS);
             
                fRad.push_back(log10(fRadiation));
                
                //if( log10(frho_c) < 0. )
                //    printf("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",s,log10(frho_c),fRadiation,fHI_c,cell_width_cos_theta,fabs( cos( ( _PI_ * s ) / fLength ) ), fDensityDifference);
                
                ++num_elements;
                
            }
        }
        
        previous_s = s;
    }
    
    fclose(infile);
    
    // Output the column density vs heating into the heat profile
    ofstream outfile("Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat");
    outfile << num_elements << endl;
    outfile.precision(10);
    //outfile << -10.0 << "\t" << -1.901435 << endl;
    for( i=0; i<num_elements; ++i)
    {
        if( fRho[i] >= 0.0 ) fRad[i] = -300.0;
        outfile << fRho[i] << "\t" << fRad[i] << endl;
    }
    outfile.close();
    
    return 0;
}
