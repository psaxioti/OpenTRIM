#include "corteo.h"

#include <iostream>
#include <fstream>

#include "xs.h"

using std::cout;
using std::endl;
using std::ofstream;

int gencorteo4bit()
{
    unsigned long nThetaErr = 0;
    reducedXS_quad< screening<reducedXS::ZBL> > xs_quad_zbl;
    ofstream ofs("corteo4bit.cpp");
    int k = 0;
    int klast = corteo4bit::rows * corteo4bit::cols - 1;

    cout << "Computing 4-bit corteo scattering table";
    // compute matrix for each reduced energy, reduced impact parameter pair

    for(corteo4bit::e_index ie; ie!=ie.end(); ie++) {

        if (ie % (corteo4bit::rows/10)==0) {
            cout << ".";
            cout.flush();
        }

        for(corteo4bit::s_index is; is!=is.end(); is++) {
            // calculations (summation) made using double to decrease numerical noise
            float sin2ThetaBy2 = xs_quad_zbl.sin2Thetaby2(*ie, *is);
            if(std::isnan(sin2ThetaBy2)) nThetaErr++;

            // store in matrix sin^2(theta_CM/2) as float
            // matrix[is+ie*s_index::dim] = (float)(sin2ThetaBy2);
            ofs << sin2ThetaBy2;
            if (k!=klast) ofs << ',';
            k++;
            if (k % 8 ==0) ofs << endl;
        }
    }

    cout << " done.\n";
    cout.flush();

    if(nThetaErr)
        cout << "ERROR: " << nThetaErr << " error(s) in theta evaluation.\n";

    return !nThetaErr;
}
