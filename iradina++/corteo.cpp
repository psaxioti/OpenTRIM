#include "corteo.h"
#include "xs.h"

#include <stdexcept>

corteo::corteo() :
    matrix(epsilon_index::dim * s_index::dim)  // (DIME*DIMS)
{
}

/************** Adapted from corteo20130715 *********************/
/* Length of float and int must be 4 byte, otherwise indexing doesn't
work correctly. */
int corteo::check_type_representation()
{
    // compile-time checks
    static_assert(sizeof(unsigned int)==4, "size of int not 4");
    static_assert(sizeof(float)==4, "size of float not 4");

    // run-time check of correct float/int conversion / IEEE-754 conformance
    float f = 3328.625f;
	if(*(unsigned int *)&f != 0x45500a00)  {
        std::runtime_error(
            "This machine does not follow IEEE-754 standard for binary representation of float.");
        return -1;
	}
	return 0;
}

/*************** Functions adapted from corteo.c ***********************/

int corteo::calcMatrix(const xs_base& xs, std::ostream *os)
{
    unsigned long nThetaErr = 0;

    if (os) *os << "Computing scattering matrix ";
    // compute matrix for each reduced energy, reduced impact parameter pair

    for(epsilon_index ie; ie!=ie.end(); ie++) {

        if(os && (int(ie) % (epsilon_index::dim/10)==0) ) {
            *os << ".";
            os->flush();
        }

        for(s_index is; is!=is.end(); is++) {
            // calculations (summation) made using double to decrease numerical noise
            double sin2ThetaBy2 = xs.sin2Thetaby2(*ie, *is);
            if(std::isnan(sin2ThetaBy2)) nThetaErr++;

            // store in matrix sin^2(theta_CM/2) as float
            matrix[int(is)+int(ie)*s_index::dim] = (float)(sin2ThetaBy2);
        }
    }

    if(os) {
        *os << " done.\n";
        os->flush();
    }

    if(nThetaErr && os)
        *os << "ERROR: " << nThetaErr << " error(s) in theta evaluation.\n";

    return !nThetaErr;
}



