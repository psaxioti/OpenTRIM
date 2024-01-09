#include "random_vars.h"

#include <algorithm>

/* prime numbers for lists length */
#define MAXLOGLIST 104729
#define MAXAZILIST 72211   //  (1 << 17)
#define MAXRANLIST 1000003

/* generate and randomize lists randomlist,sqrtloglist, sinAzimAngle & cosAzimAngle   */
void random_vars_tbl::initTables()
{
    randomlist.resize(MAXRANLIST);
    sqrtrandomlist.resize(MAXRANLIST);
    /* generate a uniformly spaced list of values between .5/MAXRANLIST and 1-.5/MAXRANLIST */
    for(int i=0; i<randomlist.size(); i++)
        randomlist[i] = (i+0.5f)/randomlist.size();
    std::shuffle(randomlist.begin(), randomlist.end(), urbg); /* put the list in random order */

    for(int i=0; i<randomlist.size(); i++)
        sqrtrandomlist[i] = std::sqrt(randomlist[i]); /* also compute sqrt of these values */

    /* produce a list of MAXLOGLIST sqrt(-log(x)) values for x uniformly distributed */
    /* ...between x=.5/MAXLOGLIST and x=1-.5/MAXLOGLIST */
    sqrtloglist.resize(MAXLOGLIST);
    sqrtloglist1.resize(MAXLOGLIST);
    for(int i=0; i<sqrtloglist.size(); i++)
        sqrtloglist[i] = std::sqrt(-std::log((i+.5)/sqrtloglist.size()));
    std::shuffle(sqrtloglist.begin(), sqrtloglist.end(), urbg); /* put the list in random order */

    /* precompute 1/sqrtloglist[] */
    for(int i=0; i<sqrtloglist.size(); i++)
        sqrtloglist1[i] = 1.0f/sqrtloglist[i];

    /* produce a list uniformly distributed but randomly ordered azimutal angles  */
    cosAzimAngle.resize(MAXAZILIST);
    sinAzimAngle.resize(MAXAZILIST);
    for(int i=0; i<cosAzimAngle.size(); i++)
        /* cosAzimAngle temporarly contains angles */
        cosAzimAngle[i]= 2*M_PI*i/cosAzimAngle.size();
    std::shuffle(cosAzimAngle.begin(), cosAzimAngle.end(), urbg); /* put the list in random order */

    for(int i=0; i<cosAzimAngle.size(); i++) {
        /* compute the cos and sine of these angles */
        // float phi = 2*M_PI*i/cosAzimAngle.size();
        sinAzimAngle[i] = std::sin(cosAzimAngle[i]);
        cosAzimAngle[i] = std::cos(cosAzimAngle[i]);
    }
}
