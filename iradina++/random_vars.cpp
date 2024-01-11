#include "random_vars.h"

#include <algorithm>

/* prime numbers for lists length */
#define MAXLOGLIST 104729
#define MAXAZILIST 72211   //  (1 << 17)
#define MAXRANLIST 1000003
#define MAXERFLIST 79999

/*
 * Single precision implementation of erf^{-1}(x)
 *
 * From https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
 * Prof. Mike Giles
 *
 * Inverse cum. PDF of Normal = \sqrt{2} \erfinv( 2*u - 1 )
 * where u = [0, 1)
 *
 */
float MBG_erfinv(float x)
{
    float w, p;
    w = - std::log((1.0f-x)*(1.0f+x));
    if ( w < 5.000000f ) {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p*w;
        p = -3.5233877e-06f + p*w;
        p = -4.39150654e-06f + p*w;
        p = 0.00021858087f + p*w;
        p = -0.00125372503f + p*w;
        p = -0.00417768164f + p*w;
        p = 0.246640727f + p*w;
        p = 1.50140941f + p*w;
    }
    else {
        w = sqrtf(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p*w;
        p = 0.00134934322f + p*w;
        p = -0.00367342844f + p*w;
        p = 0.00573950773f + p*w;
        p = -0.0076224613f + p*w;
        p = 0.00943887047f + p*w;
        p = 1.00167406f + p*w;
        p = 2.83297682f + p*w;
    }
    return p*x;
}

/* generate and randomize lists randomlist,sqrtloglist, sinAzimAngle & cosAzimAngle   */
template<class _E>
void random_vars_tbl<_E>::initTables()
{

    iUniform = iSqrtLog = iAzimuth = iInvErf = 0;

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

    /*
     * Generate a list of evenly distributed but randomly ordered
     * values of inverse erf of -1+1/MAXERFLIST to 1-1/MAXERFLIST
     *
     * Multiply by sqrt(2) to make suitable for Normal distr. sampling
     *
     */
    inverflist.resize(MAXERFLIST);
    const float sqrt2 = std::sqrt(2.f);
    for(int i=0; i<MAXERFLIST; i++) {
        float x = 2.f*(i+1)/(MAXERFLIST+1) - 1;
        inverflist[i] = sqrt2*MBG_erfinv(x);
    }
    std::shuffle(inverflist.begin(), inverflist.end(), urbg);

}

template void random_vars_tbl<std::mt19937>::initTables();
template void random_vars_tbl<std::minstd_rand>::initTables();
