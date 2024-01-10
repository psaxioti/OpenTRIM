#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include "urbg.h"
#include <vector>

struct random_vars_base
{

};

struct random_vars
{
    URBG& urbg;
    std::normal_distribution<float> N;
    std::uniform_real_distribution<float> X;

    random_vars(URBG& g) : urbg(g)
    {}

    /*
     * Used for impact factor sampling
     *
     * $f p \propto \sqrt{u} $f
     *
     */
    virtual void u_sqrtu(float& u, float& sqrtu)
    {
        u = urbg.u01lopen();
        sqrtu = std::sqrt(u);
    }

    /*
     * Used for flight path sampling
     *
     * P(\ell) = e^{-\ell/\ell_0}
     *
     * \ell \propto -log(u)
     *
     * The \sqrt{\ell} is also needed for the impact factor
     */
    virtual void poisson(float& u, float& sqrtu)
    {
        u = urbg.u01open();
        u = -std::log(u);
        sqrtu = std::sqrt(u);
    }
    virtual void poisson(float& u)
    {
        u = urbg.u01open();
        u = -std::log(u);
    }
    /*
     * Used for sampling of straggling energy
     */
    virtual float normal()
    {
        return N(urbg);
    }
    /*
     * Samples uniformly a random 2D direction
     * nx, ny are the direction cosines
     *
     * Used for sampling the polar angle after collision
     */
    virtual void azimuth(float& nx, float&ny)
    {
        float s1,s2,r2;
        do
        {
            nx = 2.f*urbg.u01() - 1.f;
            ny = 1.f*urbg.u01();
            s1 = nx*nx;
            s2 = ny*ny;
            r2 = s1 + s2;
        } while (r2>1.f || r2==0.f);
        ny = 2*nx*ny/r2;
        nx = (s1-s2)/r2;
    }
};

struct random_vars_tbl : public random_vars
{

    random_vars_tbl(URBG& g) : random_vars(g)
    {
        initTables();
    }

    virtual void u_sqrtu(float& u, float& sqrtu) override
    {
        u = randomlist[iUniform];
        sqrtu = sqrtrandomlist[iUniform++];
        iUniform %= randomlist.size();
    }
    virtual void poisson(float& u) override
    {
        u = sqrtloglist[iSqrtLog++];
        u *= u;
        iSqrtLog %= sqrtloglist.size();
    }
    virtual void poisson(float& u, float& sqrtu) override
    {
        sqrtu = sqrtloglist[iSqrtLog++];
        u = sqrtu*sqrtu;
        iSqrtLog %= sqrtloglist.size();
    }
    virtual void azimuth(float& nx, float&ny) override
    {
        ny = sinAzimAngle[iAzimuth];
        nx = cosAzimAngle[iAzimuth++];
        iAzimuth %= sinAzimAngle.size();

//        const int mask = (1 << 17) - 1;
//        int i = urbg() & mask;
//        ny = sinAzimAngle[i];
//        nx = cosAzimAngle[i];
    }

private:
    void initTables();

    /*************** Adapted from corteo.h ***********************/
    std::vector<float> randomlist;     /* list of evenly distributed but randomly ordered values between 0 and 1 */
    std::vector<float> sqrtrandomlist; /* sqrt of randomlist */
    std::vector<float> sqrtloglist;    /* list of evenly distributed but randomly ordered values of sqrt of -log of 1/MAXLOGLIST to 1 */
    std::vector<float> sqrtloglist1;   /* 1/sqrtloglist */
    std::vector<float> sinAzimAngle;   /* list cos and sin components of angles... */
    std::vector<float> cosAzimAngle;   /*   ...this angle are evenly distributed but randomly ordered between 0 and 2*PI */

    int iUniform, iSqrtLog, iAzimuth;
};

#endif // RANDOM_VARS_H
