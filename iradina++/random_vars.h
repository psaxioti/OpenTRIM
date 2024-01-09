#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include "urbg.h"
#include <vector>

struct random_vars
{
    URBG& urbg;
    std::normal_distribution<float> N;

    random_vars(URBG& g) : urbg(g)
    {}

    virtual void uniform01(float& u, float& sqrtu)
    {
        u = urbg.u01no0();
        sqrtu = std::sqrt(u);
    }
    virtual void sqrtLog(float& u, float& invu)
    {
        u = urbg.u01no0no1();
        u = std::sqrt(-std::log(u));
        invu = 1.f/u;
    }
    virtual void poisson(float& u, float& sqrtu)
    {
        u = urbg.u01no0no1();
        u = -std::log(u);
        sqrtu = std::sqrt(u);
    }
    virtual void poisson(float& u)
    {
        u = urbg.u01no0no1();
        u = -std::log(u);
    }
    virtual float normal()
    {
        return N(urbg);
    }
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
    {}

    virtual void uniform01(float& u, float& sqrtu) override
    {
        u = randomlist[iUniform];
        sqrtu = sqrtrandomlist[iUniform++];
        iUniform %= randomlist.size();
    }
    virtual void sqrtLog(float& u, float& invu) override
    {
        u = sqrtloglist[iSqrtLog];
        invu = sqrtloglist1[iSqrtLog++];
        iSqrtLog %= sqrtloglist.size();
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
