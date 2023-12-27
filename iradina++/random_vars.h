#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include <random>
#include <vector>

typedef std::mt19937 URBG;

struct random_vars
{
    URBG& urbg;

    random_vars(URBG& g) : urbg(g)
    {}

    virtual void uniform01(float& u, float& sqrtu)
    {
        u = 1.f*urbg()/urbg.max();
        sqrtu = std::sqrt(u);
    }
    virtual void sqrtLog(float& u, float& invu)
    {
        do u = 1.f*urbg()/urbg.max(); while(u==0.f);
        u = std::sqrt(-std::log(u));
        invu = 1.f/u;
    }
    virtual void azimuth(float& nx, float&ny)
    {
        float s1,s2,r2;
        do
        {
            nx = 2.f*urbg()/urbg.max() - 1.f;
            ny = 1.f*urbg()/urbg.max();
            s1 = nx*nx;
            s2 = ny*ny;
            r2 = s1 + s2;
        } while (r2>1);
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
        sqrtu = sqrtrandomlist[iUniform];
        iUniform++;
        if (iUniform == randomlist.size()) iUniform = 0;
    }
    virtual void sqrtLog(float& u, float& invu) override
    {
        u = sqrtloglist[iSqrtLog];
        invu = sqrtloglist1[iSqrtLog];
        iSqrtLog++;
        if (iSqrtLog == sqrtloglist.size()) iSqrtLog = 0;
    }
    virtual void azimuth(float& nx, float&ny) override
    {
        ny = sinAzimAngle[iAzimuth];
        nx = cosAzimAngle[iAzimuth];
        iAzimuth++;
        if (iAzimuth == sinAzimAngle.size()) iAzimuth = 0;
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
