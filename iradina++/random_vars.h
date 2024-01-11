#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include <cassert>
#include <random>
#include <vector>

/*
 * Generate random distributions needed in Iradina++
 *
 * Interface definition with pure virtual funcs
 *
 */
struct random_vars_base
{
    virtual ~random_vars_base()
    {}

    /*
     * Used for impact factor sampling
     *
     * $f p \propto \sqrt{u} $f
     *
     */
    virtual void u_sqrtu(float& u, float& sqrtu) = 0;

    /*
     * Used for flight path sampling
     *
     * P(\ell) = e^{-\ell/\ell_0}
     *
     * \ell \propto -log(u)
     *
     * The \sqrt{\ell} is also needed for the impact factor
     */
    virtual void poisson(float& u, float& sqrtu) = 0;
    virtual void poisson(float& u) = 0;

    /*
     * Used for sampling of straggling energy
     *
     * Returns a number distributed as N(0,1)
     */
    virtual float normal() = 0;

    /*
     * Samples uniformly a random 2D direction
     *
     * nx, ny are the direction cosines
     *
     * Used for sampling the azimuthal angle after collision
     */
    virtual void azimuth(float& nx, float&ny) = 0;

};


/*
 * Basic 32bit float random number generator (RNG)
 *
 * Uses internally an integral RNG engine from std lib
 *
 * Tested engines :
 *
 * 1. mt19937 : good overall engine but slower
 *
 * 2. minstd_rand : "minimal standard", low resolution but fast
 *
 */
template<class _E>
struct URBG_ : public _E
{
    typedef _E rng_engine;

    explicit URBG_(_E& e) : _E(e)
    {}
    URBG_() : _E()
    {}

    // [0, 1)
    float u01ropen() {
        return std::generate_canonical<float, 32>(*this);
    }
    // (0, 1)
    float u01open() {
        float u;
        do u = u01ropen(); while(u==0.f);
        return u;
    }
    // (0, 1]
    float u01lopen() {
        return 1.f - u01ropen();
    }
    // [0, 1]
    // experimental - needs testing
    float u01() { return 1.f*(*this)()/rng_engine::max(); }
};

typedef URBG_< std::mt19937 > URBGmt;
typedef URBG_< std::minstd_rand > URBGmsrand;

/*
 * This object use pure random generation
 */
template<class _E>
struct random_vars : public random_vars_base
{
    typedef URBG_<_E> _my_urbg_t;

    _my_urbg_t& urbg;

    std::normal_distribution<float> N;

    random_vars(_my_urbg_t& e) : urbg(e)
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

/*
 * This object use tabulated random numbers
 */
template<class _E>
struct random_vars_tbl : public random_vars_base
{
    typedef URBG_<_E> _my_urbg_t;

    _my_urbg_t& urbg;

    random_vars_tbl(_my_urbg_t& e) : urbg(e)
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
    virtual float normal() override
    {
        float v = inverflist[iInvErf++];
        iInvErf %= inverflist.size();
        return v;
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
    std::vector<float> inverflist;     /* list of evenly distributed but randomly ordered values of inverse erf of -1+1/MAXERFLIST to 1-1/MAXERFLIST */


    int iUniform, iSqrtLog, iAzimuth, iInvErf;
};

#endif // RANDOM_VARS_H
