#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include <random>

/*
 * Generate random samples of 32bit float numbers
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
template<class rng_engine>
struct random_vars : public rng_engine
{
    std::normal_distribution<float> N;

    explicit random_vars(rng_engine& e) : rng_engine(e)
    {}

    explicit random_vars() : rng_engine()
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
    float u01() {
        return 1.f*(*this)()/rng_engine::max();
    }

    /*
     * Used for sampling of straggling energy
     */
    float normal()
    {
        return N(*this);
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
            nx = 2.f*u01() - 1.f;
            ny = 1.f*u01();
            s1 = nx*nx;
            s2 = ny*ny;
            r2 = s1 + s2;
        } while (r2>1.f || r2==0.f);
        ny = 2*nx*ny/r2;
        nx = (s1-s2)/r2;
    }
};

typedef random_vars< std::mt19937 > rng_mt;
typedef random_vars< std::minstd_rand > rng_msrand;



#endif // RANDOM_VARS_H
