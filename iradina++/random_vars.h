#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include <random>
#include <array>
#include <cstdint>
// #include "../extern/Xoshiro-cpp/XoshiroCpp.hpp"

// xoshiro128+
// Output: 32 bits
// Period: 2^128 - 1
// Footprint: 16 bytes
// Original implementation: http://prng.di.unimi.it/xoshiro128plus.c
// Version: 1.0
class Xoshiro128Plus
{
    static constexpr std::uint32_t RotL(const std::uint32_t x, const int s) noexcept
    {
        return (x << s) | (x >> (32 - s));
    }

public:

    using state_type	= std::array<std::uint32_t, 4>;
    using result_type	= std::uint32_t;


    explicit Xoshiro128Plus(result_type seed = DefaultSeed) noexcept
        : m_state()
    {
        std::mt19937 mt(seed);

        for (auto& state : m_state)
        {
            state = static_cast<std::uint32_t>(mt());
        }
    }

    explicit constexpr Xoshiro128Plus(const state_type& state) noexcept
        : m_state(state)
    {}

    void seed(result_type s)
    {
        std::mt19937 mt(s);

        for (auto& state : m_state)
        {
            state = static_cast<std::uint32_t>(mt());
        }
    }

    constexpr result_type operator()() noexcept
    {
        const std::uint32_t result = m_state[0] + m_state[3];
        const std::uint32_t t = m_state[1] << 9;
        m_state[2] ^= m_state[0];
        m_state[3] ^= m_state[1];
        m_state[1] ^= m_state[2];
        m_state[0] ^= m_state[3];
        m_state[2] ^= t;
        m_state[3] = RotL(m_state[3], 11);
        return result;
    }

    // This is the jump function for the generator. It is equivalent
    // to 2^64 calls to next(); it can be used to generate 2^64
    // non-overlapping subsequences for parallel computations.
    constexpr void jump() noexcept
    {
        constexpr std::uint32_t JUMP[] = { 0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b };

        std::uint32_t s0 = 0;
        std::uint32_t s1 = 0;
        std::uint32_t s2 = 0;
        std::uint32_t s3 = 0;

        for (std::uint32_t jump : JUMP)
        {
            for (int b = 0; b < 32; ++b)
            {
                if (jump & UINT32_C(1) << b)
                {
                    s0 ^= m_state[0];
                    s1 ^= m_state[1];
                    s2 ^= m_state[2];
                    s3 ^= m_state[3];
                }
                operator()();
            }
        }

        m_state[0] = s0;
        m_state[1] = s1;
        m_state[2] = s2;
        m_state[3] = s3;
    }


    // This is the long-jump function for the generator. It is equivalent to
    // 2^96 calls to next(); it can be used to generate 2^32 starting points,
    // from each of which jump() will generate 2^32 non-overlapping
    // subsequences for parallel distributed computations.
    constexpr void longJump() noexcept
    {
        constexpr std::uint32_t LONG_JUMP[] = { 0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662 };

        std::uint32_t s0 = 0;
        std::uint32_t s1 = 0;
        std::uint32_t s2 = 0;
        std::uint32_t s3 = 0;

        for (std::uint32_t jump : LONG_JUMP)
        {
            for (int b = 0; b < 32; ++b)
            {
                if (jump & UINT32_C(1) << b)
                {
                    s0 ^= m_state[0];
                    s1 ^= m_state[1];
                    s2 ^= m_state[2];
                    s3 ^= m_state[3];
                }
                operator()();
            }
        }

        m_state[0] = s0;
        m_state[1] = s1;
        m_state[2] = s2;
        m_state[3] = s3;
    }


    inline constexpr result_type min() noexcept
    {
        return std::numeric_limits<result_type>::lowest();
    }

    inline constexpr result_type max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    inline const state_type& state() const noexcept
    {
        return m_state;
    }

    inline constexpr void state(const state_type& s) noexcept
    {
        m_state = s;
    }

    friend bool operator ==(const Xoshiro128Plus& lhs, const Xoshiro128Plus& rhs) noexcept
    {
        return (lhs.m_state == rhs.m_state);
    }

    friend bool operator !=(const Xoshiro128Plus& lhs, const Xoshiro128Plus& rhs) noexcept
    {
        return (lhs.m_state != rhs.m_state);
    }

private:
    static constexpr result_type DefaultSeed = 0x61884152;

    state_type m_state;
};

/**
 * @brief The random_vars class is used for generating random quantities needed in the simulation
 *
 * The class uses internally an integral random number engine from the
 * C++ standard library.
 * Tested engines are:
 * 1. std::mt19937 : 32-bit Mersenne Twister - good overall engine
 * 2. std::minstd_rand : "Minimal standard", lower statistical quality but slightly faster engine
 *
 * @tparam
 *   rng_engine the std library random number engine used underneath
 *
 * @ingroup MC
 *
 */
template<class rng_engine>
class random_vars : public rng_engine
{
    std::normal_distribution<float> N_;

    // max value returned by engine
    constexpr static const float maxVal_ = float(rng_engine::max());

public:
    /// default constructor
    random_vars() : rng_engine()
    {}
    /// Quasi copy constructor, using an already defined engine
    explicit random_vars(rng_engine& e) : rng_engine(e)
    {}

    /// Return random value in [0, 1)
    float u01ropen() {
        // This does not return 1, garantied
        // return std::generate_canonical<float,32>(*this);
        std::uint32_t u = (*this)();
        return (u >> 8) * 0x1.0p-24f;
    }
    float u01() { return u01ropen(); }
    /// Return random value in (0, 1]
    float u01lopen() {
        return 1.f - u01ropen();
    }
    /// Return random value in (0, 1)
    float u01open() {
        float u;
        do u = u01ropen(); while(u==0.f);
        return u;
    }
    /// Return a random value distributed as N(0,1)
    float normal()
    {
        return N_(*this);
    }
    /**
     * @brief Generate a random azimuthal direction
     *
     * Equivalent to:
     * @code{.cpp}
     * phi = 2*pi*u; nx = sin(phi); ny = cos(phi);
     * @endcode
     *
     * The employed algorithm samples (nx, ny) without calling
     * trigonometric functions.
     *
     * @param nx sin(phi)
     * @param ny cos(phi)
     */
    void azimuth(float& nx, float& ny)
    {
        float s1,s2,r2;
        do
        {
            nx = 2.f*u01ropen() - 1.f;
            ny = u01ropen();
            s1 = nx*nx;
            s2 = ny*ny;
            r2 = s1 + s2;
        } while (r2>1.f || r2==0.f);
        ny = 2*nx*ny/r2;
        nx = (s1-s2)/r2;
    }
};

/// Predefined random_vars with std::mt19937 (32-bit "Mersenne Twister") as base generator
typedef random_vars< std::mt19937 > rng_mt;
/// Predefined random_vars with std::minstd_rand ("Minimum standard") as base generator
typedef random_vars< std::minstd_rand > rng_msrand;
/// Predefined random_vars with std::minstd_rand ("Minimum standard") as base generator
typedef random_vars< Xoshiro128Plus > rng_xoshiro;


#endif // RANDOM_VARS_H
