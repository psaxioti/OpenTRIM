#ifndef _RANDOM_VARS_H_
#define _RANDOM_VARS_H_

#include <random>
#include <array>
#include <cstdint>

/**
 * \defgroup RNG Random numbers
 * @brief Classes for generating random quantities needed in the Monte-Carlo simulation.
 *
 * @{
 *
 * @ingroup MC
 *
 * @}
 *
 *
 *
 */


/**
 * @brief Xoshiro128+ random number generator
 *
 * Adopted from the original by David Blackman and Sebastiano Vigna (vigna@acm.org) 2018
 *
 * http://prng.di.unimi.it/xoshiro128plus.c
 *
 * > This is xoshiro128+ 1.0, our best and fastest 32-bit generator for 32-bit
 * > floating-point numbers. We suggest to use its upper bits for
 * > floating-point generation, as it is slightly faster than xoshiro128**.
 * > It passes all tests we are aware of except for
 * > linearity tests, as the lowest four bits have low linear complexity, so
 * > if low linear complexity is not considered an issue (as it is usually
 * > the case) it can be used to generate 32-bit outputs, too.
 * >
 * > We suggest to use a sign test to extract a random Boolean value, and
 * > right shifts to extract subsets of bits.
 * >
 * > The state must be seeded so that it is not everywhere zero."
 *
 *  The C++ implementation here  qualifies as a std::uniform_random_bit_generator
 *  and can be used with std lib funtions.
 *
 *  - Output type: 32bit
 *  - Period: 2^128 - 1
 *  - State: 16 bytes
 *  - Version 1.0
 *
 *  @ingroup RNG
 */
class Xoshiro128Plus
{
    static constexpr std::uint32_t RotL(const std::uint32_t x, const int s) noexcept
    {
        return (x << s) | (x >> (32 - s));
    }

public:

    using state_type	= std::array<std::uint32_t, 4>;
    using result_type	= std::uint32_t;

    /**
     * @brief Xoshiro128Plus default constructor
     *
     * The argument is used to initialize the state by calling seed()
     *
     * a std::mt19937 object, which is
     * then used to produce the initial state.
     *
     * If no seed is given then a default seed (0x61884152) is used.
     *
     * @param seed the seed used initialize the state
     */
    explicit Xoshiro128Plus(result_type aseed = DefaultSeed) noexcept
        : m_state()
    {
        seed(aseed);
    }
    /**
     * @brief Construct with a given state
     * @param state the data to initialize the internal state
     */
    explicit constexpr Xoshiro128Plus(const state_type& state) noexcept
        : m_state(state)
    {}

    /**
     * @brief Set the internal state by a seed s
     *
     * s is used to initialize a std::mt19937 object, which is
     * then used to produce the internal state data.
     *
     * @param s the seed used initialize the state
     */
    void seed(result_type s)
    {
        std::mt19937 mt(s);

        for (auto& state : m_state)
        {
            state = static_cast<std::uint32_t>(mt());
        }
    }

    /// Advances the state and returns the generated value
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

    /**
     * @brief Advances the state by 2^64 steps
     *
     * Can be used to generate 2^64 non-overlapping subsequences for parallel computations.
     */
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

    /**
     * @brief Advances the state by 2^96 steps
     *
     * Can be used to generate 2^32 non-overlapping subsequences for parallel computations.
     */
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

    /// Smallest output value (0)
    static inline constexpr result_type min() noexcept
    {
        return result_type(0);
    }

    /// Largest output value (2^32 - 1)
    static inline constexpr result_type max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    /// Return a const reference to the internal state
    inline const state_type& state() const noexcept
    {
        return m_state;
    }

    /// Set the internal state to a new value
    inline constexpr void state(const state_type& s) noexcept
    {
        m_state = s;
    }

    /// Return true if the two objects have equal states
    friend bool operator ==(const Xoshiro128Plus& lhs, const Xoshiro128Plus& rhs) noexcept
    {
        return (lhs.m_state == rhs.m_state);
    }
    /// Return true if the two objects have unequal states
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
 * It is based on the xoshiro128+ random number generator.
 *
 * @ingroup RNG
 *
 */
class random_vars : public Xoshiro128Plus
{
    typedef Xoshiro128Plus rng_engine;

    std::normal_distribution<float> N_;

public:
    /// default constructor
    random_vars() : rng_engine()
    {}
    /// Quasi copy constructor, using an already defined engine
    explicit random_vars(rng_engine& e) : rng_engine(e.state())
    {}

    /// Return random value in [0, 1)

    /**
     * @brief Return a random 32bit float in [0, 1)
     *
     * The function uses the upper 24 bits of the output
     * from xoshiro128+ divided by 2^24.
     *
     * This creates a 32bit float uniformly distributed in [0, 1)
     *
     * The value of 1 is never returned. The maximum value that can be returned
     * is 1.f - std::limits<float>::epsilon()/2
     *
     */
    float u01ropen() {
        std::uint32_t u = (*this)();
        return (u >> 8) * 0x1.0p-24f;
    }
    /// Same as u01ropen()
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

#endif // RANDOM_VARS_H
