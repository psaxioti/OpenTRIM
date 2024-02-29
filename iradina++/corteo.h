#ifndef _CORTEO_H_
#define _CORTEO_H_

#include <cassert>
#include <stdexcept>

/**
 * \defgroup CorteoIdx Corteo index definitions
 * @{
 *
 * This is a fast indexing method for log-spaced ranges of real numbers based
 * on the widely used IEEE-754 standard for floating-point number representation
 * in computer systems.
 *
 * It was originally proposed by Yuan et al NIMB83(1993) p.413 for tabulating
 * scattering
 * cross sections over wide ranges of energy and impact parameter.
 *
 * The implentation here is based on the program Corteo, written by Francois Schiettekatte
 * and released under the GNU GPL.
 * The original corteo source code can be found here:
 * http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
 *
 * The Corteo indexing specifies a correspondence between an integer sequence
 * \f$ I = 0, 1, 2, \dots D \f$
 * and a range of real numbers \f$ X_I \f$, \f$ 2^{E_{min}} \leq X_I \leq 2^{E_{max}} \f$
 * where
 * \f{eqnarray*}{
 *   X_I &=& 2^{\alpha_I + E_{min}} \times \left[ 1 + \beta_I \right] \\
 *   \alpha_I &=& \frac{I}{2^N} - \beta_I \\
 *   \beta_I &=& \frac{I \bmod 2^N}{2^N} \\
 *   D &=& 2^N\left(E_{max} - E_{min}\right)
 * \f}
 *
 * The form of \f$ X_I \f$ is very similar to the IEEE-754 representation of real numbers used in most computers. Thus, \f$ I \f$ can be converted
 * to \f$ X_I \f$ and vice-versa by fast bit manipulation operations.
 *
 * \f$ N \f$ is the number of bits used to represent the mantissa
 * (this is the \f$ 1 + \beta_I \f$ part).
 *
 * IEEE-754 32-bit floats have \f$ N=23 \f$.
 *
 * Corteo cross-section tables typically use \f$ N=4 \f$ or 6.
 *
 */

/**
 * @brief The corteo_index struct template implements corteo indexing in C++
 *
 * The corteo_index object can be used like a c++ iterator to loop over
 * the corteo range and provides simultaneously access to the index
 * and the floating point value.
 *
 * For example, a typical for-loop over the corteo range can be implemented
 * as follows:
 *
 * \code
 * corteo_index<4, -5, 5> i; // i initialized to 0
 * for(; i!=i.end(); i++) {
 *     A[i] = ... // array indexing OK - implicit conversion to int
 *     float f = *i; // dereference operator returns the floating point number
 * }
 * \endcode
 * 
 * @tparam _Nb # of mantissa bits
 * @tparam _minExp Exponent of the min value = 2^minExp
 * @tparam _maxExp Exponent of the max value = 2^maxExp
 * @tparam iT integral type for the index, defaults to int
 */
template<int _Nb, int _minExp, int _maxExp, class iT = int>
struct corteo_index 
{
    // these constants are computed at compile time
    // this is ensured by the constexpr

    /// @brief # of mantissa bits
    constexpr static const int Nbits = _Nb;
    /// @brief minimum float value =  2^minExp
    constexpr static const float minVal = (_minExp>=0) ? 1 << _minExp : 1.f/(1 << -_minExp);
    /// @brief maximum float value = 2^maxExp
    constexpr static const float maxVal = (_maxExp>=0) ? 1 << _maxExp : 1.f/(1 << -_maxExp);
    /// @brief size of indexing range 0 ... dim
    constexpr static const int dim = (_maxExp - _minExp)*(1 << _Nb);   
    /// @brief exponent bias
    constexpr static const int bias = (127 + _minExp)*(1 << _Nb);
    /// @brief exponent shift
    constexpr static const int shift = (23 - _Nb);



private:

    /**
     * @brief Convert float value to index
     * 
     * The returned value is truncated to within the corteo range
     * 
     * @param val a floating point value
     * @return unsigned int the corresponding index
     */
    static iT val2idx(float val) {
        // assert(val >= minVal);
        if (val <= minVal) return 0;
        iT ll = *reinterpret_cast<iT*>(&val);
        ll = (ll >> shift)-bias;
        if (ll >= dim) return dim-1;
        // assert(ll < dim);
        return ll;
    }

    /**
     * @brief Convert an index to a floating point value
     * 
     * @param index 
     * @return float 
     */
    static float idx2val(iT index) {
        float v;
        iT ll;

        assert(index >= 0 && index<dim);
        ll = (index+bias) << shift;
        v = *reinterpret_cast<float*>(&ll);
        return v;
    }

public:
    /**
     * @brief Constructor with index initializer.
     * @param i The value of the idex, defaults to 0
     */
    corteo_index(iT i = 0) : i_(i) {}

    /**
     * @brief Constructs an index that corresponds to float v
     * @param v The floating point number
     */
    explicit corteo_index(const float& v) : i_(val2idx(v)) {}

    /**
     * @brief Constructs an index that corresponds to float(v)
     * @param v The double floating point number
     */
    explicit corteo_index(const double& v) : i_(val2idx((float)v)) {}

    /**
     * @brief Helper function to convert a real number to index
     * @param v is the floating point number
     * @return the corresponding index
     */
    static corteo_index fromValue(float v) { return corteo_index(val2idx(v)); }

    /**
     * @brief Advance the index by one
     * @return ref to the index
     */
    corteo_index& operator++() { i_++; return *this; }
    corteo_index operator++(int) { corteo_index retval = *this; ++(*this); return retval; }

    /**
     * @brief Compares two indexes
     * @param other is the 2nd index
     * @return true if the indexes are equal
     */
    bool operator==(corteo_index other) const { return i_ == other.i_; }
    bool operator!=(corteo_index other) const { return !(*this == other); }

    /**
     * @brief Returns an index to the start of the range, i.e. 0
     * @return A corteo_index pointing to 0
     */
    constexpr corteo_index begin() const { return corteo_index(iT(0));   }

    /**
     * @brief Returns an index to 1 past the end of the range, i.e. D
     * @return A corteo_index pointing to 1 past the end
     */
    constexpr corteo_index end()   const { return corteo_index(iT(dim)); }

    /**
     * @brief Returns an index to the start of the range, i.e. 0
     * @return A corteo_index pointing to 0
     */
    constexpr corteo_index rbegin() const { return corteo_index(iT(dim-1));   }

    /**
     * @brief Returns an index to 1 past the end of the range, i.e. D
     * @return A corteo_index pointing to 1 past the end
     */
    constexpr corteo_index rend()   const { return corteo_index(iT(-1)); }

    /**
     * @brief operator * returns the corresponding float
     * @return
     */
    float operator*() { return idx2val(i_); }

    /**
     * @brief Implicit conversion to integral type iT (typically int)
     */
    operator iT() const { return i_; }

    float value() const { return idx2val(i_); }
    iT index() const {return i_; }

    /*
     * *** Adapted from corteo20130715 ***
     * Length of float and int must be 4 byte, otherwise indexing doesn't
     * work correctly.
     */
    /**
     * @brief A function for checking machine compatibility
     *
     * This function performs some compile-time and run-time checks to
     * see if the hardware is compatible, i.e. if it follows IEEE-754 std.
     *
     * @return true if the machine is compatible
     */
    static int check_type_representation()
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
private:
    iT i_;
};

/**
 * @brief The corteo4bit struct provides 4-bit indexing for tabulated cross-sections
 *
 * The struct defines 2 corteo indexes:
 *   - corteo4bit::e_index for the reduced energy with 640 points (rows) from \f$ 2^{-19} \f$ to \f$ 2^{21} \f$
 *   - corteo4bit::s_index for the reduced impact parameter with 512 points (columns) from \f$ 2^{-26} \f$ to \f$ 2^{6} \f$
 *
 */
struct corteo4bit {
    /// Corteo 4-bit index for the reduced energy
    typedef corteo_index<4, -19, 21> e_index;
    /// Corteo 4-bit index for the reduced impact parameter
    typedef corteo_index<4, -26,  6> s_index;
    /// number of rows (energy values) = 640
    const static int rows = e_index::dim;
    /// number of columns (impact parameter values) = 512
    const static int cols = s_index::dim;
    /**
     * @brief Calculate the table index for given energy and impact parameter
     *
     * Return the index to the memory location where the value that corresponds
     * to given energy and impact parameter is stored.
     *
     * Tables are stored in C-style column-major order.
     *
     * @param e the reduced energy
     * @param s the reduced impact parameter
     * @return the table index
     */
    static int table_index(const float& e, const float& s) {
        return e_index(e)*cols + s_index(s);
    }
};

/**
 * @brief The corteo6bit struct provides 6-bit indexing for tabulated cross-sections
 *
 * The struct defines 2 corteo indexes:
 *   - corteo6bit::e_index for the reduced energy with 2560 points (rows) from \f$ 2^{-19} \f$ to \f$ 2^{21} \f$
 *   - corteo6bit::s_index for the reduced impact parameter with 2048 points (columns) from \f$ 2^{-26} \f$ to \f$ 2^{6} \f$
 */
struct corteo6bit {
    /// Corteo 6-bit index for the reduced energy
    typedef corteo_index<6, -19, 21, int> e_index;
    /// Corteo 6-bit index for the reduced impact parameter
    typedef corteo_index<6, -26, 6, int> s_index;
    /// number of rows (energy values)
    const static int rows = e_index::dim;
    /// number of columns (impact parameter values)
    const static int cols = s_index::dim;
    /**
     * @brief Calculate the table index for given energy and impact parameter
     *
     * Return the index to the memory location where the value that corresponds
     * to given energy and impact parameter is stored.
     *
     * Tables are stored in C-style column-major order.
     *
     * @param e the reduced energy
     * @param s the reduced impact parameter
     * @return the table index
     */
    static int table_index(const float& e, const float& s) {
        return e_index(e)*cols + s_index(s);
    }
};



/**@}*/

#endif
