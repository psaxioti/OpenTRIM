#ifndef _CORTEO_H_
#define _CORTEO_H_

#include <cassert>
#include <stdexcept>

/*********************************************************************/
/* This file contains definitions adapted from Corteo (version
   Corteo20090527).
   The index values defined in this file can be used with the scattering
   matrix indexed with 4 mantissa bits in the reduced energy and the
   reduced impact factor.

   Corteo was written by Francois Schiettekatte.
   Corteo was released under the GNU General Public License as
   published by the Free Software Foundation, version 3.
   
   You may obtain the original corteo source code from:
   http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
*/
/*********************************************************************/

    /*************** Adapted from corteoindex.c ***********************/

    /* Xindex: compute an integer index using the binary representation of a floating point value val/MIN
    (X is E for reduced energy, S for reduced impact parameter or D for energy of stopping).
    These indexes are used to access sin/cos tables of the scattering angle and stopping power tables.
    The division by MIN is actually carried out by subtracting the appropriate BIAS from the resulting index

    Input: val, with val >= MIN, (MIN = MINE, MINS, or MIND)

    Returns an integer equal to 16 times the exponent base 2 of val 
    + 1 for each 1/16 interval between each power of 2
    (Index of input value = MIN is 0)

    examples
    val/MIN    index
    1           0
    1.0625      1
    1.0626      1
    1.1249      1
    1.1250      2
    1.2         3
    1.25        4
    2           16
    5           40

    BEWARE: architecture dependent!
    Assuming 32-bit 'unsigned long' integer and float
    Assuming IEEE Standard 754 representation for float:
    bit 31: sign
    bit 30-23: exponent (base 2), biased by 127 (i.e. value of exponent for 2^0 is 127)
    bit 22-0: mantissa
    */

/* matrix index calculation parameters */
//#define MINE  0.0000019073486328125f  //  minimum reduced energy = 2^(-19)
//#define MAXE  2097152.f // maximum reduced energy = 2^21
//#define DIME  640    // matrix dimension: (21-(-19))*2^4 (4 mantissa bits: *2^4)
//#define BIASE 1728   // (127-19)*2^4: exponent bias correction while performing division by 2^(-4)
//#define SHIFTE 19    // keep 4 of the 23 mantissa bits (23-4=19)

//#define MINS  0.00000001490116119384765625f // mimimum reduced impact parameter = 2^(-26)
//#define MAXS  64.f   // maximum reduced impact parameter = 2^6
//#define DIMS  512    // (6-(-26))*2^4
//#define BIASS 1616   // (127-26)*2^4: exponent bias correction while performing division by 2^(-16)
//#define SHIFTS 19    // keep 4 of the 23 mantissa bits (23-4=19)




/**
 * @brief A fast indexing method for log-spaced ranges of real numbers
 *
 * This is a fast indexing method for log-spaced tables of real numbers based
 * on the widely used floating number representation in the IEEE-754 standard.
 *
 * It was originally proposed by Yuan et al NIMB83(1993)p.413 for tabulating scattering
 * cross sections.
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
 * \f[
 *   X_I = 2^{\alpha_I + E_{min}} \times \left[ 1 + \beta_I \right] \\
 *   \alpha_I = \frac{I}{2^N} - \beta_I \\
 *   \beta_I = \frac{I \bmod 2^N}{2^N} \\
 *   D = 2^N\left(E_{max} - E_{min}\right)
 * \f]
 *
 * The form of \f$ X_I \f$ is very similar to the IEEE-754 representation of real numbers used in most computers. Thus, \f$ I \f$ can be converted
 * to \f$ X_I \f$ and vice-versa by fast bit manipulation operations.
 *
 * \f$ N \f$ is the number of bits used to represent the mantissa (this is the \f$ 1 + \beta_I \f$ part). IEEE-754 32-bit floats have \f$ N=23 \f$.
 *
 * The corteo_index object can be used like a c++ iterator to loop over the corteo range and provides access to the index and the floating point value.
 *
 * For example, a typical for-loop could be
 * \code
 * corteo_index<4, -5, 5> i; // i initialized to 0
 * for(; i!=i.end(); i++) {
 *     A[i] = ... // array indexing OK - implicit conversion to int
 *     float f = *i; // dereference operator returns the floating point number
 * }
 * \endcode
 *
 * 
 * @tparam _Nb # of bits  for the mantissa
 * @tparam _minExp Exponent of the min value = 2^minExp
 * @tparam _maxExp Exponent of the max value = 2^maxExp
 * @tparam iT integral type for the index, typically int
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
        ll = (index+bias+1) << shift;
        v += *reinterpret_cast<float*>(&ll);

        return v*0.5f;
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
    corteo_index begin() const { return corteo_index(iT(0));   }

    /**
     * @brief Returns an index to 1 past the end of the range, i.e. D
     * @return A corteo_index pointing to 1 past the end
     */
    corteo_index end()   const { return corteo_index(iT(dim)); }

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
     * @return true is the machine is compatible
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

struct corteo6bit {
    // Corteo 6-bit indexes
    // reduced energy index 6bit
    typedef corteo_index<6, -19, 21, int> e_index;
    // reduced impact parameter index 6bit
    typedef corteo_index<6, -26, 6, int> s_index;
    const static int rows = e_index::dim;
    const static int cols = s_index::dim;
    static int table_index(const float& e, const float& s) {
        return e_index(e)*cols + s_index(s);
    }
};

struct corteo4bit {
    // Corteo 4-bit indexes
    // reduced energy index
    typedef corteo_index<4, -19, 21> e_index;
    // reduced impact parameter index
    typedef corteo_index<4, -26,  6> s_index;
    const static int rows = e_index::dim;
    const static int cols = s_index::dim;
    static int table_index(const float& e, const float& s) {
        return e_index(e)*cols + s_index(s);
    }
};



/* These  have been checked for compile time evaluation & correctness
static_assert(reduced_energy_index::minVal == MINE, "values are not equal!");
static_assert(reduced_energy_index::maxVal == MAXE, "values are not equal!");
static_assert(reduced_energy_index::dim == DIME, "values are not equal!");
static_assert(reduced_energy_index::bias == BIASE, "values are not equal!");
static_assert(reduced_energy_index::shift == SHIFTE, "values are not equal!");

static_assert(reduced_s_index::minVal == MINS, "values are not equal!");
static_assert(reduced_s_index::maxVal == MAXS, "values are not equal!");
static_assert(reduced_s_index::dim == DIMS, "values are not equal!");
static_assert(reduced_s_index::bias == BIASS, "values are not equal!");
static_assert(reduced_s_index::shift == SHIFTS, "values are not equal!");

static_assert(energy_index::minVal == MIND, "values are not equal!");
static_assert(energy_index::maxVal == MAXD, "values are not equal!");
static_assert(energy_index::dim == DIMD, "values are not equal!");
static_assert(energy_index::bias == BIASD, "values are not equal!");
static_assert(energy_index::shift == SHIFTD, "values are not equal!");
*/

#endif
