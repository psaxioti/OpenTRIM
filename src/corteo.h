#ifndef _CORTEO_H_
#define _CORTEO_H_

#include <cassert>
#include <limits>
#include <cfloat>

/**
 * \defgroup CorteoIdx Corteo indexing
 *
 * \brief A fast indexing method for log-spaced ranges.
 * @{
 *
 * It was originally proposed by Yuan et al NIMB83(1993) p.413 for indexing
 * of tabulated
 * scattering cross sections over wide ranges of energy and impact parameter.
 *
 * The implentation here is based on the program Corteo, written by Francois Schiettekatte
 * and released under the GNU GPL.
 * The original corteo source code can be found here:
 * http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
 *
 * The indexing method utilizes the binary floating-point number representation
 * used in most computer systems which is based on the IEC559/IEEE754 standard.
 *
 * The IEEE754 float has the form
 * \f[
 * (-1)^s \times m \times 2^E
 * \f]
 * where \f$ s \f$ is the sign bit, \f$ m \f$ is the "mantissa" and \f$ E \f$ the exponent in the
 * range \f$ [1-E_{max}, E_{max}] \f$.
 *
 * For the subset of "normal" floats
 * \f[
 * m = 1 + k/2^n : 0 \leq k < 2^n
 * \f]
 * where \f$ n \f$ is the number of mantissa bits.
 *
 * The 32-bit IEEE754 float has \f$ E_{max}=\f$
 * FLT_MAX_EXP - 1 = 127  and \f$ n=\f$
 * FLT_MANT_DIG - 1 = 23.
 *
 * The Corteo indexing specifies a correspondence between an integer sequence
 * \f$ I = 0, 1, 2, \dots D \f$
 * and a range of real numbers \f$ X_I \f$, \f$ 2^{E_{min}} \leq X_I \leq 2^{E_{max}} \f$,
 * where
 * \f{eqnarray*}{
 *   X_I &=& m \times 2^E \\
 *   m &=& 1 + (I \bmod 2^N) \times 2^{-N}  \\
 *   E &=& E_{min} + (I \div 2^N) \\
 *   D &=& 2^N\left(E_{max} - E_{min}\right)
 * \f}
 *
 * The form of \f$ X_I \f$ is very similar to the IEEE754 representation with \f$ N \f$ mantissa bits.
 * Thus, \f$ I \f$ can be converted
 * to \f$ X_I \f$ and vice-versa by fast bit manipulation operations.
 *
 * An example with \f$ N=2 \f$, \f$ E_{min}=0 \f$, \f$ E_{max}=2 \f$ is shown in the table:
 *
 *
 * | 2-bit corteo range from 2^0 to 2^4 ||||||||||
 * | :- | :-: | :-:  | :-:  | :-:  | :-:  | :-:  | :-:  | :-:  | :-: |
 * | \f$ I \f$  | 0   | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8   |
 * | \f$ X_I \f$ | 1   | 1.25 | 1.5 | 1.75 | 2 | 2.5 | 3 | 3.5 | 4 |
 *
 *
 * Thus, it is evident that this type of indexing can be used for fast lookup of logarithmic tables, e.g., in interpolation schemes, without the need
 * to call the log() function.
 *
 * Tabulating a scattering cross-section over a wide range of particle energies spanning several orders of magnitude
 * is such an application.
 *
 * @}
 *
 */

// get templated function for the # of mantissa digits
template<class T> struct num_detail;
template<>
struct num_detail<float>
{
    constexpr static int mantissa_digits() { return FLT_MANT_DIG; }
};
template<>
struct num_detail<double>
{
    constexpr static int mantissa_digits() { return DBL_MANT_DIG; }
};
template<>
struct num_detail<long double>
{
    constexpr static int mantissa_digits() { return LDBL_MANT_DIG; }
};


/**
 * @brief The corteo_index structure template implements corteo indexing in C++
 *
 * The corteo_index object can be used like a C++ iterator to loop over
 * a log-spaced range.
 *
 * At the same time, it provides access to the floating point value by means
 * of the dereference operator.
 *
 * For example, a typical for-loop over the corteo range can be implemented
 * as follows:
 *
 * \code
 * corteo_index<float, int, 4, 0, 1> i; // i initialized to 0
 * for(; i<i.end(); i++) {
 *     float f = *i; // dereference operator returns the floating point number
 *     std::cout << f << std::endl;
 *     A[i] = ... // array indexing OK - implicit conversion to int
 * }
 * \endcode
 *
 * This iterates over the defined corteo range and prints
 * all 2^4+1=17 float values from 2^0 = 1 to 2^1 = 2.
 *
 * The following can also be used:
 * \code
 * corteo_index<float, int, 4, 0, 2> R;
 * for(const float& f : R) {
 *     std::cout << f << std::endl;
 * }
 * \endcode
 *
 * The value of std::numeric_limits::is_iec559
 * defined in the <limits> standard library header
 * is checked at compile-time for compatibility
 * of the underlying implementation with the IEEE-754 std.
 *
 * The real and integral types must be of the same size. Thus, the
 * combinations float-int and double-long will generally work.
 *
 * @todo
 * Include a compile-time C++ check for endianess.
 * Currently gcc in Ubuntu 22.04 does not have std::endian
 *
 *
 * @tparam RealType floating-point type
 * @tparam IntType signed integral type for the index
 * @tparam _Nb # of mantissa bits
 * @tparam _minExp Exponent of the min value = 2^minExp
 * @tparam _maxExp Exponent of the max value = 2^maxExp
 *
 * @ingroup CorteoIdx
 */
template<class _RealType,
         class _IntType,
         _IntType _Nb,
         _IntType _minExp,
         _IntType _maxExp>
struct corteo_index
{
    typedef _RealType RealType;
    typedef _IntType IntType;

    // compile-time checks
    static_assert(std::numeric_limits<RealType>::is_iec559,
                  "floating-point type is not IEC559/IEEE754 conformant");
    static_assert(_maxExp > _minExp,
                  "_maxExp must be larger than _minExp");
    static_assert(sizeof(RealType) == sizeof(IntType),
                  "real & int types must have the same size");

    // these constants are computed at compile time
    // this is ensured by the constexpr

    /// @brief Number of mantissa bits
    constexpr static const IntType Nbits = _Nb;
    /// @brief Minimum real value in the range =  2^minExp
    constexpr static const RealType minVal =
        (_minExp >= IntType(0)) ? IntType(1) << _minExp : RealType(1)/(IntType(1) << -_minExp);
    /// @brief Maximum real value in the range = 2^maxExp
    constexpr static const RealType maxVal =
        (_maxExp>=IntType(0)) ? IntType(1) << _maxExp : RealType(1)/(IntType(1) << -_maxExp);
    /// @brief Size of indexing range
    constexpr static const IntType size = (_maxExp - _minExp)*(IntType(1) << _Nb) + IntType(1);


private:
    /// @brief Exponent bias
    constexpr static const IntType bias =
        (IntType(std::numeric_limits<RealType>::max_exponent - 1) + _minExp)*(IntType(1) << _Nb);
    /// @brief Exponent shift
    constexpr static const IntType shift =
        (IntType( num_detail<RealType>::mantissa_digits() - 1 ) - _Nb);
    /// @brief Size of indexing range
    constexpr static const IntType dim = (_maxExp - _minExp)*(IntType(1) << _Nb);


    /**
     * @brief Convert float value to index
     * 
     * The returned value is truncated to within the corteo range
     * 
     * @param val a floating point value
     * @return unsigned int the corresponding index
     */
    static IntType val2idx(RealType val) {
        if (val <= minVal) return IntType(0);
        if (val >= maxVal) return dim;
        IntType ll = *reinterpret_cast<IntType*>(&val);
        ll = (ll >> shift)-bias;
        return ll;
    }

    /**
     * @brief Convert an index to a floating point value
     * 
     * @param index 
     * @return float 
     */
    static RealType idx2val(IntType index) {

        assert(index >= 0 && index <= dim);
        IntType ll = (index+bias) << shift;
        return *reinterpret_cast<RealType*>(&ll);
    }

public:
    /**
     * @brief Constructor with index initializer.
     * @param i The initial value of the index, defaults to 0
     */
    corteo_index(IntType i = IntType(0)) : i_(i) {}

    /**
     * @brief Constructs an index that corresponds to float v
     *
     * The index corresponds to fromValue() called with argument v
     *
     * @param v The floating point number
     */
    explicit corteo_index(const RealType& v) : i_(val2idx(v)) {}

    /**
     * @brief Helper function to convert a real number to index
     *
     * The function finds the real number X_I within the corteo range
     * which is closest to v and returns the corresponding index I.
     *
     * See \ref CorteoIdx for more details.
     *
     * If v is smaller than corteo_index::minVal the function returns begin().
     *
     * If v is larger than corteo_index::maxVal the function returns end().
     *
     * @param v is the floating-point number
     * @return the corresponding index
     */
    static corteo_index fromValue(RealType v) { return corteo_index(val2idx(v)); }

    /**
     * @brief Advance the index by one
     * @return A reference to the index
     */
    corteo_index& operator++() { i_++; return *this; }
    corteo_index operator++(int) { corteo_index retval = *this; ++(*this); return retval; }

    /**
     * @brief Returns an index to the start of the range, i.e. 0
     * @return A corteo_index pointing to 0
     */
    constexpr corteo_index begin() const { return corteo_index(IntType(0));   }

    /**
     * @brief Returns an index to one past the last point of the range, i.e. corteo_index::size
     * @return A corteo_index pointing to one past the end of the range
     */
    constexpr corteo_index end()   const { return corteo_index(IntType(size)); }

    /**
     * @brief Returns an index to the last point of the range, i.e. corteo_index::size-1
     * @return A corteo_index pointing to the last point
     */
    constexpr corteo_index rbegin() const { return corteo_index(IntType(dim));   }

    /**
     * @brief Returns an index to the position before the first point, i.e. -1
     */
    constexpr corteo_index rend()   const { return corteo_index(IntType(-1)); }

    /**
     * @brief The dereference operator * returns the corresponding real value
     * @return The real value corresponding to the corteo index
     */
    RealType operator*() { return idx2val(i_); }

    /**
     * @brief Implicit conversion to integral type IntType (typically int)
     */
    operator IntType() const { return i_; }

    /**
     * @brief Returns the corresponding real value
     * @return The real value corresponding to the corteo index
     */
    RealType value() const { return idx2val(i_); }

    /**
     * @brief Returns the index as an IntType
     */
    IntType index() const {return i_; }

private:
    IntType i_;
};

#endif
