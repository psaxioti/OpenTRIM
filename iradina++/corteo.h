#ifndef _CORTEO_H_
#define _CORTEO_H_

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

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
 * @brief An integer indexing of IEEE floating point numbers
 * 
 * @tparam _Nb # of bits  for the mantissa
 * @tparam _minExp Exponent of the min value = 2^minExp
 * @tparam _maxExp Exponent of the max value = 2^maxExp
 */
template<int _Nb, int _minExp, int _maxExp, class iT> 
struct corteo_index 
{
    // these constants are computed at compile time
    // this is ensured by the constexpr

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

    static const int Nbits = _Nb;

private:
    /* get index from float value */
    /**
     * @brief Convert float value to index
     * 
     * No checking is performed wether val is within the range
     * 
     * @param val a floating point value
     * @return unsigned int the corresponding index
     */
    static iT val2idx(float val) {
        assert(val >= minVal);
        iT ll = *reinterpret_cast<iT*>(&val);
        ll = (ll >> shift)-bias;
        assert(ll < dim);
        return ll;
    }

    /* get value from index */
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

    corteo_index(iT i = 0) : i_(i) {}
    explicit corteo_index(float v) : i_(val2idx(v)) {}
    static corteo_index fromValue(float v) { return corteo_index(val2idx(v)); }
    corteo_index& operator++() { i_++; return *this; }
    corteo_index operator++(int) { corteo_index retval = *this; ++(*this); return retval; }
    bool operator==(corteo_index other) const { return i_ == other.i_; }
    bool operator!=(corteo_index other) const { return !(*this == other); }

    corteo_index begin() const { return corteo_index(iT(0)); }
    corteo_index end() const { return corteo_index(iT(dim)); }

    float operator*() { return idx2val(i_); }

    operator iT() const { return i_; }

    float value() const { return idx2val(i_); }
    iT index() const {return i_; }
private:
    iT i_;
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

struct xs_base;

class corteo {

public:
#ifdef CORTEO_6_BIT
    // reduced energy index
    typedef corteo_index<6, -19, 21, int> epsilon_index;
    // reduced impact parameter index
    typedef corteo_index<6, -26, 6, int> s_index;
#else
    // reduced energy index
    typedef corteo_index<4, -19, 21, int> epsilon_index;
    // reduced impact parameter index
    typedef corteo_index<4, -26, 6, int> s_index;
#endif

protected:

    /* scattering matrix: sin^2(theta_CM/2) */
    std::vector<float> matrix;

public:
    corteo();

    static int check_type_representation();

    float operator()(const epsilon_index& ie, const s_index& is) const
    {
        return matrix[int(is)+int(ie)*s_index::dim];
    }
//    float operator()(int ie, int is) const
//    {
//        return matrix[is+ie*s_index::dim];
//    }
    float operator()(float epsilon, float s) const
    {
        return (*this)(epsilon_index(epsilon),s_index(s));
    }

    static unsigned int table_index(float epsilon, float s)
    {
        //return (int)(s_index(s)) +  ((int)epsilon_index(epsilon))*s_index::dim;
        return s_index(s) +  epsilon_index(epsilon)*s_index::dim;
    }

    const float* data() const { return matrix.data(); }

    int size_epsilon() const { return epsilon_index::dim; }
    int size_s() const { return s_index::dim; }
    int size() const { return s_index::dim * epsilon_index::dim; }

    /* compute all the elements of the matrix
    user sets showProgress!=0 to display the progress of this (long) calculation to the console
    return 1 if successful, 0 if not able to write file */
    /*
     * compute all the elements of the xs matrix
     *
     * if a ostream pointer is passed, msgs are written showing the progress
     *
     * return 1 if successful, 0 if not able to write file
     *
     */

    int calcMatrix(const xs_base& xs, std::ostream* os = NULL);

};







#endif
