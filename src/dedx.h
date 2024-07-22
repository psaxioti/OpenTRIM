#ifndef _DEDX_H_
#define _DEDX_H_

#include "corteo.h"

/**
 * \defgroup dedx libiondedx shared library
 *
 * @brief Tables of ion electronic energy loss
 *
 * @{
 *
 * The tables have been compiled with the program
 * `SRmodule.exe`
 * distributed with SRIM-2013 (http://www.srim.org).
 *
 * \f$dE/dx\f$ is given as a function of ion energy on a log-spaced corteo range
 * defined by the class \ref dedx_index.
 * 
 * Tables are provided for all projectile (\f$Z=Z_1\f$) / target (\f$Z=Z_2\f$) compinations with \f$ 1 \leq Z_1, Z_2 \leq 92\f$ 
 *
 * The data are compiled into a dynamic library (libdedx.so or .dll).
 * Access is provided by the function \ref dedx().
 *
 * @}
 *
 */

/**
 * @brief A 4-bit corteo_index for tables of ion energy loss dE/dx
 *
 * This index provides fast access to a log-spaced table of ion energy values.
 * They are intended for interpolating tables of energy loss, dE/dx.
 *
 * The number of tabulated values is (30-4)*2^4 + 1 = 417,
 * in the range \f$ 2^4 = 16 \leq E \leq 2^{30} \sim 10^9 \f$ in [eV].
 *
 * @ingroup dedx
 *
 */
typedef corteo_index<float, int, 4, 4, 30> dedx_index;

/**
 * @brief Return a table of electronic energy loss for a projectile (atomic number Z1) moving inside a target (atomic number Z2)
 *
 * The energy loss is calculated at log-spaced ion energy values given by \ref dedx_index.
 *
 * The tables have been generated for all \f$ (Z_1,Z_2) \f$ combinations,
 * with \f$  1 \leq Z_{i=1,2} \leq 92 \f$,
 * using the free program SRModule provided by SRIM-2013.
 *
 * The energy loss is in units of eV / [10^15 at/cm^2] = 10^-15 eV-cm^2.
 *
 * Multiply by the target atomic density to obtain the stopping power dE/dx in eV/cm.
 *
 * @param Z1 projectile atomic number
 * @param Z2 target atomic number
 * @return pointer to the first data point in the table
 *
 * @ingroup dedx
 */
const float* dedx(int Z1, int Z2);

#endif // DEDX_H
