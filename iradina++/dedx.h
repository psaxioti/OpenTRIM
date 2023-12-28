#ifndef _DEDX_H_
#define _DEDX_H_

#include "corteo.h"

//
// Adopted from the following CORTEO code
//#define MIND  16.f         // minimum energy (eV) for stopping power tables = 2^4
//#define MAXD  1073741824.f // maximum energy (eV) for stopping power tables = 2^30 ~ 1 GeV
//#define DIMD  416          //  (30-4)*2^4
//#define BIASD 2096         // (127+4)*2^4: exponent bias correction while performing division by 2^10
//#define SHIFTD 19          // keep 4 of the 23 mantissa bits (23-4=19)

// dedx energy index
typedef corteo_index<4, 4, 30, int> dedx_index;


const float* dedx(int Z1, int Z2);

#endif // DEDX_H
