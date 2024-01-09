#include "dedx.h"
#include "elements.h"

extern const float** dedx_data[];

const float* dedx(int Z1, int Z2)
{
    return (Z1 && Z2 &&
            Z1<=elements::max_atomic_num &&
            Z2<=elements::max_atomic_num) ?
            dedx_data[Z1][Z2] : nullptr;
}




