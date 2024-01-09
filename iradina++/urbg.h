#ifndef URBG_H
#define URBG_H

#include <random>

struct URBG : public std::mt19937
{
    float u01() { return 1.f*(*this)()/(*this).max(); }
    float u01no1() {
        float u;
        do u = u01(); while(u==1.f);
        return u;
    }
    float u01no0no1() {
        float u;
        do u = u01(); while(u==1.f || u==0.f);
        return u;
    }
    float u01no0() {
        float u;
        do u = u01(); while(u==0.f);
        return u;
    }
};

#endif // URBG_H
