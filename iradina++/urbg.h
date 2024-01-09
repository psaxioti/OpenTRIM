#ifndef URBG_H
#define URBG_H

#include <random>

struct URBG : public  std::minstd_rand // 1: std::mt19937 2: std::minstd_rand (20 % faster)
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
