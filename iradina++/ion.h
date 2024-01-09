#ifndef _ION_H_
#define _ION_H_

#include <cmath>

#include "geometry.h"

class atom;


class ion
{
    vector3 pos_; // position = x,y,z, cm
    ivector3 icell_;
    int cellid_;

public:

    vector3 dir; // direction cosines
    float erg; // energy in eV??
    float t; // time in us??
    int ion_id;
    int recoil_id;
    const atom* atom_;

    const grid3D& grid;

    ion(const grid3D& g) : pos_(0.f,0.f,0.f),
        dir(0.f,0.f,1.f),
        erg(1.f),
        t(0.f),
        ion_id(0), recoil_id(0), atom_(0),
        grid(g)
    {}

    ion(const ion& i) : pos_(i.pos_),
        icell_(i.icell_), cellid_(i.cellid_),
        dir(i.dir),
        erg(i.erg), t(i.t),
        ion_id(i.ion_id), recoil_id(i.recoil_id),
        atom_(i.atom_),
        grid(i.grid)
    {}

    ion& operator=(const ion& i) {
        pos_ = i.pos_; icell_ = i.icell_; cellid_ = i.cellid_;
        dir = i.dir; erg = i.erg; t = i.t;
        ion_id = i.ion_id; recoil_id = i.recoil_id;
        atom_ = i.atom_;
        return *this;
    }

    void deflect(vector3 n);

    float distance2boundary(const box3D& b);

    void propagate(const float& s);
    int setPos(vector3 x);
    const vector3& pos() const { return pos_; }
    const ivector3& icell() const {return icell_; }
    int cellid() const { return cellid_; }
};

struct invSqrt {
    virtual float operator()(float val) const { return 1.f/std::sqrt(val); }
};
struct invSqrt_tbl : public invSqrt {
    invSqrt_tbl()
    { mySqrtTableFill(); }
    // 1/sqrt of val is the product of the 1/sqrt of its mantissa and 1/sqrt of its exponent
    // WARNING: approximate solution, precise to 0.0015%, use when precision is not critical
    virtual float operator()(float val) const override
    {
        if ( (*reinterpret_cast<unsigned int *>(&val))==0 )
            return 1.0f/val; // prevent division by 0
        return myInvSqrtTable[ ((*reinterpret_cast<unsigned int *>(&val)) >> 7) & 0xFFFF ]*
               myInvSqrtTableExp[ ((*reinterpret_cast<unsigned int *>(&val)) >> 23) & 0xFF];
    }
private:
    static const std::size_t sz1 = 1 << 8;
    static const std::size_t sz2 = 1 << 16;
    float myInvSqrtTableExp[sz1];
    float mySqrtTableExp[sz1];
    float myInvSqrtTable[sz2];
    float mySqrtTable[sz2];
    // fill tables
    void mySqrtTableFill();
};


#endif // ION_H
