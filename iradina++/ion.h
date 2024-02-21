#ifndef _ION_H_
#define _ION_H_

#include <cmath>

#include "geometry.h"

class atom;


class ion
{
    vector3 pos_; // position = x,y,z in nm
    vector3 dir_; // direction cosines
    float erg_; // energy in eV
    ivector3 icell_;
    int cellid_;
    int ion_id_;
    int recoil_id_;
    const atom* atom_;
    const grid3D* grid_;

public:

    ion() :
        pos_(0.f,0.f,0.f),
        dir_(0.f,0.f,1.f),
        erg_(1.f),
        icell_(), cellid_(-1),
        ion_id_(0), recoil_id_(0),
        atom_(nullptr), grid_(nullptr)
    {}

    ion(const ion& i) :
        pos_(i.pos_), dir_(i.dir_), erg_(i.erg_),
        icell_(i.icell_), cellid_(i.cellid_),
        ion_id_(i.ion_id_), recoil_id_(i.recoil_id_),
        atom_(i.atom_), grid_(i.grid_)
    {}

    ion& operator=(const ion& i) {
        pos_ = i.pos_; dir_ = i.dir_; erg_ = i.erg_;
        icell_ = i.icell_; cellid_ = i.cellid_;
        ion_id_ = i.ion_id_; recoil_id_ = i.recoil_id_;
        atom_ = i.atom_; grid_ = i.grid_;
        return *this;
    }

    const vector3& pos() const { return pos_; }
    const vector3& dir() const { return dir_; }
    float erg() const { return erg_; }
    const ivector3& icell() const { return icell_; }
    int cellid() const { return cellid_; }
    int ion_id() const { return ion_id_; }
    int recoil_id() const { return recoil_id_; }
    const atom* myAtom() const { return atom_; }

    int setPos(vector3 x);
    vector3& dir() { return dir_; }
    float& erg() { return erg_; }
    int& ion_id() { return ion_id_; }
    int& recoil_id() { return recoil_id_; }
    const atom*& myAtom() { return atom_; }
    void setGrid(const grid3D* g) { grid_ = g; }

    void deflect(vector3 n);
    float distance2boundary(const box3D& b);
    void propagate(const float& s);
};

#endif // ION_H
