#include "ion.h"


int ion::setPos(vector3 x)
{
    pos_ = x;
    assert(grid.contains(x));
    icell_ = grid.pos2cell(x);
    assert(!grid3D::isNull(icell_));
    cellid_ = grid.cellid(icell_);
    return 0;
}

void ion::propagate(const float& s)
{
    pos_ += s*dir;
    assert(pos_.allFinite());
    grid.applyBC(pos_);
    if (grid.periodic_contains(pos_)) {
        assert(grid.contains(pos_));
        if (!grid.contains(icell_,pos_)) {
            icell_ = grid.pos2cell(pos_); // TODO: improve algorithm (e.g. search nn cells)
            assert(!grid3D::isNull(icell_));
            cellid_ = grid.cellid(icell_);
        }
    } else {
        icell_ = grid3D::nullcell();
        cellid_ = -1;
    }
}

/*
 * Deflect the ion by angles theta, phi
 *
 * n = (cos(phi)*sin(th), sin(phi)*sin(th), cos(th))
 *
 */
void ion::deflect(vector3 n)
{
    // if the ion moves parallel to the z-axis
    // then the new direction is n
    float smz = 1-dir.z()*dir.z();
    if (smz <= 0.f) {
        dir = n;
        return;
    }

    vector3 m = dir;
    smz = std::sqrt(smz);
    dir.x() = m.x()*n.z() + (m.x()*m.z()*n.x() - m.y()*n.y())/smz;
    dir.y() = m.y()*n.z() + (m.y()*m.z()*n.x() + m.x()*n.y())/smz;
    dir.z() = m.z()*n.z() - n.x()*smz;

    //if (!dir.allFinite()) {
    //    dir  = n;
    //}

}

/*
 * Find distance traveled by ion until it reaches
 * the rectangular box boundary
 *
 * Assumes that ion is in the box
 *
 * box.min() is the vector at the box origin
 *
 * For each dimension the equation for the distance l is
 * (e.g. for x)
 *
 * lx = [Lx - (x - xo)]/vx  if vx>0
 *
 * OR
 *
 * lx = -(x - xo)/vx  if vx<0
 *
 * Finally l = min{lx, ly, lz}
 *
 */
float ion::distance2boundary(const box3D& b)
{
    vector3 L = b.sizes();
    for(int i=0; i<3; i++) if (dir[i]<0.f) L[i]=0.f;
    L -= pos_;
    L += b.min();
    L.array() /= dir.array();
    return L.minCoeff();
}

// fill sqrt tables
void invSqrt_tbl::mySqrtTableFill() {
    //  int safe;
    unsigned long i, j; // 1<<16: mantissa tables contain 65536 values
    float val;

    // mantissa from 0 to 2 in 65536 steps, giving off precision!
    for(i=0;i<sz2;i++) {
        // generate a float values
        j = (i<<7) + ((1<<23)*127);
        //    for(safe=0;safe<10;safe++);
        val = *reinterpret_cast<float*>(&j);
        // store its sqrt and 1/sqrt in tables
        mySqrtTable[i]    = std::sqrt(val);
        myInvSqrtTable[i] = 1.0f/mySqrtTable[i];
    }

    // exponent from 2^-128 to 2^128
    for(i=0;i<sz1;i++) {
        // generate a float values
        j = i<<23;
        val = *reinterpret_cast<float*>(&j);
        // store its sqrt and 1/sqrt in tables
        mySqrtTableExp[i]    = std::sqrt(val);
        myInvSqrtTableExp[i] = 1.0f/mySqrtTableExp[i];
    }
}
