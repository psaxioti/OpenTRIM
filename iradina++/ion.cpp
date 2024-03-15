#include "ion.h"


int ion::setPos(vector3 x)
{
    pos_ = x;
    assert(grid_->contains(x));
    icell_ = grid_->pos2cell(x);
    assert(!grid3D::isNull(icell_));
    cellid_ = grid_->cellid(icell_);
    return 0;
}

void ion::propagate(const float& s)
{
    pos_ += s*dir_;  // advance ion position
    if (grid_->contains_with_bc(pos_)) { // is ion still inside the target ?

        if (!grid_->contains(icell_,pos_)) { // did the ion change cell ?
            icell_ = grid_->pos2cell(pos_); // TODO: improve algorithm (e.g. search nn cells)
            //assert(!grid3D::isNull(icell_));
            cellid_ = grid_->cellid(icell_);
        }
    } else {
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
    float smz = 1-dir_.z()*dir_.z();
    if (smz <= 0.f) {
        dir_ = n;
        return;
    }

    vector3 m = dir_;
    smz = std::sqrt(smz);
    dir_.x() = m.x()*n.z() + (m.x()*m.z()*n.x() - m.y()*n.y())/smz;
    dir_.y() = m.y()*n.z() + (m.y()*m.z()*n.x() + m.x()*n.y())/smz;
    dir_.z() = m.z()*n.z() - n.x()*smz;

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
    for(int i=0; i<3; i++) if (dir_[i]<0.f) L[i]=0.f;
    L -= pos_;
    L += b.min();
    L.array() /= dir_.array();
    return L.minCoeff();
}


