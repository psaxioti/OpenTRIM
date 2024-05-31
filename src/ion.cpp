#include "ion.h"


int ion::setPos(const vector3 &x)
{
    pos_ = x;
    assert(grid_->contains(x));
    icell_ = grid_->pos2cell(x);
    assert(!grid3D::isNull(icell_));
    cellid_ = grid_->cellid(icell_);
    return 0;
}

/**
 * @brief Propagate the ion for a distance s [nm] taking care of boundary crossing
 *
 * The function first checks if the new ion position would be within the
 * simulation volume (taking periodic bc into account).
 *
 * If not, then it propagates the ion to just beyond its cell boundary and checks
 * if indeed the ion escapes the simulation. If yes then the ion exits. Otherwise, it
 * just changes cell.
 *
 * If the new position is within the simulation volume, the function checks
 * if it is in the same cell.
 * If not, the new cell is found and
 * the index vector and cell id of the ion are updated.
 *
 * In all cases of boundary crossing, s is adjusted to the minimum travel needed to just
 * overcome the boundary.
 *
 * @param s the distance to propagate the ion [nm]
 * @return the type of boundary crossing (none, internal (cell change), external (ion left simulation))
 */
BoundaryCrossing ion::propagate(float& s)
{
    vector3 x = pos_ + s*dir_;  // calc new ion position
    if (grid_->contains_with_bc(x)) { // is the ion still inside the target ?
        if (!grid_->contains(icell_,x)) { // does the ion exit the cell ?
            box3D b = grid_->box(icell_);
            s = distance2boundary(b);
            float ds = std::max(s,1.f)*std::numeric_limits<float>::epsilon()*10;
            ivector3 ix;
            do {
                s += ds;
                x = pos_ + s*dir_;
                ix = grid_->pos2cell(x);
            } while (icell_ == ix); // TODO: improve algorithm (e.g. search nn cells)
            pos_ = x;
            icell_ = ix;
            cellid_ = grid_->cellid(icell_);
            return BoundaryCrossing::Internal;
        } else {
            pos_ = x;
            return BoundaryCrossing::None;
        }
    } else { // ion is bound to exit simulation
        // 1. Reduce s to just cross the boundary
        box3D b = grid_->box(icell_);
        s = distance2boundary(b);
        float ds = std::max(s,1.f)*std::numeric_limits<float>::epsilon()*10;
        s += ds;
        x = pos_ + s*dir_;
        // 2. still exiting ?
        if (!grid_->contains_with_bc(x)) {
            pos_ = x;
            cellid_ = -1;
            return BoundaryCrossing::External;
        } else { // ion is still in target. just change the cell
            ivector3 ix = grid_->pos2cell(x);
            while (icell_ == ix) {
                s += ds;
                x = pos_ + s*dir_;
                ix = grid_->pos2cell(x);
            }
            pos_ = x;
            icell_ = ix; // TODO: improve algorithm (e.g. search nn cells)
            cellid_ = grid_->cellid(icell_);
            return BoundaryCrossing::Internal;
        }
    }
}

/**
 * @brief Changes the direction of the ion
 *
 * The argument is a vector
 * \f[
 * \mathbf{n} = (\cos\phi\,\sin\theta, \sin\phi\,\sin\theta,\cos\theta)
 * \f]
 * where \f$ (\theta,\phi) \f$ are the polar and azimuthal scattering
 *  angles.
 *
 * The new direction is obtained by a well-known result
 * \f{eqnarray*}{
 * m'_x &=& m_x\, n_z + \frac{m_x\, m_z\, n_x - m_y\, n_y}{\sqrt{1-m_z^2}} \\
 * m'_y &=& m_y\, n_z + \frac{m_y\, m_z\, n_x - m_x\, n_y}{\sqrt{1-m_z^2}}  \\
 * m'_z &=& m_z\, n_z - n_x\, \sqrt{1-m_z^2}
 * \f}
 * where \f$ \mathbf{m} \f$ is the ion's vector of direction cosines.
 *
 */
void ion::deflect(const vector3 &n)
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

    // normalize
    dir_.normalize();

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


