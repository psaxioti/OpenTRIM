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
 * @brief Propagate the ion for a distance s [nm]
 *
 * The function first advances the position of the ion by
 * a distance s along its direction of motion.
 *
 * It then checks if the ion is still inside the simulation volume.
 * If not, the cell id is set to -1 (invalid) and the function returns.
 *
 * Otherwise, the function checks if the ion is
 * still in the same cell.
 * If not, the new cell is found and
 * the index vector and cell id of the ion are updated.
 *
 * @param s the distance to propagate the ion [nm]
 */
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


