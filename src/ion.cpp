#include "ion.h"


int ion::setPos(const vector3 &x)
{
    pos_ = x;
    assert(x.allFinite());
    assert(grid_->contains(x));
    icell_ = grid_->pos2cell(x);
    assert(!grid3D::isNull(icell_));
    cellid_ = grid_->cellid(icell_);
    return 0;
}

void ion::reset_counters()
{
    vac_=impl_=repl_=ncoll_=0;
    path_=ioniz_=phonon_=recoil_=0.f;
}

/**
 * @brief Propagate the ion for a distance s [nm] taking care of boundary crossings
 *
 * The function first checks if the new ion position would be within the
 * simulation volume (taking periodic bc into account).
 *
 * If not, then the ion is propagated to just beyond its cell boundary and the
 * function checks again
 * if indeed the ion escapes the simulation. If yes then the ion exits. Otherwise, it
 * just changes cell.
 *
 * If the ion remains inside the simulation volume, the function checks
 * if it is in the same cell.
 * 
 * If not, then the ion is propagated to just beyond its cell boundary,
 * the new cell is found
 * and the index vector and cell id are updated.
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
            // get our distance to the boundary
            float s1 = distance2boundary(grid_->box(icell_), pos_, dir_);
            // slightly increase s so that we get just outside the box
            float ds = std::max(s1,1.f)*std::numeric_limits<float>::epsilon()*100;
            s1 += ds;
            x = pos_ + s1*dir_;
            grid_->apply_bc(x);
            ivector3 ix = grid_->pos2cell(x);
            while (ix == icell_ && s1 < s) {
                s1 += ds;
                x = pos_ + s1*dir_;
                grid_->apply_bc(x);
                ix = grid_->pos2cell(x);
            };
            // new pos and cell
            s = s1;
            path_ += s;
            pos_ = x;
            if (ix != icell_) {
                icell_ = ix;
                prev_cellid_ = cellid_;
                cellid_ = grid_->cellid(icell_);
                return BoundaryCrossing::Internal;
            } else {
                /* This is a very rare case where although grid_->contains(icell_,x)
                 * returned false in the end the particle does not change cell.
                 * This can happen if the following 2 conditions hold simultaneously:
                 *   A) x_i = (pos + s*dir)_i == L1, where L1 is the boundary of the simulation
                 *      volume along the i-th direction
                 *   B) We have periodic boundary conditions along the i-th direction
                 *
                 * Thus, the call to grid_->contains_with_bc(x) returns true, due to the periodic BCs,
                 * BUT the call to grid_->contains(icell_,x) returns false because the ion is exactly at
                 * the boundary and this is considered "outside the cell" (condition= "x0 <= x < x1").
                 * Note that contains(icell_,x) does not check for periodic BCs.
                 */
                return BoundaryCrossing::None;
            }
        } else { // we remain in the cell
            path_ += s;
            pos_ = x;
            return BoundaryCrossing::None;
        }
    } else { // ion is bound to exit simulation
        // 1. Reduce s to just cross the boundary
        s = distance2boundary(grid_->box(icell_),pos_,dir_);
        float ds = std::max(s,1.f)*std::numeric_limits<float>::epsilon()*100;
        s += ds;
        x = pos_ + s*dir_;
        grid_->apply_bc(x);
        while (grid_->contains(icell_,x)) {
            s += ds;
            x = pos_ + s*dir_;
            grid_->apply_bc(x);
        };
        // 2. still exiting ?
        if (!grid_->contains_with_bc(x)) {
            pos_ = x;
            path_ += s;
            prev_cellid_ = cellid_;
            cellid_ = -1;
            return BoundaryCrossing::External;
        } else { // ion is still in target. just change the cell
            // get new pos and cell
            pos_ = x;
            path_ += s;
            prev_cellid_ = cellid_;
            icell_ = grid_->pos2cell(x);
            cellid_ = grid_->cellid(icell_);
            return BoundaryCrossing::Internal;
        }
    }
}






