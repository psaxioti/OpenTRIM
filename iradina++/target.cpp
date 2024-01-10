#include "target.h"
#include "xs.h"
#include "dedx.h"

#include <sstream>
#include <stdexcept>



atom::atom(class inventory* i, class material* m, int id) :
    id_(id), mat_(m), inv_(i)
{}

material::material(inventory *i, const char* name, const float& density, int id) :
    id_(id), name_(name), massDensity_(density), inv_(i)
{}

atom* material::addAtom(int Z, float M, float X,
                         float Ed, float El, float Es, float Er)
{
    atom* a = new atom(inv_, this, inv_->atoms_.size());
    inv_->atoms_.push_back(a);
    this->atoms_.push_back(a);
    a->Z_ = Z;
    a->M_ = M;
    a->X_ = X;
    a->Ed_ = Ed;
    a->El_ = El;
    a->Es_ = Es;
    a->Er_ = Er;
    return a;
}

void material::init()
{
    cumX_.resize(atoms_.size());
    // prepare concentrations
    float tot = 0;
    for(int i=0; i<atoms_.size(); i++) tot += atoms_[i]->X_;
    meanZ_ = 0;
    meanM_ = 0;
    for(int i=0; i<atoms_.size(); i++) {
        float & x = atoms_[i]->X_;
        x /= tot;
        meanZ_ += x*atoms_[i]->Z_;
        meanM_ += x*atoms_[i]->M_;
    }

    // make cumulative (for random sel)
    cumX_[0] = atoms_[0]->X_;
    for(int i=1; i<atoms_.size(); i++)
        cumX_[i] += atoms_[i]->X_ + cumX_[i-1];

    // TODO
    /*
     * For testing purposes and comparisons for SRIM,
     * we calculate the mean screening length and
     * energy reduction factor (as in ZBL85)
     */
    int ionZ = inv_->projectile()->Z();
    float ionM = inv_->projectile()->M();
    meanA_ = screeningZBL::screeningLength(ionZ, meanZ_); // nm
    meanF_ = meanA_ * meanM_ / ( ionZ * meanZ_ * (ionM + meanM_) * E2 );
    meanMinRedTransfer_ = dedx_index::minVal * meanF_ * (ionM + meanM_)*(ionM + meanM_) / (4*ionM*meanM_) ;

    //static const float AvogadroNum = 6.02214076e23;
    atomicDensityNM_ = 6.02214076e2f * massDensity_ / meanM_; // at / nm^3
    atomicDensity_ = atomicDensityNM_*1e21f; // at/cm^3
    atomicDistance_ = 1.0/std::pow(4.0*M_PI*atomicDensityNM_/3, 1.0/3); /* for conversion from cm to nm */
    layerDistance_ = 1.0/std::pow(atomicDensityNM_, 1.0/3); /* for conversion from cm to nm */
    sqrtAtomicDistance_ = std::sqrt(atomicDistance_);
    sqrtRecFlDensity_   = 1.0 / std::sqrt(M_PI * atomicDensityNM_); // TODO  x 1 / sqrt(flight_length_constant)

    meanImpactPar_ = 1.0 / std::sqrt(M_PI*atomicDensityNM_*atomicDistance_);

}

inventory::inventory()
{
    // create the projectile ion placehoder
    atoms_.push_back(new atom(this,0,0));
}

inventory::~inventory()
{

    for(material* m : materials_) delete m;
    for(atom* a : atoms_) delete a;
}

void inventory::setProjectile(int Z, float M)
{
    atom* i = atoms_[0];
    i->Z_ = Z;
    i->M_ = M;
}

void inventory::init() {
    for(material* m : materials_) m->init();
}

material* inventory::addMaterial(const char* name, const float& density)
{
    materials_.push_back(new material(this, name, density, materials_.size()));
    return materials_.back();
}

target::target()
{}

void target::fill(const box3D& box, const material* m)
{
    if (cells_.isNull() || cells_.size()!=grid_.ncells())
        cells_ = Array3D<const material*>(grid_.x().size()-1, grid_.y().size()-1, grid_.z().size()-1);

    box3D realbox = grid_.box().intersection(box);
    box1D bx, by, bz;
    bx.min() = box1D::VectorType(realbox.min().x());
    bx.max() = box1D::VectorType(realbox.max().x());
    by.min() = box1D::VectorType(realbox.min().y());
    by.max() = box1D::VectorType(realbox.max().y());
    bz.min() = box1D::VectorType(realbox.min().z());
    bz.max() = box1D::VectorType(realbox.max().z());
    ibox1D rx = grid_.x().range(bx);
    ibox1D ry = grid_.y().range(by);
    ibox1D rz = grid_.z().range(bz);
    for(int i=rx.min().x(); i<rx.max().x(); i++)
        for(int j=ry.min().x(); j<=ry.max().x(); j++)
            for(int k=rz.min().x(); k<=rz.max().x(); k++)
                cells_[i][j][k] = m;
}




