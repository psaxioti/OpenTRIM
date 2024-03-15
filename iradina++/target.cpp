#include "target.h"
#include "xs.h"
#include "dedx.h"

#include <sstream>
#include <stdexcept>



atom::atom(target *t, class material* m, int id) :
    id_(id), mat_(m), target_(t)
{}

float atom::LSS_Tdam(float recoilE) const
{
    float x = lss_Efact_ * recoilE;
    x = x + 3.4008f * std::pow (x, 1.0f / 6.0f) + 0.40244 * std::pow (x, 0.75f);
    return recoilE / (1.f + lss_Kd_*x);
}
float atom::NRT(float Tdam) const
{
    if (Tdam < p_.Ed) return 0.f;
    float v = Tdam/nrt_L_;
    return (v < 1.f) ? 1.f : v;
}

material::material(target *t, const char* name, int id) :
    id_(id), name_(name), target_(t)
{}

atom* material::addAtom(const atom::parameters& p, float x)
{
    atom* a = new atom(target_, this, target_->atoms_.size());
    target_->atoms_.push_back(a);
    this->atoms_.push_back(a);
    a->p_ = p;
    a->X_ = x;
    X_.push_back(x);
    return a;
}

void material::init()
{
    cumX_.resize(atoms_.size());
    // prepare concentrations
    float tot = 0;
    for(int i=0; i<atoms_.size(); i++) tot += X_[i];
    meanZ_ = 0;
    meanM_ = 0;
    for(int i=0; i<atoms_.size(); i++) {
        X_[i] /= tot;
        atoms_[i]->X_ /= tot;
        meanZ_ += X_[i]*atoms_[i]->Z();
        meanM_ += X_[i]*atoms_[i]->M();
    }

    // make cumulative (for random sel)
    cumX_[0] = X_[0];
    for(int i=1; i<atoms_.size(); i++)
        cumX_[i] += X_[i] + cumX_[i-1];

    /*
     * For testing purposes and comparisons for SRIM,
     * we calculate the mean screening length and
     * energy reduction factor (as in ZBL85)
     */
    int ionZ = target_->projectile()->Z();
    float ionM = target_->projectile()->M();
    meanA_ = screening_function<ScreeningZBL>::screeningLength(ionZ, meanZ_); // nm
    meanF_ = meanA_ * meanM_ / ( ionZ * meanZ_ * (ionM + meanM_) * E2 );
    meanMinRedTransfer_ = dedx_index::minVal * meanF_ * (ionM + meanM_)*(ionM + meanM_) / (4*ionM*meanM_) ;

    static const float AvogadroNum = 6.02214076e2f; // note the nm^3
    if (atomicDensityNM_ <= 0) {
        atomicDensityNM_ = AvogadroNum * massDensity_ / meanM_;
    } else if (massDensity_ <= 0) {
        massDensity_ = atomicDensityNM_ * meanM_ / AvogadroNum;
    }
    //static const float AvogadroNum = 6.02214076e23;
    // atomicDensity_ = atomicDensityNM_*1e21f; // at/cm^3
    atomicDistance_ = 1.0/std::pow(4.0*M_PI*atomicDensityNM_/3, 1.0/3); /* for conversion from cm to nm */
    layerDistance_ = 1.0/std::pow(atomicDensityNM_, 1.0/3); /* for conversion from cm to nm */
    sqrtAtomicDistance_ = std::sqrt(atomicDistance_);
    sqrtRecFlDensity_   = 1.0 / std::sqrt(M_PI * atomicDensityNM_); // TODO  x 1 / sqrt(flight_length_constant)

    meanImpactPar_ = 1.0 / std::sqrt(M_PI*atomicDensityNM_*atomicDistance_);

    // LSS & NRT
    lss_Kd_ = 0.1334f * std::pow (1.f*meanZ_, 2.0f / 3.0f ) / std::sqrt( meanM_ );
    lss_Efact_ = 0.01014f * std::pow (1.f*meanZ_, -7.0f / 3.0f);
    // Ghoniem and Chou JNM1988 effective Ed
    nrt_Ed_ = 0;
    for(int i=0; i<atoms_.size(); i++) nrt_Ed_ += X_[i]/atoms_[i]->Ed();
    nrt_Ed_ = 1.f/nrt_Ed_;
    nrt_L_ = 5*nrt_Ed_/2;
    // atoms
    for(int i=0; i<atoms_.size(); i++) {
        atom* a = atoms_[i];
        a->lss_Kd_ = 0.1334f * std::pow (1.f*a->Z(), 2.0f / 3.0f ) / std::sqrt( a->M() );
        a->lss_Efact_ = 0.01014f * std::pow (1.f*a->Z(), -7.0f / 3.0f);
        a->nrt_L_ = 5*a->Ed()/2;
    }

}

material::material_desc_t material::getDescription() const
{
    material_desc_t md;
    md.density = massDensity_;
    md.isMassDensity = true;
    for(const atom* a : atoms_) {
        md.Z.push_back(a->Z());
        md.M.push_back(a->M());
        md.Ed.push_back(a->Ed());
        md.El.push_back(a->El());
        md.Es.push_back(a->Es());
        md.Er.push_back(a->Er());
    }
    md.X = X_;
    return md;
}

float material::LSS_Tdam(float recoilE) const
{
    float x = lss_Efact_ * recoilE;
    x = x + 3.4008f * std::pow (x, 1.0f / 6.0f) + 0.40244 * std::pow (x, 0.75f);
    return recoilE / (1.f + lss_Kd_*x);
}
float material::NRT(float Tdam) const
{
    if (Tdam < nrt_Ed_) return 0.f;
    float v = Tdam/nrt_L_;
    return (v < 1.f) ? 1.f : v;
}


target::target()
{
    // create the projectile ion placehoder
    atoms_.push_back(new atom(this,0,0));
}

target::target(const target& t) :
    grid_(t.grid_),
    cells_(t.cells_.copy())
{
    // create projectile and copy from t
    atoms_.push_back(new atom(this,0,0));
    atoms_[0]->p_ = t.atoms_[0]->p_;

    // create all materials as in t
    if (!t.materials_.empty()) {
        for(const material* m : t.materials_) {
            material* m1 = addMaterial(m->name().c_str());
            for(const atom* a : m->atoms_)
                m1->addAtom(a->p_, a->X());
            m1->massDensity_ = m->massDensity_;
            m1->X_ = m->X_;
            m1->cumX_ = m->cumX_;
            m1->atomicDistance_ = m->atomicDistance_;
            m1->sqrtAtomicDistance_ = m->sqrtAtomicDistance_;
            m1->layerDistance_ = m->layerDistance_;
            m1->atomicDensityNM_ = m->atomicDensityNM_;
            m1->sqrtRecFlDensity_ = m->sqrtRecFlDensity_;
            m1->meanZ_ = m->meanZ_;
            m1->meanM_ = m->meanM_;
            m1->meanF_ = m->meanF_;
            m1->meanA_ = m->meanA_;
            m1->meanMinRedTransfer_ = m->meanMinRedTransfer_;
            m1->meanImpactPar_ = m->meanImpactPar_;
        }
    }
}

target::~target()
{
    for(material* m : materials_) delete m;
    for(atom* a : atoms_) delete a;
}

void target::setProjectile(int Z, float M)
{
    atom* i = atoms_[0];
    i->p_.Z = Z;
    i->p_.M = M;
}

void target::init() {
    for(material* m : materials_) m->init();
}

material* target::addMaterial(const char* name)
{
    materials_.push_back(new material(this, name, materials_.size()));
    return materials_.back();
}

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
    for(int i=rx.min().x(); i<=rx.max().x(); i++)
        for(int j=ry.min().x(); j<=ry.max().x(); j++)
            for(int k=rz.min().x(); k<=rz.max().x(); k++)
                cells_[i][j][k] = m;

    region_desc_t rd;
    rd.material_id = m->name();
    rd.min = realbox.min();
    rd.max = realbox.max();
    regions_.push_back( rd );
}

target::target_desc_t target::getDescription() const
{
    target_desc_t td;

    for(const material* m : materials_)
        td.materials.push_back(m->name());
    int k = 1;
    for(const region_desc_t& rd : regions_) {
        std::string nm = "R" + std::to_string(k++);
        td.regions.push_back(nm);
    }
    td.cell_count = {(int)grid_.x().size() - 1, (int)grid_.y().size() - 1, (int)grid_.z().size() - 1};
    td.periodic_bc = {grid_.x().periodic(), grid_.y().periodic(), grid_.z().periodic()};
    td.cell_size = {grid_.x()[1] - grid_.x()[0],
                    grid_.y()[1] - grid_.y()[0],
                    grid_.z()[1] - grid_.z()[0]
    };

    return td;
}

void target::getMaterialDescriptions(std::vector< material::material_desc_t >& mds)
{
    for(const material* m : materials_) {
        mds.push_back(m->getDescription());
    }
}





