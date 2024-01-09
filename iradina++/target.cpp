#include "target.h"
#include "xs.h"
#include "dedx.h"
#include "elements.h"

#include <sstream>
#include <stdexcept>

void calcStraggling(const float* dedx, const float* dedx1, int Z1, const float& M1, int Z2,
                    const float& Density2, StragglingModel model, float* strag);

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

void material::init(const reducedXS& xs)
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
    meanA_ = xs.screeningLength(ionZ, meanZ_); // nm
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
    int natoms = atoms_.size();
    for(int z1 = 0; z1<natoms; z1++)
        for(int z2 = 1; z2<natoms; z2++)
        {
            scatteringXS* s = scattering_matrix_[z1][z2];
            if (s) delete s;
        }
    for(material* m : materials_) delete m;
    for(atom* a : atoms_) delete a;
}

void inventory::setProjectile(int Z, float M)
{
    atom* i = atoms_[0];
    i->Z_ = Z;
    i->M_ = M;
}

material* inventory::addMaterial(const char* name, const float& density)
{
    materials_.push_back(new material(this, name, density, materials_.size()));
    return materials_.back();
}

void inventory::init(const reducedXS &xs, StragglingModel smodel)
{
    for(material* m : materials_) m->init(xs);

    /*
     * create a scattering matrix for all ion compinations
     * # of combinations =
     * (all target atoms + projectile ) x (all target atoms)
     */
    int natoms = atoms_.size();
    //scattering_matrix_.resize(natoms * (natoms-1));
    scattering_matrix_ = Array2D<scatteringXS*>(natoms, natoms);
    for(int z1 = 0; z1<natoms; z1++)
        for(int z2 = 1; z2<natoms; z2++)
        {
            scatteringXS* sc = new scatteringXS(xs);
            sc->init(atoms_[z1]->Z(), atoms_[z1]->M(),
                     atoms_[z2]->Z(), atoms_[z2]->M());
            scattering_matrix_[z1][z2] = sc;
        }

    /*
     * create dedx tables for all ion - material
     * combinations, # =
     * (all target atoms + projectile ) x (all target materials)
     * For each combi, get a corteo dedx table
     */
    int nmat = materials_.size();
    int nerg = dedx_index::dim;
    dedx_ = Array3Df(natoms,nmat,nerg);
    for(atom* at1 : atoms_) {
        float amuRatio = elements::mostAbundantIsotope(at1->Z())/at1->M();
        int iat1 = at1->id();
        for(const material* mat : materials_)
        {
            int im = mat->id();
            float* p = dedx_[iat1][im];
            for(atom* at2 : mat->atoms_)
            {
                /*
                 * SRIM dedx tables are stored in eV / 10^15 (at/cm^2)
                 * They are multiplied by the atomicDensity [at/nm^3]
                 * factor 0.1 needed so that resulting dedx
                 * unit is eV/nm
                 */
                const float* q = dedx(at1->Z(), at2->Z());
                float w = (at2->X()) * (mat->atomicDensity()) * 0.1;
                for(dedx_index i; i!=i.end(); i++) {
                    float erg = *i;
                    dedx_index j = dedx_index::fromValue(erg * amuRatio);
                    p[j] += q[j]*w;
                }


                /* The compound correction needs to be added here !!! */
                /* following copied from corteo:
                   compound correction according to Zeigler & Manoyan NIMB35(1998)215, Eq.(16) (error in Eq. 14)
                   if(compoundCorr!=1.0f)
                   for(k=0; k<DIMD; k++) {
                   f = d2f( 1.0/(1.0+exp( 1.48*( sqrt(2.*Dval(k)/projectileMass/25e3) -7.0) )) );
                   spp[k]*=f*(compoundCorr-1.0f)+1.0f;
                   }
                */
            }
        }
    }

    dedx1 = Array2Df(nmat,nerg);
    for(const material* mat : materials_)
    {
        int im = mat->id();
        float* p1 = dedx1[im];
        for(atom* at2 : mat->atoms_)
        {
            /*
                 * SRIM dedx tables are stored in eV / 10^15 (at/cm^2)
                 * They are multiplied by the atomicDensity [at/nm^3]
                 * factor 0.1 needed so that resulting dedx
                 * unit is eV/nm
                 */
            const float* q1 = dedx(1, at2->Z());
            float w = (at2->X()) * (mat->atomicDensity()) * 0.1;
            for(dedx_index i; i!=i.end(); i++) {
                p1[i] += q1[i]*w;
            }


            /* The compound correction needs to be added here !!! */
            /* following copied from corteo:
                   compound correction according to Zeigler & Manoyan NIMB35(1998)215, Eq.(16) (error in Eq. 14)
                   if(compoundCorr!=1.0f)
                   for(k=0; k<DIMD; k++) {
                   f = d2f( 1.0/(1.0+exp( 1.48*( sqrt(2.*Dval(k)/projectileMass/25e3) -7.0) )) );
                   spp[k]*=f*(compoundCorr-1.0f)+1.0f;
                   }
                */
        }
    }

    /*
     * create straggling tables for all ion - material
     * combinations, # =
     * (all target atoms + projectile ) x (all target materials)
     * For each combi, get a corteo table
     *
     *
     */
    de_strag_ = Array3Df(natoms,nmat,nerg);
    for(int z1 = 0; z1<natoms; z1++) {
        int Z1 = atoms_[z1]->Z();
        float M1 = atoms_[z1]->M();
        for(const material* mat : materials_)
        {
            int im = mat->id();
            const float* dedx = dedx_[z1][im];
            const float* dedxH = dedx1[im];
            float* p = de_strag_[z1][im];
            float Nl0 = mat->atomicDensity()*mat->atomicDistance(); // at/nm^2
            for(const atom* z2 : mat->atoms())
                calcStraggling(dedx,dedxH,Z1,M1,z2->Z(),
                               Nl0*z2->X(),
                               smodel,p);
            for(dedx_index ie; ie!=ie.end(); ie++)
                p[ie] = std::sqrt(p[ie]);
        }
    }
}

const scatteringXS* inventory::getScatteringXS(const atom* z1, const atom* z2) const
{
    assert(z2->id()>0);
    return scattering_matrix_[z1->id()][z2->id()];
}

/**
 * @brief Get stopping & straggling energy change
 *
 * Returns the energy change from electronic stopping & straggling
 * of an ion Z1 traveling in material m with energy erg.
 *
 * The values are returned [eV/nm] (stopping) and [eV/nm^(1/2)] from precalculated tables
 *
 * @param z1 Pointer to atom struct describing Z1
 * @param m Pointer to material struct
 * @param erg The kinetic energy of Z1
 * @param dedx (out) Stoping dE/dx [eV/nm]
 * @param de_stragg (out) Straggling dE [eV/nm^(1/2)]
 * @return 0 on succes
 */
int inventory::getDE(const atom* z1, const material* m, const float& erg,
                       float &dedx, float &de_stragg) const
{
    int ia = z1->id();
    int im = m->id();
    int ie = dedx_index::fromValue(erg);

    dedx = dedx_[ia][im][ie];
    de_stragg = de_strag_[ia][im][ie];
    return 0;
}

int inventory::getDEtables(const atom* z1, const material* m,
                     const float* &dedx, const float* &de_stragg) const
{
    // float amuRatio = elements::mostAbundantIsotope(z1->Z())/z1->M();
    int ia = z1->id();
    int im = m->id();
    dedx = dedx_[ia][im];
    de_stragg = de_strag_[ia][im];
    return 0;
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




