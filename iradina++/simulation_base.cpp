#include "simulation.h"
#include "dedx.h"
#include "elements.h"

#include <iostream>

void calcStraggling(const float* dedx, const float* dedx1, int Z1, const float& M1,
                    int Z2, const float& Ns,
                    simulation_base::straggling_model_t model, float* strag);

simulation_base::simulation_base(const char *name) :
    name_(name),
    simulation_type(FullCascade),
    flight_path_type(Poisson),
    straggling_model(YangStraggling),
    random_var_type(Sampled),
    energy_cutoff_(1.f)
{}

simulation_base::~simulation_base()
{
    while (!ion_buffer_.empty()) {
        ion* i = ion_buffer_.front();
        ion_buffer_.pop();
        delete i;
    }
}

void simulation_base::setProjectile(int Z, float M, float E0)
{
    inventory_.setProjectile(Z,M);
    source_.setProjectile(inventory_.projectile(),E0);
}

int simulation_base::init() {

    inventory_.init();

    // calculate sqrt{l/l0} in each material for constant flight path
    auto materials = inventory_.materials();
    sqrtfp_const.resize(materials.size());
    for(int i=0; i<materials.size(); i++)
        sqrtfp_const[i] = std::sqrt(flight_path_const / materials[i]->atomicDistance());

    /*
     * create dedx tables for all ion - material
     * combinations, # =
     * (all target atoms + projectile ) x (all target materials)
     * For each combi, get a corteo dedx table
     */
    auto atoms = inventory_.atoms();
    int natoms = atoms.size();
    int nmat = materials.size();
    int nerg = dedx_index::dim;
    dedx_ = Array3Df(natoms,nmat,nerg);
    for(atom* at1 : atoms) {
        float amuRatio = elements::mostAbundantIsotope(at1->Z())/at1->M();
        int iat1 = at1->id();
        for(const material* mat : materials)
        {
            int im = mat->id();
            float* p = dedx_[iat1][im];
            for(atom* at2 : mat->atoms())
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
    for(const material* mat : materials)
    {
        int im = mat->id();
        float* p1 = dedx1[im];
        for(atom* at2 : mat->atoms())
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
        int Z1 = atoms[z1]->Z();
        float M1 = atoms[z1]->M();
        for(const material* mat : materials)
        {
            int im = mat->id();
            const float* dedx = dedx_[z1][im];
            const float* dedxH = dedx1[im];
            float* p = de_strag_[z1][im];
            float Nl0 = mat->atomicDensity()*mat->atomicDistance(); // at/nm^2
            for(const atom* z2 : mat->atoms())
                calcStraggling(dedx,dedxH,Z1,M1,z2->Z(),
                               Nl0*z2->X(),
                               straggling_model,p);
            for(dedx_index ie; ie!=ie.end(); ie++)
                p[ie] = std::sqrt(p[ie]);
        }
    }

    /*
     * Allocate Tally Memory
     *
     * In FullCascade we keep all info (I,V,R,Ei, Ep)
     * for every separate atom type (including the ion)
     *
     * In KP mode :
     *
     *   Ei, Ep store only 1 col for total ion + recoils
     *
     *   V tally stores ion V (= PKAs)
     *
     *   Special KP Tally with 3 cols:
     *     Av. Tdam, for Ed < Tdam < 2.5 Ed
     *     Av. Tdam for Tdam > 2.5 Ed
     *
     *   I, R are not available
     */
    int nCells = target_.grid().ncells();
    int nAtoms = inventory_.atoms().size();
    if (simulation_type == FullCascade) {
        InterstitialTally = Array2Dui(nAtoms,nCells);
        ReplacementTally = Array2Dui(nAtoms,nCells);
        VacancyTally = Array2Dui(nAtoms,nCells);
        IonizationEnergyTally = Array2Dd(nAtoms,nCells);
        PhononEnergyTally = Array2Dd(nAtoms,nCells);
    } else {
        KPTally = Array2Dd(3,nCells);
        VacancyTally = Array2Dui(1,nCells);
        IonizationEnergyTally = Array2Dd(1,nCells);
        PhononEnergyTally = Array2Dd(1,nCells);
    }

    // reset counters
    ion_histories_ = 0;

    return 0;
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
int simulation_base::getDEtables(const atom* z1, const material* m,
                           const float* &dedx, const float* &de_stragg) const
{
    // float amuRatio = elements::mostAbundantIsotope(z1->Z())/z1->M();
    int ia = z1->id();
    int im = m->id();
    dedx = dedx_[ia][im];
    de_stragg = de_strag_[ia][im];
    return 0;
}

void simulation_base::tallyIonizationEnergy(const ion* i, const float& v) {
    int ic = i->cellid();
    int ia = i->atom_->id();
    IonizationEnergyTally(ia, ic) += v;
}
void simulation_base::tallyPhononEnergy(const ion* i, const float& v) {
    int ic = i->cellid();
    int ia = i->atom_->id();
    double& d = PhononEnergyTally(ia, ic);
    d += v;
}
void simulation_base::tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed)
{
    int ic = i->cellid();
    int ie = (Tdam < 2.5f*Ed) ? 0 : 1;
    KPTally(ie,ic) += Tdam;
    KPTally( 2,ic) += nv;
}



ion* simulation_base::new_ion(const ion *parent)
{
    ion* i;
    if (ion_buffer_.empty()) i = new ion(target_.grid());
    else { i = ion_buffer_.front(); ion_buffer_.pop(); }
    if (parent) *i = *parent;
    return i;
}

ion* simulation_base::new_recoil(const ion* proj, const atom *target, const float& recoil_erg,
                            const vector3& dir0, const float &mass_ratio)
{
    ion* j = new_ion(proj);
    assert(grid().contains(j->pos()));
    float f1, f2;
    f1 = proj->erg / recoil_erg;
    f2 = std::sqrt(f1);
    f1 = std::sqrt(f1+1);
    j->dir = (f1*dir0 - f2*(proj->dir))*mass_ratio;
    j->ion_id = proj->ion_id;
    j->recoil_id = proj->recoil_id + 1;
    j->atom_ = target;
    j->erg = recoil_erg;
    push_ion(j);
    return j;
}

void simulation_base::free_ion(ion* i)
{
    ion_buffer_.push(i);
}
ion* simulation_base::pop_ion()
{
    ion* i = nullptr;
    if (!ion_queue_.empty()) {
        i = ion_queue_.front();
        ion_queue_.pop();
    }
    return i;
}
void simulation_base::push_ion(ion* i)
{
    ion_queue_.push(i);
}

float simulation_base::LSS_Tdam(int Z, float M, float T)
{
    float k_dr = 0.1334f * std::pow (1.f*Z, 2.0f / 3.0f ) / std::sqrt( M );
    float e_d = 0.01014f * std::pow (1.f*Z, -7.0f / 3.0f) * T;
    float g_ed = 3.4008f * std::pow (e_d, 1.0f / 6.0f) + 0.40244 * std::pow (e_d, 0.75f) + e_d;
    return T / (1.0 + k_dr * g_ed);

}

float simulation_base::NRT(float Ed, float T)
{
    if (T<Ed) return 0.f;
    float v = 2*T/5/Ed;
    return (v < 1.f) ? 1.f : v;
}




