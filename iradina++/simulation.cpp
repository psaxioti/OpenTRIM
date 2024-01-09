#include "simulation.h"
#include "random_vars.h"
#include "dedx.h"
#include "xs.h"
#include "out_file.h"

#include <iostream>

simulation::simulation(const char *name, const reducedXS &x) :
    name_(name), xs_(x),
    rnd(new random_vars(urbg)),
    flight_path_type(Poisson),
    simulation_type(FullCascade),
    straggling_model(YangStraggling),
    energy_cutoff_(1.f)
{

}

simulation::~simulation()
{
    if (rnd) delete rnd;
    while (!ion_buffer_.empty()) {
        ion* i = ion_buffer_.front();
        ion_buffer_.pop();
        delete i;
    }
}

void simulation::setProjectile(int Z, float M, float E0)
{
    inventory_.setProjectile(Z,M);
    source_.setProjectile(inventory_.projectile(),E0);
}

int simulation::init() {

    inventory_.init(xs_, straggling_model);

    // calculate sqrt{l/l0} in each material for constant flight path
    auto materials = inventory_.materials();
    sqrtfp_const.resize(materials.size());
    for(int i=0; i<materials.size(); i++)
        sqrtfp_const[i] = std::sqrt(flight_path_const / materials[i]->atomicDistance());

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
    return 0;
}

void simulation::tallyIonizationEnergy(const ion* i, const float& v) {
    int ic = i->cellid();
    int ia = i->atom_->id();
    IonizationEnergyTally(ia, ic) += v;
}
void simulation::tallyPhononEnergy(const ion* i, const float& v) {
    int ic = i->cellid();
    int ia = i->atom_->id();
    double& d = PhononEnergyTally(ia, ic);
    d += v;
}
void simulation::tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed)
{
    int ic = i->cellid();
    int ie = (Tdam < 2.5f*Ed) ? 0 : 1;
    KPTally(ie,ic) += Tdam;
    KPTally( 2,ic) += nv;
}


int simulation::run(int count, const char *outfname)
{
    out_file_ = new out_file(this);
    if (out_file_->open(outfname)!=0) return -1;

    for(int k=0; k<count; k++) {
        ion* i = new_ion();
        source_.source_ion(urbg, target_, *i);
        transport(i);
        while (!ion_queue_.empty())
            transport(pop_ion());
    }

    out_file_->save();
    out_file_->close();
    delete out_file_;
    out_file_ = nullptr;

    return 0;
}
int simulation::transport(ion* i)
{
    const material* mat = target_.cell(i->icell());
    const float* de_stopping_tbl = nullptr;
    const float* de_straggling_tbl = nullptr;
    if (mat) inventory_.getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

    int loop_count = 0;

    while (i->erg > 1.f) {

        loop_count++;

        float sqrtfp, pmax;
        float fp = flightPath(i, mat, sqrtfp, pmax);

        if (mat) { // calc stopping + straggling

            int ie = dedx_index::fromValue(i->erg);
            float de_stopping = fp * de_stopping_tbl[ie];
            float de_straggling = de_straggling_tbl[ie] * rnd->normal() * sqrtfp;

            /*
             * Due to gaussian distribution, the straggling can in some cases
             * get so large that the projectile gains energy or suddenly looses
             * a huge amount of energy. Is this realistic? This might actually
             * happen. However, in the simulation, ions may have higher energy
             * than the initial energy.
             * We will therefore limit the |straggling| to |stopping|.
             * Furthermore, with hydrogen, the straggling is often so big,
             * that ions gain huge amount of energy, and the phononic system
             * would actually gain energy.
             */
            if (std::abs(de_straggling)>de_stopping) {
                if(de_straggling<0){
                    de_straggling = -de_stopping;
                } else {
                    de_straggling =  de_stopping;
                }
            }
            de_stopping += de_straggling;

            /*
             * The stopping tables have no values below minVal = 16 eV.
             * Therefore, we do simple linear downscaling of electronic
             * stopping below 16 eV.
             */
            if (i->erg < dedx_index::minVal)
                de_stopping *= i->erg/dedx_index::minVal;

            if (de_stopping != 0.f) {
                i->erg -= de_stopping;
                tallyIonizationEnergy(i,de_stopping);
            }

        }

        ivector3 icell0 = i->icell();
        i->propagate(fp); // apply any periodic boundary
        if (grid3D::isNull(i->icell())) {

            // TODO: exit event

            break; // ion left the target, history ends!
        }

        if (i->icell() != icell0) { // did the ion change cell ?

            mat = target_.cell(i->icell());
            if (mat) inventory_.getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

            /*
             * TODO: if material changed we should correct stopping.
             * However, since the flightlengths are significantly smaller
             * than the cell dimensions, this can only cause very small
             * errors and only in cases where the material changes and
             * also the stopping powers are very different!
             */
        }

        if (mat && (i->erg >= energy_cutoff_)) { // collision

            const atom* z2 = mat->selectAtom(urbg);
            float p = impactPar(i,mat,sqrtfp,pmax);

            /*
             * IRADINA
             * If impact parameter much larger than interatomic distance:
             * assume collision has missed
             * Is this valid? Yes: tests have shown, that increasing the
             * limit to 10 times the atomic distance yields practically
             * the same distribution of implanted ions and damage!
             *
             * This is not done for SRIM KP mode
             */
            if (flight_path_type != SRIMlike &&
                p > 2.f*mat->layerDistance()) {
                // TODO: register missed collision
                continue;
            }

            const scatteringXS* xs = inventory_.getScatteringXS(i->atom_,z2);
            float T; // recoil energy
            float sintheta, costheta; // scattering angle in Lab sys
            assert(i->erg > 0);
            xs->scatter(i->erg, p, T, sintheta, costheta);
            float nx, ny; // azimuthial dir
            rnd->azimuth(nx,ny);

            vector3 dir0 = i->dir; // store initial dir
            i->deflect(
                vector3(nx*sintheta, ny*sintheta, costheta)
                );
            i->erg -= T;

            if (simulation_type == FullCascade) { // We need to consider the recoil

                // TODO: special treatment of surface effects (sputtering etc.)

                if (T >= z2->Ed()) { // displacement

                    // TODO: count displacements in matrix [target_atom][cell]

                    // Create recoil and store in ion queue
                    ion* j = new_recoil(i,z2,T,dir0,xs->sqrtMassRatio());
                    // Energy of recoil has to be reduced by lattice binding energy
                    // add El to phonon energy
                    j->erg -= z2->El();
                    tallyPhononEnergy(j, z2->El());

                    /*
                     * Now check whether the projectile might replace the recoil
                     *
                     * Check if Z1==Z2 and E < Er
                     *
                     * This is different from Iradina, where it is required
                     * Z1==Z2 && M1==M2
                     */
                    if ((i->atom_->Z() == z2->Z()) && (i->erg < z2->Er())) {
                        // Count replacement, energy goes to phonons {???}
                        ReplacementTally(i->atom_->id(), i->cellid())++;
                        tallyPhononEnergy(i, i->erg);
                        // TODO: event: ion stops - replacement
                        break; // end of ion history
                    } else { /* projectile and target not the same */
                        // vacancy is created
                        // store the vacancy
                        VacancyTally(i->atom_->id(), i->cellid())++;
                    }

                } else { // E<Ed, recoil cannot be displaced
                    // energy goes to phonons
                    tallyPhononEnergy(i, T);
                }

            } else if (simulation_type == KP1) {

                if (T >= z2->Ed()) { // displacement
                    /*
                     * CROC
                     * Pure elemental solid after collisions apparently corresponding
                     * to  SRIM (page 7-28 book by ZBL)
                     */
                    float Tdam = LSS_Tdam(z2->Z(), z2->M(), T);
                    float nv = NRT(z2->Ed(), Tdam);

                    VacancyTally(0,i->cellid())++;
                    tallyKP(i, Tdam, nv, z2->Ed());

                    /* CROC
                     * outputs the ballistic energy and the electronic energy of the
                     * cascade as estimated by the above formulas
                     */
                    // store energy: Tdam->phonons, T-Tdam -> ionization
                    // store vacancies
                    tallyIonizationEnergy(i, i->erg - Tdam);
                    tallyPhononEnergy(i, Tdam);

                } else { // E<Ed, recoil cannot be displaced
                    // energy goes to phonons
                    tallyPhononEnergy(i, i->erg);
                }
            }

        }

        /* Check what happens to the projectile after possible collision: */
        if(i->erg < energy_cutoff_){ /* projectile has to stop. Store as implanted ion or recoil */
            InterstitialTally(i->atom_->id(), i->cellid())++;
            // energy goes to phonons
            tallyPhononEnergy(i, i->erg);
            // TODO: event: ion stops
            break; // history ends
        } /* else: Enough energy to advance to next collision site */

    }

    if (i->atom_->id()==0) {
        std::cout << "Ion: " << i->ion_id << std::endl;
    }

    // ion finished
    free_ion(i);
    return 0;
}

float simulation::flightPath(const ion* i, const material* m, float &sqrtfp, float& pmax)
{
    if (!m) return 0.3f; // Vacuum. TODO: change this! ion should go to next boundary {???}

    float fp;
    switch (flight_path_type) {
    case Poisson:
        rnd->poisson(fp,sqrtfp); // get a poisson distributed value u. temp = u^(1/2)
        fp = fp * m->atomicDistance();
        break;
    case AtomicSpacing:
        fp = m->atomicDistance();
        sqrtfp = 1;
        break;
    case Constant:
        fp = flight_path_const;
        sqrtfp = sqrtfp_const[m->id()];
        break;
    case SRIMlike:
        {
            float epsilon = i->erg * m->meanF();
            float xsi = std::sqrt(epsilon * m->meanMinRedTransfer());
            float bmax = 1.f/(xsi + std::sqrt(xsi) + std::pow(xsi,.1f));
            pmax = bmax * m->meanA();
            fp = 1./(M_PI * m->atomicDensity() * pmax * pmax);
        }
        break;
    }
    return fp;
}

float simulation::impactPar(const ion* i, const material* m, const float& sqrtfp, const float& pmax)
{
    float p, d;
    switch (flight_path_type) {
    case Poisson:
    case AtomicSpacing:
    case Constant:
        rnd->uniform01(d,p); // get a sqrt(u) TODO: make clear naming of rnd vars
        p *= m->meanImpactPar()/sqrtfp;
        break;
    case SRIMlike:
        rnd->uniform01(d,p);
        p *= pmax;
    break;
    }
    return p;
}

ion* simulation::new_ion(const ion *parent)
{
    ion* i;
    if (ion_buffer_.empty()) i = new ion(target_.grid());
    else { i = ion_buffer_.front(); ion_buffer_.pop(); }
    if (parent) *i = *parent;
    return i;
}

ion* simulation::new_recoil(const ion* proj, const atom *target, const float& recoil_erg,
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

void simulation::free_ion(ion* i)
{
    ion_buffer_.push(i);
}
ion* simulation::pop_ion()
{
    ion* i = nullptr;
    if (!ion_queue_.empty()) {
        i = ion_queue_.front();
        ion_queue_.pop();
    }
    return i;
}
void simulation::push_ion(ion* i)
{
    ion_queue_.push(i);
}

float simulation::LSS_Tdam(int Z, float M, float T)
{
    float k_dr = 0.1334f * std::pow (1.f*Z, 2.0f / 3.0f ) / std::sqrt( M );
    float e_d = 0.01014f * std::pow (1.f*Z, -7.0f / 3.0f) * T;
    float g_ed = 3.4008f * std::pow (e_d, 1.0f / 6.0f) + 0.40244 * std::pow (e_d, 0.75f) + e_d;
    return T / (1.0 + k_dr * g_ed);

}

float simulation::NRT(float Ed, float T)
{
    if (T<Ed) return 0.f;
    float v = 2*T/5/Ed;
    return (v < 1.f) ? 1.f : v;
}




