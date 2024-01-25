#include "simulation.h"
#include "random_vars.h"
#include "dedx.h"
#include "event_stream.h"

#include <iostream>
#include <type_traits>

template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::simulation(const char* t) :
    simulation_base(), rnd(nullptr)
{
    if (t) par_.title = t;

    if (std::is_same<_XScm, XS_zbl_magic>::value)
        par_.scattering_calculation =  ZBL_MAGICK;
    else if (std::is_same<_XScm, XS_corteo4bit>::value)
        par_.scattering_calculation =  Corteo4bit;
    else if (std::is_same<_XScm, XS_corteo6bit>::value)
        par_.scattering_calculation =  Corteo6bit;

    if (std::is_same<_RNG_E, std::mt19937>::value)
        par_.random_generator_type =  MersenneTwister;
    else if (std::is_same<_RNG_E, std::minstd_rand>::value)
        par_.random_generator_type =  MinStd;

}

template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::simulation(const parameters &p) :
    simulation_base(p), rnd(nullptr)
{
    /*
     * create random variables object
     */
    createRandomVars();
}

template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::simulation(const _Myt& S) :
    simulation_base(S), rnd(nullptr),
    scattering_matrix_(S.scattering_matrix_)
{
    /*
     * create random variables object
     */
    createRandomVars();
}


template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::~simulation()
{
    if (rnd) delete rnd;
    if (ref_count_.use_count() == 1) {
        scatteringXSlab** xs = scattering_matrix_.data();
        for (int i=0; i<scattering_matrix_.size(); i++) {
            delete xs[i];
        }
    }
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::init() {

    simulation_base::init();

    /*
     * create a scattering matrix for all ion compinations
     * # of combinations =
     * (all target atoms + projectile ) x (all target atoms)
     */
    auto atoms = target_->atoms();
    int natoms = atoms.size();
    scattering_matrix_ = Array2D<scatteringXSlab*>(natoms, natoms);
    for(int z1 = 0; z1<natoms; z1++)
        for(int z2 = 1; z2<natoms; z2++)
        {
            scattering_matrix_[z1][z2] = new scatteringXSlab;
            scattering_matrix_[z1][z2]->init(atoms[z1]->Z(), atoms[z1]->M(),
                     atoms[z2]->Z(), atoms[z2]->M());
        }

    createRandomVars();

    return 0;
}

template<class _XScm, class _RNG_E>
void simulation<_XScm, _RNG_E>::createRandomVars()
{
    if (rnd) delete rnd;
    switch (par_.random_var_type) {
    case Sampled:
        rnd = new random_vars< _RNG_E >(urbg);
        break;
    case Tabulated:
        rnd = new random_vars_tbl< _RNG_E >(urbg);
        break;

    }
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::run()
{
    ion i0(target_->grid());

    pka.setNatoms(target_->atoms().size()-1);
    if (out_opts_.store_pka)
        pka_stream_.open(outFileName("pka").c_str(),pka.size());

    for(int k=0; k<par_.max_no_ions; k++) {

        // generate ion
        ion* i = q_.new_ion(i0);
        i->ion_id = ++tally_.Nions() + count_offset_;
        i->recoil_id = 0;
        source_->source_ion(urbg, *target_, *i);

        // transport the ion
        transport(i);

        // transport all PKAs
        ion* j;
        while ((j = q_.pop_pka()) != nullptr) {
            // transport PKA
            pka.init(j->ion_id, j->atom_->id(), j->cellid(), j->erg);
            transport_recoil(j,pka);
            tally_.Npkas()++;
            tally_.Nrecoils()++; // a PKA is also a recoil
            // transport all secondary recoils
            ion* k;
            while ((k = q_.pop_recoil())!=nullptr) {
                transport_recoil(k,pka);
                tally_.Nrecoils()++;
            }
            pka_stream_.write(&pka);
        }

        nion_thread_++;
    }

    pka_stream_.close();

    return 0;
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::transport(ion* i)
{
    // std::cout << "ion " << i->ion_id << ", recoil " << i->recoil_id << std::endl;

    const material* mat = target_->cell(i->icell());
    const float* de_stopping_tbl = nullptr;
    const float* de_straggling_tbl = nullptr;
    if (mat)
        getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

    int loop_count = 0;

    while (i->erg > 1.f) {

        loop_count++;

        float sqrtfp, pmax;
        float fp = flightPath(i, mat, sqrtfp, pmax);

        if (mat) { // calc stopping + straggling

            int ie = dedx_index::fromValue(i->erg);
            float de_stopping = fp * de_stopping_tbl[ie];
            if (par_.straggling_model != NoStraggling) {

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
            }

            /*
             * The stopping tables have no values below minVal = 16 eV.
             * Therefore, we do simple linear downscaling of electronic
             * stopping below 16 eV.
             */
            if (i->erg < dedx_index::minVal)
                de_stopping *= i->erg/dedx_index::minVal;

            if (de_stopping != 0.f) {
                i->erg -= de_stopping;

                tally_.ionization(i->atom_->id(), i->cellid()) += de_stopping;
            }

        }

        int cellid0 = i->cellid();
        i->propagate(fp); // apply any periodic boundary
        int cellid1 = i->cellid();
        if (cellid1 < 0) {

            // TODO: exit event

            break; // ion left the target, history ends!
        }

        if (cellid1 != cellid0) { // did the ion change cell ?

            mat = target_->cell(cellid1);
            if (mat) getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

            /*
             * TODO: if material changed we should correct stopping.
             * However, since the flightlengths are significantly smaller
             * than the cell dimensions, this can only cause very small
             * errors and only in cases where the material changes and
             * also the stopping powers are very different!
             */
        }

        if (mat && (i->erg >= par_.min_energy)) { // collision

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
            if (par_.flight_path_type != SRIMlike &&
                p > 2.f*mat->layerDistance()) {
                // TODO: register missed collision
                continue;
            }

            auto xs = scattering_matrix_[i->atom_->id()][z2->id()];
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

            if (par_.simulation_type == FullCascade) { // We need to consider the recoil

                // TODO: special treatment of surface effects (sputtering etc.)

                if (T >= z2->Ed()) { // displacement

                    // TODO: count displacements in matrix [target_atom][cell]

                    // Create recoil and store in ion queue
                    ion* j = new_recoil(i,z2,T,dir0,xs->sqrtMassRatio());
                    // Energy of recoil has to be reduced by lattice binding energy
                    // add El to phonon energy
                    j->erg -= z2->El();
                    tally_.phonons(j->atom_->id(), j->cellid()) += z2->El();

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
                        tally_.replacements(i->atom_->id(), i->cellid())++;
                        tally_.Nrepl()++;
                        tally_.phonons(i->atom_->id(), i->cellid()) += i->erg;
                        //tallyPhononEnergy(i, i->erg);
                        // TODO: event: ion stops - replacement
                        break; // end of ion history
                    }
                    /* projectile and target not the same */
                    // vacancy is created
                    // store the vacancy
                    tally_.vacancies(z2->id(), i->cellid())++;
                    tally_.Nvac()++;


                } else { // E<Ed, recoil cannot be displaced
                    // energy goes to phonons
                    // tallyPhononEnergy(i, T);
                    tally_.phonons(i->atom_->id(), i->cellid()) += T;
                }

            } else if (par_.simulation_type == KP1) {

                // TODO

                if (T >= z2->Ed()) { // displacement
                    /*
                     * CROC
                     * Pure elemental solid after collisions apparently corresponding
                     * to  SRIM (page 7-28 book by ZBL)
                     */
                    float Tdam = LSS_Tdam(z2->Z(), z2->M(), T);
                    float nv = NRT(z2->Ed(), Tdam);

                    //VacancyTally(0,i->cellid())++;
                    //tallyKP(i, Tdam, nv, z2->Ed());

                    /* CROC
                     * outputs the ballistic energy and the electronic energy of the
                     * cascade as estimated by the above formulas
                     */
                    // store energy: Tdam->phonons, T-Tdam -> ionization
                    // store vacancies
                    //tallyIonizationEnergy(i, i->erg - Tdam);
                    //tallyPhononEnergy(i, Tdam);

                } else { // E<Ed, recoil cannot be displaced
                    // energy goes to phonons
                    //tallyPhononEnergy(i, i->erg);
                }
            }

        }

        /* Check what happens to the projectile after possible collision: */
        if(i->erg < par_.min_energy){ /* projectile has to stop. Store as implanted ion or recoil */
            tally_.implantations(i->atom_->id(), i->cellid())++;
            tally_.Nimpl()++;
            // energy goes to phonons
            tally_.phonons(i->atom_->id(), i->cellid()) += i->erg;
            // tallyPhononEnergy(i, i->erg);
            // TODO: event: ion stops
            break; // history ends
        } /* else: Enough energy to advance to next collision site */

    }

    // ion finished
    q_.free_ion(i);
    return 0;
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::transport_recoil(ion* i, pka_event& pka)
{

    const material* mat = target_->cell(i->icell());
    const float* de_stopping_tbl = nullptr;
    const float* de_straggling_tbl = nullptr;
    if (mat)
        getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

    int loop_count = 0;

    while (i->erg > 1.f) {

        loop_count++;

        float sqrtfp, pmax;
        float fp = flightPath(i, mat, sqrtfp, pmax);

        if (mat) { // calc stopping + straggling

            int ie = dedx_index::fromValue(i->erg);
            float de_stopping = fp * de_stopping_tbl[ie];
            if (par_.straggling_model != NoStraggling) {

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
            }

            /*
             * The stopping tables have no values below minVal = 16 eV.
             * Therefore, we do simple linear downscaling of electronic
             * stopping below 16 eV.
             */
            if (i->erg < dedx_index::minVal)
                de_stopping *= i->erg/dedx_index::minVal;

            if (de_stopping != 0.f) {
                i->erg -= de_stopping;

                tally_.ionization(i->atom_->id(), i->cellid()) += de_stopping;
            }

        }

        int cellid0 = i->cellid();
        i->propagate(fp); // apply any periodic boundary
        int cellid1 = i->cellid();
        if (cellid1 < 0) {

            // TODO: exit event

            break; // ion left the target, history ends!
        }

        if (cellid1 != cellid0) { // did the ion change cell ?

            mat = target_->cell(cellid1);
            if (mat) getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

            /*
             * TODO: if material changed we should correct stopping.
             * However, since the flightlengths are significantly smaller
             * than the cell dimensions, this can only cause very small
             * errors and only in cases where the material changes and
             * also the stopping powers are very different!
             */
        }

        if (mat && (i->erg >= par_.min_energy)) { // collision

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
            if (par_.flight_path_type != SRIMlike &&
                p > 2.f*mat->layerDistance()) {
                // TODO: register missed collision
                continue;
            }

            auto xs = scattering_matrix_[i->atom_->id()][z2->id()];
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

                // TODO: special treatment of surface effects (sputtering etc.)

                if (T >= z2->Ed()) { // displacement

                    // TODO: count displacements in matrix [target_atom][cell]

                    // Create recoil and store in ion queue
                    ion* j = new_recoil(i,z2,T,dir0,xs->sqrtMassRatio());
                    // Energy of recoil has to be reduced by lattice binding energy
                    // add El to phonon energy
                    j->erg -= z2->El();
                    tally_.phonons(j->atom_->id(), j->cellid()) += z2->El();
                    pka.Tdam() += z2->El();

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
                        tally_.replacements(i->atom_->id(), i->cellid())++;
                        tally_.Nrepl()++;
                        pka.addRepl(i->atom_->id()-1);
                        tally_.phonons(i->atom_->id(), i->cellid()) += i->erg;
                        pka.Tdam() += i->erg;
                        // TODO: event: ion stops - replacement
                        break; // end of ion history
                    }
                    /* projectile and target not the same */
                    // vacancy is created
                    // store the vacancy
                    tally_.vacancies(z2->id(), i->cellid())++;
                    tally_.Nvac()++;
                    pka.addVac(z2->id()-1);

                } else { // E<Ed, recoil cannot be displaced
                    // energy goes to phonons
                    // tallyPhononEnergy(i, T);
                    tally_.phonons(i->atom_->id(), i->cellid()) += T;
                    pka.Tdam() += T;
                }


        }

        /* Check what happens to the projectile after possible collision: */
        if(i->erg < par_.min_energy){ /* projectile has to stop. Store as implanted ion or recoil */
            tally_.implantations(i->atom_->id(), i->cellid())++;
            tally_.Nimpl()++;
            pka.addImpl(i->atom_->id()-1);
            // energy goes to phonons
            tally_.phonons(i->atom_->id(), i->cellid()) += i->erg;
            pka.Tdam() += i->erg;
            // TODO: event: ion stops
            break; // history ends
        } /* else: Enough energy to advance to next collision site */

    }

    // ion finished
    q_.free_ion(i);
    return 0;
}

template<class _XScm, class _RNG_E>
float simulation<_XScm,  _RNG_E>::flightPath(const ion* i, const material* m, float &sqrtfp, float& pmax)
{
    if (!m) return 0.3f; // Vacuum. TODO: change this! ion should go to next boundary {???}

    float fp;
    switch (par_.flight_path_type) {
    case Poisson:
        rnd->poisson(fp,sqrtfp); // get a poisson distributed value u. temp = u^(1/2)
        fp = fp * m->atomicDistance();
        break;
    case AtomicSpacing:
        fp = m->atomicDistance();
        sqrtfp = 1;
        break;
    case Constant:
        fp = par_.flight_path_const;
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

template<class _XScm, class _RNG_E>
float simulation<_XScm, _RNG_E>::impactPar(const ion* i, const material* m, const float& sqrtfp, const float& pmax)
{
    float p, d;
    switch (par_.flight_path_type) {
    case Poisson:
    case AtomicSpacing:
    case Constant:
        rnd->u_sqrtu(d,p); // get a sqrt(u) TODO: make clear naming of rnd vars
        p *= m->meanImpactPar()/sqrtfp;
        break;
    case SRIMlike:
        rnd->u_sqrtu(d,p);
        p *= pmax;
        break;
    }
    return p;
}

// explicit instantiation of all variants
template class simulation< XS_zbl_magic,   std::mt19937 >;
template class simulation< XS_corteo4bit,  std::mt19937 >;
// template class simulation< XS_corteo6bit,  std::mt19937 >;

template class simulation< XS_zbl_magic,   std::minstd_rand >;
template class simulation< XS_corteo4bit,  std::minstd_rand >;
// template class simulation< XS_corteo6bit,  std::minstd_rand >;

simulation_base* simulation_base::fromParameters(const parameters& par)
{
    simulation_base* S = nullptr;
    switch (par.scattering_calculation) {
    case Corteo4bit:
        switch (par.random_generator_type) {
        case MinStd:
            S = new SimCorteo4bit_MSRAND(par);
            break;
        case MersenneTwister:
            S = new SimCorteo4bit_MT(par);
            break;
        }
        break;
    case Corteo6bit:
//        switch (par.random_generator_type) {
//        case MinStd:
//            S = new SimCorteo4bit_MSRAND(par);
//            break;
//        case MersenneTwister:
//            S = new SimCorteo4bit_MT(par);
//            break;
//        }
//        break;
    case ZBL_MAGICK:
        switch (par.random_generator_type) {
        case MinStd:
            S = new SimZBLMagic_MSRAND(par);
            break;
        case MersenneTwister:
            S = new SimZBLMagic_MT(par);
            break;
        }
        break;
    }
    return S;
}



