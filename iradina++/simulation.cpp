#include "simulation.h"
#include "random_vars.h"
#include "dedx.h"
#include "out_file.h"

#include <iostream>

template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::simulation(const char *name) :
    simulation_base(name), rnd(nullptr)
{}

template<class _XScm, class _RNG_E>
simulation<_XScm,  _RNG_E>::~simulation()
{
    int natoms = inventory_.atoms().size();
    for(int z1 = 0; z1<natoms; z1++)
        for(int z2 = 1; z2<natoms; z2++)
        {
            scatteringXSlab* s = scattering_matrix_[z1][z2];
            if (s) delete s;
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
    auto atoms = inventory_.atoms();
    int natoms = atoms.size();
    scattering_matrix_ = Array2D<scatteringXSlab*>(natoms, natoms);
    for(int z1 = 0; z1<natoms; z1++)
        for(int z2 = 1; z2<natoms; z2++)
        {
            scatteringXSlab* sc = new scatteringXSlab;
            sc->init(atoms[z1]->Z(), atoms[z1]->M(),
                     atoms[z2]->Z(), atoms[z2]->M());
            scattering_matrix_[z1][z2] = sc;
        }

    /*
     * create random variables object
     */
    if (rnd) delete rnd;
    rnd = nullptr;
    switch (random_var_type) {
    case Sampled:
        rnd = new random_vars< _RNG_E >(urbg);
        break;
    case Tabulated:
        rnd = new random_vars_tbl< _RNG_E >(urbg);
        break;

    }

    return 0;
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::run(int count, const char *outfname)
{
    out_file_ = new out_file(this);
    if (out_file_->open(outfname)!=0) return -1;

    // TIMING
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    for(int k=0; k<count; k++) {

        // generate ion
        ion* i = new_ion();
        source_.source_ion(urbg, target_, *i);

        // transport the ion
        transport(i);

        // transport all ion recoils
        while (!ion_queue_.empty())
            transport(pop_ion());

        // count history
        ion_histories_++;

        if (ion_histories_ % 100 ==0) {
            std::cout << "Ion: " << ion_histories_ << std::endl;
        }
    }

    // CALC TIME/ion
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    ms_per_ion_ = (end.tv_sec - start.tv_sec) * 1.e3 / count;
    ms_per_ion_ += 1.e-6*(end.tv_nsec - start.tv_nsec) / count;

    out_file_->save();
    out_file_->close();
    delete out_file_;
    out_file_ = nullptr;

    return 0;
}

template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::transport(ion* i)
{
    const material* mat = target_.cell(i->icell());
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
            if (straggling_model != NoStraggling) {

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
            if (mat) getDEtables(i->atom_, mat, de_stopping_tbl, de_straggling_tbl);

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

            const scatteringXSlab* xs = scattering_matrix_[i->atom_->id()][z2->id()];
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

    // ion finished
    free_ion(i);
    return 0;
}

template<class _XScm, class _RNG_E>
float simulation<_XScm,  _RNG_E>::flightPath(const ion* i, const material* m, float &sqrtfp, float& pmax)
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

template<class _XScm, class _RNG_E>
float simulation<_XScm, _RNG_E>::impactPar(const ion* i, const material* m, const float& sqrtfp, const float& pmax)
{
    float p, d;
    switch (flight_path_type) {
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
template class simulation< XS_corteo6bit,  std::mt19937 >;

template class simulation< XS_zbl_magic,   std::minstd_rand >;
template class simulation< XS_corteo4bit,  std::minstd_rand >;
template class simulation< XS_corteo6bit,  std::minstd_rand >;




