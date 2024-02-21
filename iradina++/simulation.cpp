#include "simulation.h"
#include "random_vars.h"
#include "dedx.h"
#include "event_stream.h"

#include <iostream>
#include <type_traits>

template<class CI, class array>
float interp1d(float x, CI i, array data)
{
    if (i==i.rbegin()) return *i; // x out of range
    float y1 = data[i], x1 = *i++;
    float y2 = data[i], x2 = *i;
    return y1 + (y2-y1)*(x-x1)/(x2-x1);
}

float calcDE(float x, float dx, float E0, const float* dedx)
{
    float E(E0);
    while (x>dx) {
        x -= dx;
        E -= dedx[dedx_index(E)]*dx;
    }
    E -= dedx[dedx_index(E)]*x;
    return E0-E;

}

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
    scattering_matrix_(S.scattering_matrix_),
    max_impact_par_(S.max_impact_par_)
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
        if (!scattering_matrix_.isNull()) {
            scatteringXSlab** xs = scattering_matrix_.data();
            for (int i=0; i<scattering_matrix_.size(); i++) delete xs[i];
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

    /*
     * create arrays of sig(E) - "total" cross-section vs E
     * for each ion in each materials
     * # of combinations =
     * (all target atoms + projectile ) x (all materials)
     */
    auto materials = target_->materials();
    int nmat = materials.size();
    max_impact_par_ = Array3Df(natoms, nmat, dedx_index::dim);
    float Tmin = 0.01f; // TODO: this should be user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            const material* m = materials[im];
            float* p = max_impact_par_[z1][im];
            for(const atom* a : m->atoms())
            {
                int z2 = a->id();
                float x=a->X();
                for(dedx_index ie; ie!=ie.end(); ie++) {
                    float d = scattering_matrix_[z1][z2]->impactPar(*ie, Tmin);
                    p[ie] += x*d*d;
                }
            }
            for(dedx_index ie; ie!=ie.end(); ie++) p[ie] = std::sqrt(p[ie]);
        }
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

    pka.setNatoms(target_->atoms().size()-1);
    if (out_opts_.store_pka)
        pka_stream_.open(outFileName("pka").c_str(),pka.size());
    if (out_opts_.store_transmitted_ions)
        exit_stream_.open(outFileName("exit").c_str(),exit_ev.size());

    for(int k=0; k<par_.max_no_ions; k++) {

        // generate ion
        ion* i = q_.new_ion();
        i->ion_id() = ++tally_.Nions() + count_offset_;
        i->recoil_id() = 0;
        source_->source_ion(urbg, *target_, *i);

        // transport the ion
        transport(i);

        // transport all PKAs
        ion* j;
        while ((j = q_.pop_pka()) != nullptr) {
            // transport PKA
            pka.init(j);
            tally_.Npkas()++;
            tally_.Nrecoils()++; // a PKA is also a recoil

            // calc LSS/NRT - QC type damage
            // based on material / average Ed, Z, M
//            const material* mat = target_->cell(pka.cellid());
//            float Tdam_LSS_ = mat->LSS_Tdam(pka.recoilE());
//            tally_.Tdam_LSS(pka.cellid()) += Tdam_LSS_;
//            tally_.Vnrt_LSS(pka.cellid()) += mat->NRT(Tdam_LSS_);

            const atom* z2 = j->myAtom();
            float Tdam_LSS_ = z2->LSS_Tdam(pka.recoilE());
            tally_.Tdam_LSS(pka.cellid()) += Tdam_LSS_;
            tally_.Vnrt_LSS(pka.cellid()) += z2->NRT(Tdam_LSS_);

            // FullCascade or CascadesOnly
            if (par_.simulation_type != IonsOnly) {
                transport(j,&pka);
                // transport all secondary recoils
                ion* k;
                while ((k = q_.pop_recoil())!=nullptr) {
                    transport(k,&pka);
                    tally_.Nrecoils()++;
                }
                // store total Tdam & NRT vacancies
                tally_.Tdam(pka.cellid()) += pka.Tdam();
                //tally_.Vnrt(pka.cellid()) += mat->NRT(pka.Tdam());
                tally_.Vnrt(pka.cellid()) += z2->NRT(pka.Tdam());
            }
            pka_stream_.write(&pka);
        }
        nion_thread_++;
    }

    pka_stream_.close();
    exit_stream_.close();

    return 0;
}

template<class _XScm, class _RNG_E>
void simulation<_XScm,  _RNG_E>::doDedx(ion* i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl)
{
    dedx_index ie(i->erg());
    float de_stopping = fp * interp1d(i->erg(), ie, stopping_tbl);
    if (par_.straggling_model != NoStraggling) {

        float de_straggling = straggling_tbl[ie] * rnd->normal() * sqrtfp;

        /* IRADINA
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
        if (std::abs(de_straggling)>de_stopping)
            de_straggling = (de_straggling<0) ? -de_stopping : de_stopping;

        de_stopping += de_straggling;
    }

    /* IRADINA
             * The stopping tables have no values below minVal = 16 eV.
             * Therefore, we do simple linear downscaling of electronic
             * stopping below 16 eV.
             */
    if (i->erg() < dedx_index::minVal)
        de_stopping *= i->erg()/dedx_index::minVal;

    if (de_stopping != 0.f) {
        i->erg() -= de_stopping;

        tally_.ionization(i->myAtom()->id(), i->cellid()) += de_stopping;
    }
}

/**
 * @brief Transport an ion through the target
 * @param i ion to transport
 * @return 0 if succesfull
 */
template<class _XScm, class _RNG_E>
int simulation<_XScm,  _RNG_E>::transport(ion* i, pka_event *pka)
{
    // get the material at the ion's position
    const material* mat = target_->cell(i->cellid());
    // get the corresponding dEdx & straggling table for ion/material combination
    const float* de_stopping_tbl = nullptr;
    const float* de_straggling_tbl = nullptr;
    if (mat)
        getDEtables(i->myAtom(), mat, de_stopping_tbl, de_straggling_tbl);

    while (1) {

        // decide flight path fp
        // sqrtfp, pmax variables depend on fp algorithm
        float sqrtfp, pmax;
        float fp = flightPath(i, mat, sqrtfp, pmax);

        // if we are inside a material, calc stopping + straggling
        if (mat) doDedx(i, mat, fp, sqrtfp, de_stopping_tbl, de_straggling_tbl);

        int cellid0 = i->cellid();
        i->propagate(fp); // apply any periodic boundary

        // check if ion left the target
        if (i->cellid() < 0) {

            // ion left the target - exit event
            tally_.Nexit()++;
            tally_.exits(i->myAtom()->id(), cellid0)++;
            if (out_opts_.store_transmitted_ions) {
                exit_ev.set(i,cellid0,fp);
                exit_stream_.write(&exit_ev);
            }

            break; // history ends!
        }

        // check if the ion changed cell
        if (i->cellid() != cellid0) {

            // get new material and dEdx tables
            mat = target_->cell(i->cellid());
            if (mat) getDEtables(i->myAtom(), mat, de_stopping_tbl, de_straggling_tbl);

            /* IRADINA
             * TODO: if material changed we should correct stopping.
             * However, since the flightlengths are significantly smaller
             * than the cell dimensions, this can only cause very small
             * errors and only in cases where the material changes and
             * also the stopping powers are very different!
             */
        }

        // if mat!=0 (not vacuum) and ion E>min. energy threshold -> collision
        if (mat && (i->erg() >= par_.min_energy)) {

            // select collision partner & calc p
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

            auto xs = scattering_matrix_[i->myAtom()->id()][z2->id()];
            float T; // recoil energy
            float sintheta, costheta; // scattering angle in Lab sys
            assert(i->erg() > 0);
            xs->scatter(i->erg(), p, T, sintheta, costheta);
            float nx, ny; // azimuthial dir
            rnd->azimuth(nx,ny);

            vector3 dir0 = i->dir(); // store initial dir
            i->deflect(
                vector3(nx*sintheta, ny*sintheta, costheta)
                );
            i->erg() -= T;

            // TODO: special treatment of surface effects (sputtering etc.)

            if (T >= z2->Ed()) { // displacement

                // Create recoil and store in ion queue
                new_recoil(i,z2,T,dir0,xs->sqrtMassRatio(),pka);

                /*
                     * Now check whether the projectile might replace the recoil
                     *
                     * Check if Z1==Z2 and E < Er
                     *
                     * This is different from Iradina, where it is required
                     * Z1==Z2 && M1==M2
                     */
                if ((i->myAtom()->Z() == z2->Z()) && (i->erg() < z2->Er())) {
                    // Count replacement, energy goes to phonons {???}
                    tally_.replacements(i->myAtom()->id(), i->cellid())++;
                    tally_.Nrepl()++;
                    tally_.phonons(i->myAtom()->id(), i->cellid()) += i->erg();
                    if (pka) pka->Tdam() += i->erg();
                    break; // end of ion history
                }
                /*
                 * projectile and target Z not equal or E>Er =>
                 * Z2 vacancy is created
                 */
                tally_.vacancies(z2->id(), i->cellid())++;
                tally_.Nvac()++;
                if (pka) pka->addVac(z2->id()-1);

            } else { // E<Ed, recoil cannot be displaced
                // energy goes to phonons
                tally_.phonons(i->myAtom()->id(), i->cellid()) += T;
                if (pka) pka->Tdam() += T;
            }

        }

        /* Check what happens to the projectile after possible collision: */
        if(i->erg() < par_.min_energy){ /* projectile has to stop. Store as implanted ion or recoil */
            tally_.implantations(i->myAtom()->id(), i->cellid())++;
            tally_.Nimpl()++;
            // energy goes to phonons
            tally_.phonons(i->myAtom()->id(), i->cellid()) += i->erg();
            if (pka) {
                pka->addImpl(i->myAtom()->id()-1);
                pka->Tdam() += i->erg();
            }
            // TODO: event: ion stops
            break; // history ends
        } /* else: continue proj transport */
    }

    // ion finished - free ion buffer
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
        float epsilon = i->erg() * m->meanF();
        float xsi = std::sqrt(epsilon * m->meanMinRedTransfer());
        float bmax = 1.f/(xsi + std::sqrt(xsi) + 0.125*std::pow(xsi,.1f));
        pmax = bmax * m->meanA();
        fp = 1./(M_PI * m->atomicDensity() * pmax * pmax);
    }
    case MendenhallWeller:
        pmax = max_impact_par_[i->myAtom()->id()][m->id()][dedx_index(i->erg())];
        if (pmax < m->meanImpactPar()) // TODO: check def of impact par
        {
            static const float C0 = 1.f - std::exp(-1.f);
            fp = M_PI*m->atomicDensity()*pmax*pmax;
            fp = 1.f/fp;
            pmax *= -std::log(1.f - urbg.u01open()*C0);
        } else {
            fp = m->atomicDistance();
            pmax = m->meanImpactPar()*std::sqrt(urbg.u01lopen());
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
    case MendenhallWeller:
        p = pmax; // TODO: fix flight path algorithm
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



