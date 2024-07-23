#include "mccore.h"
#include "random_vars.h"
#include "dedx.h"
#include "event_stream.h"
#include "corteo_xs.h"
#include <iostream>
#include <type_traits>

void calcStraggling(const float* dedx, const float* dedx1, int Z1, const float& M1,
                    int Z2, const float& Ns,
                    mccore::straggling_model_t model, float* strag);

mccore::mccore() :
    source_(new ion_beam),
    target_(new target),
    ref_count_(new int(0)),
    count_offset_(0),
    nion_thread_(0)
{
}

mccore::mccore(const parameters &p) :
    par_(p),
    source_(new ion_beam),
    target_(new target),
    ref_count_(new int(0)),
    count_offset_(0),
    nion_thread_(0)
{
}

mccore::mccore(const mccore &s) :
    par_(s.par_),
    source_(s.source_),
    target_(s.target_),
    ref_count_(s.ref_count_),
    nion_thread_(0),
    count_offset_(0),
    sqrtfp_const(s.sqrtfp_const), ip0(s.ip0),
    dedx_(s.dedx_),
    de_strag_(s.de_strag_),
    scattering_matrix_(s.scattering_matrix_),
    mfp_(s.mfp_), ipmax_(s.ipmax_),
    fp_max_(s.fp_max_), Tcutoff_(s.Tcutoff_),
    rng(),
    pka(s.pka)
{
    tally_.copy(s.tally_);
    dtally_.copy(s.dtally_);
    tion_.copy(s.tion_);
    tally_.clear();
    dtally_.clear();
    tion_.clear();
}

mccore::~mccore()
{
    q_.clear();
    if (ref_count_.use_count() == 1) {
        delete source_;
        delete target_;
        if (!dedx_.isNull()) {
            dedx_interp_t** p = dedx_.data();
            for (int i=0; i<dedx_.size(); i++) delete p[i];
        }
        if (!scattering_matrix_.isNull()) {
            abstract_xs_lab** xs = scattering_matrix_.data();
            for (int i=0; i<scattering_matrix_.size(); i++) delete xs[i];
        }
    }
}

int mccore::init() {

    target_->init();

    // adjust projectile codes in target & source
    target_->setProjectile(source_->ionZ(), source_->ionM());
    source_->setProjectile(target_->projectile());

    // calculate sqrt{l/l0} in each material for constant flight path
    auto materials = target_->materials();
    sqrtfp_const.resize(materials.size());
    ip0.resize(materials.size());
    for(int i=0; i<materials.size(); i++) {
        sqrtfp_const[i] = std::sqrt(par_.flight_path_const / materials[i]->atomicRadius());
        ip0[i] = materials[i]->meanImpactPar();
        if (par_.flight_path_type == Constant) ip0[i] /= sqrtfp_const[i];
    }

    /*
     * create dedx tables for all ion - material
     * combinations, # =
     * (all target atoms + projectile ) x (all target materials)
     * For each combi, get a corteo dedx table
     */
    auto atoms = target_->atoms();
    int natoms = atoms.size();
    int nmat = materials.size();
    int nerg = dedx_index::size;
    dedx_ = ArrayND< dedx_interp_t* >(natoms,nmat);
    std::vector<float> buff(nerg);
    for(atom* at1 : atoms) {
        float amuRatio = elements::mostAbundantIsotope(at1->Z())/at1->M();
        int iat1 = at1->id();
        for(const material* mat : materials)
        {
            int im = mat->id();
            //float* p = &dedx_(iat1,im,0);
            float* p = buff.data();
            for(atom* at2 : mat->atoms())
            {
                /*
                 * SRIM dedx tables are stored in eV / 10^15 (at/cm^2)
                 * They are multiplied by the atomicDensity [at/nm^3]
                 * factor 0.1 needed so that resulting dedx
                 * unit is eV/nm
                 */
                const float* q = ::dedx(at1->Z(), at2->Z());
                float w = (at2->X()) * (mat->atomicDensity()) * 0.1;
                for(dedx_index i; i<i.end(); i++) {
                    float erg = (*i) * amuRatio;
                    p[i] += w*corteo_lin_interp<float,dedx_index>(q,erg);
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
            dedx_(iat1,im) = new dedx_interp_t(p);
        }
    }

    ArrayNDf dedx1(nmat,nerg); // proton stopping (materials x energy)
    for(const material* mat : materials)
    {
        int im = mat->id();
        float* p1 = &dedx1(im,0);
        for(atom* at2 : mat->atoms())
        {
            /*
                 * SRIM dedx tables are stored in eV / 10^15 (at/cm^2)
                 * They are multiplied by the atomicDensity [at/nm^3]
                 * factor 0.1 needed so that resulting dedx
                 * unit is eV/nm
                 */
            const float* q1 = ::dedx(1, at2->Z());
            float w = (at2->X()) * (mat->atomicDensity()) * 0.1;
            for(dedx_index i; i<i.end(); i++) {
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
    de_strag_ = ArrayNDf(natoms,nmat,nerg);
    for(int z1 = 0; z1<natoms; z1++) {
        int Z1 = atoms[z1]->Z();
        float M1 = atoms[z1]->M();
        for(const material* mat : materials)
        {
            int im = mat->id();
            const dedx_interp_t* dedx = dedx_(z1,im);
            const float* dedxH = &dedx1(im,0);
            float* p = &de_strag_(z1,im,0);
            float Nl0 = mat->atomicDensity()*mat->atomicRadius(); // at/nm^2
            for(const atom* z2 : mat->atoms())
                calcStraggling(dedx->data(),
                               dedxH,
                               Z1,M1,z2->Z(),
                               Nl0*z2->X(),
                               par_.straggling_model,p);
            for(dedx_index ie; ie!=ie.end(); ie++)
                p[ie] = std::sqrt(p[ie]);
        }        
    }

    /*
     * create a scattering matrix for all ion compinations
     * # of combinations =
     * (all target atoms + projectile ) x (all target atoms)
     */
    scattering_matrix_ = ArrayND<abstract_xs_lab*>(natoms, natoms);
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int z2 = 1; z2<natoms; z2++)
        {
            switch (par_.scattering_calculation) {
            case Corteo4bitTable:
                switch (par_.screening_type) {
                case Screening::ZBL:
                    scattering_matrix_(z1,z2) = new xs_lab_zbl;
                    break;
                case Screening::LenzJensen:
                    scattering_matrix_(z1,z2) = new xs_lab_lj;
                    break;
                case Screening::KrC:
                    scattering_matrix_(z1,z2) = new xs_lab_krc;
                    break;
                case Screening::Moliere:
                    scattering_matrix_(z1,z2) = new xs_lab_moliere;
                    break;
                default:
                    scattering_matrix_(z1,z2) = new xs_lab_zbl;
                    break;
                }
                break;
            case ZBL_MAGICK:
                scattering_matrix_(z1,z2) = new xs_lab_zbl_magic;
                break;
            default:
                scattering_matrix_(z1,z2) = new xs_lab_zbl;
                break;
            }

            scattering_matrix_(z1,z2)->init(atoms[z1]->Z(), atoms[z1]->M(),
                                            atoms[z2]->Z(), atoms[z2]->M());
        }
    }

    /*
     * create tables of parameters for flight path selection, pairs of (mfp, ipmax)
     *
     * - Set a low cutoff recoil energy T0 (min_recoil_energy)
     * - T0 = min(min_recoil_energy, 0.01*Tm)
     * - Scattering events with T<T0 are not considered in the MC
     * - T0/Tm = sin(th0/2)^2, th0 : low cutoff scattering angle
     * - ipmax = ip(e,th0) : max impact parameter
     * - sig0 = pi*ipmax^2 : total xs for events T>T0
     * - mfp = 1/(N*sig0) : mean free path
     *
     * Additional conditions:
     *
     * - Electron energy loss
     *   - Δxmax for el. energy loss below 5%
     *   - mfp = min(mfp, Dxmax)
     *
     * - Allow flight path < atomic radius Rat ??
     *   - if not then mfp = max(mfp, Rat)
     *
     * - User selected max_mfp
     *   - mfp = min(mfp, max_mfp)

     * Finally:
     *
     * - ipmax = 1/(pi*mfp*N)^(1/2)
     *
     */
    mfp_   = ArrayNDf(natoms,nmat,nerg);
    ipmax_ = ArrayNDf(natoms,nmat,nerg);
    fp_max_ = ArrayNDf(natoms,nmat,nerg);
    Tcutoff_ = ArrayNDf(natoms,nmat,nerg);
    float delta_dedx = par_.max_rel_eloss;
    float Tmin = par_.min_recoil_energy;
    float mfp_ub = par_.max_mfp;
    float Tmin_rel = 0.99f; /// @TODO: make it user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            const material* m = materials[im];
            const float & N = m->atomicDensity();
            const float & Rat = m->atomicRadius();
            float mfp_lb = par_.flight_path_type == IPP && par_.allow_sub_ml_scattering ?
                               0.f : Rat;
            for(dedx_index ie; ie!=ie.end(); ie++) {
                float & mfp = mfp_(z1,im,ie);
                float & ipmax = ipmax_(z1,im,ie);
                float & fpmax = fp_max_(z1,im,ie);
                float & T0 = Tcutoff_(z1,im,ie);
                // float & dedxn = dedxn_(z1,im,ie);
                float E = *ie;
                T0 = Tmin;

                // ensure Tmin is below Tm of all target atoms
                for(const atom* a : m->atoms()) {
                    int z2 = a->id();
                    float Tm = E*scattering_matrix_(z1,z2)->gamma();
                    if (T0 > Tmin_rel*Tm) T0 = Tmin_rel*Tm;
                }

                // Calc mfp corresponding to T0, mfp = 1/(N*sig0), sig0 = pi*sum_i{ X_i * [P_i(E,T0)]^2 }
                mfp = 0.f;
                for(const atom* a : m->atoms()) {
                    int z2 = a->id();
                    float d = scattering_matrix_(z1,z2)->impactPar(E, T0);
                    mfp += a->X()*d*d;
                }
                mfp = 1/N/M_PI/mfp; // mfp = 1/(N*sig0) = 1/(N*pi*sum_i(ip_i^2))

                // ensure mfp*dEdx/E is below max_rel_eloss
                fpmax = delta_dedx*E / dedx_(z1,im)->data()[ie];
                mfp = std::min(mfp,fpmax);

                // ensure mfp not smaller than lower bound
                mfp = std::max(mfp,mfp_lb);

                // ensure mfp not larger than upper bound
                mfp = std::min(mfp,mfp_ub);

                // Find the max impact parameter ipmax = (N*pi*mfp)^(-1/2)
                ipmax = std::sqrt(1.f/M_PI/mfp/N);

                // Calc dedxn for  T<T0 = N*sum_i { X_i * Sn(E,T0) }
                // Add this to dedx
                /// @TODO: this is very slow. dedxn is very small, can be ignored
                // We need to re-calc T0 from mfp
                /// @TODO: solve ipmax^2 = sum_i { X_i * ipmax_i(e,T0) } for T0
//                int z2 = m->atoms().front()->id();
//                float s1,c1;
//                scattering_matrix_(z1,z2)->scatter(E,ipmax,T0,s1,c1);
//                dedxn = 0;
//                for(const atom* a : m->atoms()) {
//                    int z2 = a->id();
//                    dedxn += scattering_matrix_(z1,z2)->stoppingPower(E,T0) * a->X();
//                }
//                dedxn *= N;
//                if (par_.flight_path_type == MyFFP) dedx_(z1,im,ie) += dedxn;

            } // energy
        } // material
    } // Z1

    /*
     * Allocate Tally Memory
     */
    int ncells = target_->grid().ncells();
    tally_.init(natoms, ncells);
    dtally_.init(natoms, ncells);
    tion_.init(natoms, ncells);

    // prepare event buffers
    std::vector<std::string> atom_labels = target_->atom_labels();
    atom_labels.erase(atom_labels.begin());
    pka.setNatoms(natoms-1,atom_labels);

    return 0;
}

ion* mccore::new_recoil(const ion* proj, const atom *target, const float& recoil_erg,
                            const vector3& dir0, const float &sqrt_mass_ratio)
{
    // clone the projectile ion
    ion* j = q_.new_ion(*proj);
    j->init_recoil(target, recoil_erg);
    j->reset_counters();

    // calc recoil direction from momentum conserv.
    float f1, f2;
    f1 = proj->erg() / recoil_erg;
    f2 = std::sqrt(f1);
    f1 = std::sqrt(f1+1);
    j->setDir((f1*dir0 - f2*(proj->dir()))*sqrt_mass_ratio);

    // store recoil in respective queue
    if (j->recoil_id()==1) q_.push_pka(j);
    else q_.push_recoil(j);

    return j;
}

int mccore::run()
{

    for(int k=0; k<max_no_ions_; k++) {

        // generate ion
        ion* i = q_.new_ion();
        i->setId(++count_offset_);
        i->resetRecoilId();
        i->reset_counters();
        source_->source_ion(rng, *target_, *i);
        tion_(Event::NewSourceIon,*i);

        // transport the ion
        transport(i,tion_);

        // transport all PKAs
        ion* j;
        while ((j = q_.pop_pka()) != nullptr) {
            // transport PKA
            pka.init(j);
            ion j1(*j); // keep a clone ion to have initial position

            // FullCascade or CascadesOnly
            if (par_.simulation_type != IonsOnly) {
                transport(j,tion_,&pka);
                // transport all secondary recoils
                ion* k;
                while ((k = q_.pop_recoil())!=nullptr) {
                    transport(k,tion_,&pka);
                }
            }
            pka_end(&j1,tion_,pka);
            pka_stream_.write(&pka);
        } // pka loop

        // add this ion's tally to total score
        tally_ += tion_;
        dtally_.addSquared(tion_);
        tion_.clear();

        // update thread counter
        nion_thread_++;

    } // ion loop

    pka_stream_.close();
    exit_stream_.close();

    return 0;
}

float mccore::calcDedx(const ion *i, const material* m, float fp, float sqrtfp,
                     const dedx_interp_t *stopping_tbl, const float *straggling_tbl)
{
    dedx_index ie(i->erg());

    float de_stopping = fp * (*stopping_tbl)(i->erg());
    if (par_.straggling_model != NoStraggling)
    {

        float de_straggling = straggling_tbl[ie] * rng.normal() * sqrtfp;

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
        if (std::abs(de_straggling) > de_stopping)
            de_straggling = (de_straggling < 0) ? -de_stopping : de_stopping;

        de_stopping += de_straggling;
    }

    /* IRADINA
     * The stopping tables have no values below minVal = 16 eV.
     * Therefore, we do sqrt downscaling of electronic
     * stopping below 16 eV.
     */
    if (i->erg() < dedx_index::minVal)
        de_stopping *= std::sqrt(i->erg() / dedx_index::minVal);

    /*
     * This is due to some rare low-energy events (ion E<100 eV)
     * with IPP flight selection were
     * a long flight path + large straggling can cause the
     * stopping + straggling > E
     *
     * The code below changes that to almost stopping the ion,
     * E - stopping ~ 0
     */
    if (de_stopping > i->erg())
        de_stopping = 0.99*i->erg();

    return de_stopping;
}

int mccore::transport(ion* i, tally &t, pka_event *pka)
{
    // get the material at the ion's position
    const material* mat = target_->cell(i->cellid());
    // get the corresponding dEdx & straggling table for ion/material combination
    const dedx_interp_t* de_stopping_tbl = nullptr;
    const float* de_straggling_tbl = nullptr;
    std::array<const float *, 3> fp_par_tbl;
    const float* & ipmax_tbl = fp_par_tbl[0];
    const float* & mfp_tbl   = fp_par_tbl[1];
    const float* & fpmax_tbl = fp_par_tbl[2];
    float fp, ip, sqrtfp, ipmax;
    int ie;
    bool doCollision;
    if (mat) {
        getDEtables(i->myAtom(), mat, de_stopping_tbl, de_straggling_tbl);
        getMFPtables(i->myAtom(), mat, fp, sqrtfp, fp_par_tbl);
    }

    // transport loop
    while (1) {

        // Check if ion has enough energy to continue
        if(i->erg() < par_.min_energy){ 
            /* projectile has to stop. Store as implanted/interstitial atom*/
            t(Event::IonStop,*i);
            if (pka) {
                pka->addImpl(i->myAtom()->id()-1);
                pka->Tdam() += i->erg();
            }
            return 0; // history ends
        } 

        if (!mat) {  // Vacuum.
            // set a long flight path (1mm)
            // to intentionally hit a boundary
            // Then the boundary crossing algorithm
            // will take care of things
            fp = 1e6f;
            BoundaryCrossing crossing = i->propagate(fp);
            switch(crossing) {
            case BoundaryCrossing::None:
                break;
            case BoundaryCrossing::Internal:
                // register event
                t(Event::BoundaryCrossing,*i);
                i->reset_counters();
                // get new material and dEdx, mfp tables
                mat = target_->cell(i->cellid());
                if (mat) {
                    getDEtables(i->myAtom(), mat, de_stopping_tbl, de_straggling_tbl);
                    getMFPtables(i->myAtom(), mat, fp, sqrtfp, fp_par_tbl);
                }
                break;
            case BoundaryCrossing::External:
                // register ion exit event
                t(Event::IonExit,*i);
                if (exit_stream_.is_open()) {
                    exit_ev.set(i);
                    exit_stream_.write(&exit_ev);
                }
                if (pka) pka->Tdam() += i->erg();
                return 0; // ion history ends!
            }
            continue; // go to next iter
        }

        // select flight path & impact param.
        doCollision = flightPath(i,mat,fp,ip,sqrtfp,fp_par_tbl);

        // propagate ion, checking also boundary crossing
        BoundaryCrossing crossing = i->propagate(fp);

        // subtract ionization & straggling
        float de = calcDedx(i, mat, fp, sqrtfp, de_stopping_tbl, de_straggling_tbl);
        i->de_ioniz(de);

        // handle boundary
        switch(crossing) {
        case BoundaryCrossing::Internal:
            // register event
            t(Event::BoundaryCrossing,*i);
            i->reset_counters();
            // get new material and dEdx, mfp tables
            mat = target_->cell(i->cellid());
            if (mat) {
                getDEtables(i->myAtom(), mat, de_stopping_tbl, de_straggling_tbl);
                getMFPtables(i->myAtom(), mat, fp, sqrtfp, fp_par_tbl);
            }
            doCollision = false; // the collision will be in the new material
            break;
        case BoundaryCrossing::External:
            // register ion exit event
            t(Event::IonExit,*i);
            if (exit_stream_.is_open()) {
                exit_ev.set(i);
                exit_stream_.write(&exit_ev);
            }
            if (pka) pka->Tdam() += i->erg();
            return 0; // ion history ends!
        default:
            break;
        }

        // if no collision, start again
        // no collision =
        // either BoundaryCrossing::Internal OR
        // rejected by flight path algorithm
        if (!doCollision) continue;

        // Now comes the COLLISSION PART

        // select collision partner
        const atom* z2 = mat->selectAtom(rng);

        // get the cross-section and calculate scattering
        auto xs = scattering_matrix_(i->myAtom()->id(),z2->id());
        float T; // recoil energy
        float sintheta, costheta; // Lab sys scattering angle sin & cos 
        xs->scatter(i->erg(), ip, T, sintheta, costheta);
        assert(finite(T));

        // get random azimuthal dir
        float nx, ny; // nx = cos(phi), ny = sin(phi), phi: az. angle
        rng.random_azimuth_dir(nx,ny);

        // register scattering event (before changing ion data)
        // t(Event::Scattering,*i);

        // apply new ion direction & energy
        vector3 dir0 = i->dir(); // store initial dir
        i->deflect(vector3(nx*sintheta, ny*sintheta, costheta));
        i->add_coll();

        /// @TODO: special treatment of surface effects (sputtering etc.)

        // Check if recoil is displaced (T>=E_d)
        if (T >= z2->Ed()) { 

            // Create recoil (it is stored in the ion queue)
            i->de_recoil(T);
            ion* j = new_recoil(i,z2,T,dir0,xs->sqrt_mass_ratio());

            /*
             * Now check whether the projectile might replace the recoil
             *
             * Check if Z1==Z2 and E < Er
             *
             * This is different from Iradina, where it is required
             * Z1==Z2 && M1==M2
             */
            if ((i->myAtom()->Z() == z2->Z()) && (i->erg() < z2->Er())) {
                j->de_phonon(z2->El());
                // Replacement event, ion energy goes to Phonons
                t(Event::Replacement,*i,z2);
                // the energy is also added to the PKA's Tdam
                if (pka) pka->Tdam() += i->erg()+z2->El();
                return 0; // end of ion history
            }

            /*
             * Otherwise:
             * Z2 vacancy is created and ion continues
             */
            //t(Event::Vacancy,*j);
            j->de_phonon(z2->El());
            if (pka) {
                pka->Tdam() += z2->El();
                pka->addVac(z2->id()-1);
            }

        } else { // T<E_d, recoil cannot be displaced
            // energy goes to phonons
            // t(Event::Phonon,*i,&T);
            i->de_phonon(T);
            if (pka) pka->Tdam() += T;
        }

    } // main ion transport loop

    // free the ion buffer
    q_.free_ion(i);


    return 0;
}

bool mccore::flightPath(const ion* i, const material* mat, float& fp, float& ip, float& sqrtfp,
                        std::array<const float *, 3> &fp_par_tbl)
{
    bool doCollision = true;
    int ie;
    const float* & ipmax_tbl = fp_par_tbl[0];
    const float* & mfp_tbl   = fp_par_tbl[1];
    const float* & fpmax_tbl = fp_par_tbl[2];
    switch (par_.flight_path_type) {
    case AtomicSpacing:
    case Constant:
        ip = ipmax_tbl[0]*std::sqrt(rng.u01d_lopen());
        break;
    case MendenhallWeller:
        ie = dedx_index(i->erg());
        ip = ipmax_tbl[ie];
        fp = mfp_tbl[ie];
        if (ip < mat->meanImpactPar())
        {
            sqrtfp = std::sqrt(fp/mat->atomicRadius());
            ip *= std::sqrt(-std::log(rng.u01s_open()));
            doCollision = (ip <= ipmax_tbl[ie]);
        } else { // atomic spacing
            fp = mat->atomicRadius();
            sqrtfp = 1.f;
            ip = mat->meanImpactPar()*std::sqrt(rng.u01d_lopen());
        }
        break;
    case IPP:
        ie = dedx_index(i->erg());
        fp = mfp_tbl[ie]*(-std::log(rng.u01s_open()));
        doCollision = fp <= fpmax_tbl[ie];
        if (doCollision) ip = ipmax_tbl[ie]*std::sqrt(rng.u01d_lopen());
        else fp = fpmax_tbl[ie];
        sqrtfp = std::sqrt(fp/mat->atomicRadius());
        break;
    default:
        assert(false); // never get here
    }
    assert(fp>0);
    assert(finite(fp));
    return doCollision;
}

void mccore::pka_end(const ion* i, tally &t, const pka_event& pka)
{
    const atom* z2;
    const material* m;
    std::array<float, 5> dp;
    dp.fill(0.f);
    dp[0] = pka.recoilE();
    switch (par_.nrt_calculation) {
    case NRT_element:
        // NRT/LSS damage based on element
        z2 = i->myAtom();
        dp[1] = z2->LSS_Tdam(pka.recoilE());
        dp[2] = z2->NRT(dp[1]);
        // NRT damage based on element with true Tdam
        dp[3] = pka.Tdam();
        dp[4] = z2->NRT(dp[3]);

        break;
    case NRT_average:
        m = target_->cell(i->cellid());
        // NRT/LSS damage based on material / average Ed, Z, M
        dp[1] = m->LSS_Tdam(pka.recoilE());
        dp[2] = m->NRT(dp[1]);
        // NRT damage based on material with true Tdam
        dp[3] = pka.Tdam();
        dp[4] = m->NRT(dp[3]);
        break;
    default:
        break;
    }
    t(Event::CascadeComplete,*i,dp.data());
}

void mccore::merge(const mccore& other)
{
    tally_ += other.tally_;
    dtally_ += other.dtally_;
    pka_stream_.merge(other.pka_stream_);
    exit_stream_.merge(other.exit_stream_);
}

void mccore::remove_stream_files()
{
    pka_stream_.remove();
    exit_stream_.remove();
}




