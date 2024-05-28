#include "mccore.h"
#include "random_vars.h"
#include "dedx.h"
#include "event_stream.h"

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
    sqrtfp_const(s.sqrtfp_const),
    dedx_(s.dedx_),
    de_strag_(s.de_strag_),
    scattering_matrix_(s.scattering_matrix_),
    mfp_(s.mfp_), ipmax_(s.ipmax_),
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
    for(int i=0; i<materials.size(); i++)
        sqrtfp_const[i] = std::sqrt(par_.flight_path_const / materials[i]->atomicRadius());

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
    dedx_ = ArrayNDf(natoms,nmat,nerg);
    for(atom* at1 : atoms) {
        float amuRatio = elements::mostAbundantIsotope(at1->Z())/at1->M();
        int iat1 = at1->id();
        for(const material* mat : materials)
        {
            int im = mat->id();
            float* p = &dedx_(iat1,im,0);
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
                    p[i] += interp_dedx(erg,q)*w;
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
            const float* dedx = &dedx_(z1,im,0);
            const float* dedxH = &dedx1(im,0);
            float* p = &de_strag_(z1,im,0);
            float Nl0 = mat->atomicDensity()*mat->atomicRadius(); // at/nm^2
            for(const atom* z2 : mat->atoms())
                calcStraggling(dedx,dedxH,Z1,M1,z2->Z(),
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
                scattering_matrix_(z1,z2) = new xs_lab_zbl_corteo4bit;
                break;
            case Corteo6bitTable:
                scattering_matrix_(z1,z2) = new xs_lab_zbl_corteo6bit;
                break;
            case ZBL_MAGICK:
                scattering_matrix_(z1,z2) = new xs_lab_zbl_magic;
                break;
            default:
                scattering_matrix_(z1,z2) = new xs_lab_zbl_corteo4bit;
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
     *   - if not then mfp = min(mfp, Rat)
     *
     * Finally:
     *
     * - ipmax = 1/(pi*mfp*N)^(1/2)
     *
     */
    mfp_   = ArrayNDf(natoms,nmat,nerg);
    ipmax_ = ArrayNDf(natoms,nmat,nerg);
    float delta_dedx = 0.05f; /// @TODO: should be user option
    float Tmin = par_.min_recoil_energy;
    float Tmin_rel = 0.01f; /// @TODO: make it user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            const material* m = materials[im];
            const float & N = m->atomicDensity();
            const float & Rat = m->atomicRadius();
            for(dedx_index ie; ie!=ie.end(); ie++) {
                float & mfp = mfp_(z1,im,ie);
                float & ipmax = ipmax_(z1,im,ie);
                // float & dedxn = dedxn_(z1,im,ie);
                float E = *ie;
                float T0 = Tmin;

                // ensure Tmin is below Tm/2 of all target atoms
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

                // ensure mfp below 1% energy loss
                float dx = delta_dedx*E/dedx_(z1,im,ie);
                if (mfp > dx) mfp = dx;

                // ensure mfp not smaller than interatomic distance
                if (mfp < Rat) mfp = Rat;

                // this is the max impact parameter ip = (N*pi*mfp)^(-1/2)
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
    pka.setNatoms(natoms-1);

    return 0;
}

ion* mccore::new_recoil(const ion* proj, const atom *target, const float& recoil_erg,
                            const vector3& dir0, const float &sqrt_mass_ratio, tally& t, pka_event *pka)
{
    // clone the projectile ion
    ion* j = q_.new_ion(*proj);

    // calc recoil direction from momentum conserv.
    float f1, f2;
    f1 = proj->erg() / recoil_erg;
    f2 = std::sqrt(f1);
    f1 = std::sqrt(f1+1);
    j->dir() = (f1*dir0 - f2*(proj->dir()))*sqrt_mass_ratio;

    // adjust recoil atom type, energy, increase recoil generation id
    j->myAtom() = target;
    j->erg() = recoil_erg;
    j->recoil_id()++;

    // tally the recoil
    t(Event::NewRecoil,*j);

    // subtract lattice binding energy
    float El = target->El();
    j->erg() -= El;
    // add El to phonon energy tally
    t(Event::Phonon,*j,&El);

    // add lattice energy recoil to pka Tdam
    if (pka) pka->Tdam() += target->El();

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
        i->ion_id() = ++count_offset_;
        i->recoil_id() = 0;
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
            NRT(&j1,tion_,pka);
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

float mccore::doDedx(const ion *i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl)
{
    dedx_index ie(i->erg());

    float de_stopping = fp * interp_dedx(i->erg(), stopping_tbl);
    if (par_.straggling_model != NoStraggling) {

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

    return de_stopping;
}

int mccore::transport(ion* i, tally &t, pka_event *pka)
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
        //float sqrtfp, pmax;
        //float fp = flightPath(i, mat, sqrtfp, pmax);

        // decide flight path fp & impact par
        float fp, ip, sqrtfp;
        flightPath(i, mat, fp, ip, sqrtfp);
        t(Event::NewFlightPath,*i,&fp);

        // if we are inside a material, calc stopping + straggling
        if (mat) {
            float de = doDedx(i, mat, fp, sqrtfp, de_stopping_tbl, de_straggling_tbl);
            if (de > 0.f) {
                i->erg() -= de;
                t(Event::Ioniz,*i,&de);
            }
        }

        int cellid0 = i->cellid();
        i->propagate(fp); // apply any periodic boundary

        // check if ion left the target
        if (i->cellid() < 0) {

            // ion left the target - exit event
            i->cellid() = cellid0; // put the last cellid
            t(Event::IonExit,*i);
            if (exit_stream_.is_open()) {
                exit_ev.set(i,cellid0,fp);
                exit_stream_.write(&exit_ev);
            }
            if (pka) pka->Tdam() += i->erg();

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
        if (mat && std::isfinite(ip) && (i->erg() >= par_.min_energy)) {

            // select collision partner
            const atom* z2 = mat->selectAtom(rng);

            // get the cross-section and calculate scattering
            auto xs = scattering_matrix_(i->myAtom()->id(),z2->id());
            float T; // recoil energy
            float sintheta, costheta; // scattering angle in Lab sys
            xs->scatter(i->erg(), ip, T, sintheta, costheta);

            // get random azimuthal dir
            float nx, ny;
            rng.azimuth(nx,ny);
            t(Event::Scattering,*i);

            // apply to ion direction & energy
            vector3 dir0 = i->dir(); // store initial dir
            i->deflect(
                vector3(nx*sintheta, ny*sintheta, costheta)
                );
            i->erg() -= T;

            /// @TODO: special treatment of surface effects (sputtering etc.)

            if (T >= z2->Ed()) { // displacement

                // Create recoil and store in ion queue
                ion* j = new_recoil(i,z2,T,dir0,xs->sqrt_mass_ratio(),t,pka);

                /*
                     * Now check whether the projectile might replace the recoil
                     *
                     * Check if Z1==Z2 and E < Er
                     *
                     * This is different from Iradina, where it is required
                     * Z1==Z2 && M1==M2
                     */
                if ((i->myAtom()->Z() == z2->Z()) && (i->erg() < z2->Er())) {                    
                    // Replacement event, ion energy goes to Phonons
                    t(Event::Replacement,*i);
                    // energy is added to Tdam
                    if (pka) pka->Tdam() += i->erg();
                    break; // end of ion history
                }

                /*
                 * projectile and target Z not equal or E>Er =>
                 * Z2 vacancy is created
                 */
                t(Event::Vacancy,*j);
                if (pka) pka->addVac(z2->id()-1);

            } else { // E<Ed, recoil cannot be displaced
                // energy goes to phonons
                t(Event::Phonon,*i,&T);
                if (pka) pka->Tdam() += T;
            }

        }

        /* Check what happens to the projectile after possible collision: */
        if(i->erg() < par_.min_energy){ /* projectile has to stop. Store as implanted ion or recoil */
            t(Event::IonStop,*i);
            if (pka) {
                pka->addImpl(i->myAtom()->id()-1);
                pka->Tdam() += i->erg();
            }
            break; // history ends
        } /* else: continue proj transport */
    }

    // ion finished - free ion buffer
    q_.free_ion(i);
    return 0;
}

int mccore::flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp)
{
    if (!m) {  // Vacuum. TODO: change this! ion should go to next boundary {???}
        fp = 0.3f;
        ip = 0.f;
        return 0;
    }

    float d, ipmax;
    switch (par_.flight_path_type) {
    case AtomicSpacing:
        fp = m->atomicRadius();
        sqrtfp = 1.f;
        ip = m->meanImpactPar()*std::sqrt(rng.u01d_lopen());
        break;
    case Constant:
        fp = par_.flight_path_const;
        sqrtfp = sqrtfp_const[m->id()];
        ip = m->meanImpactPar()/sqrtfp_const[m->id()]*std::sqrt(rng.u01d_lopen());
        break;
    case SRIMlike:
        {
            float epsilon = i->erg() * m->meanF();
            float xsi = std::sqrt(epsilon * m->meanMinRedTransfer());
            float bmax = 1.f/(xsi + std::sqrt(xsi) + 0.125*std::pow(xsi,.1f));
            ipmax = bmax * m->meanA();
            fp = 1./(M_PI * m->atomicDensity() * ipmax * ipmax);
            sqrtfp = std::sqrt(fp/m->atomicRadius());
            ip = ipmax*std::sqrt(rng.u01d_lopen());
        }
        break;
    case MendenhallWeller:
        ipmax = ipmax_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        fp = mfp_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        if (ipmax < m->atomicRadius()) // TODO: check def of impact par
        {
            sqrtfp = std::sqrt(fp/m->atomicRadius());
            ip = ipmax*std::sqrt(-std::log(rng.u01s_open()));
            if (ip > ipmax) ip = INFINITY;
        } else { // atomic spacing
            fp = m->atomicRadius();
            sqrtfp = 1.f;
            ip = ipmax*std::sqrt(rng.u01d_lopen());
        }
        break;
    case MyFFP:
        fp = mfp_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        fp *= (-std::log(rng.u01s_open()));
        sqrtfp = std::sqrt(fp/m->atomicRadius());
        ip = ipmax_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        ip = ip*std::sqrt(rng.u01d_lopen());
        break;
    default:
        fp = ip = 0.f;
    }
    return 0;
}

void mccore::NRT(const ion* i, tally &t, const pka_event& pka)
{
    const atom* z2;
    const material* m;
    float dp[2];
    switch (par_.nrt_calculation) {
    case NRT_element:
        // NRT/LSS damage based on element
        z2 = i->myAtom();
        dp[0] = z2->LSS_Tdam(pka.recoilE());
        dp[1] = z2->NRT(dp[0]);
        t(Event::NRT_LSS_damage,*i,dp);
        // NRT damage based on element with true Tdam
        dp[0] = pka.Tdam();
        dp[1] = z2->NRT(dp[0]);
        t(Event::NRT_damage,*i,dp);
        break;
    case NRT_average:
        m = target_->cell(i->cellid());
        // NRT/LSS damage based on material / average Ed, Z, M
        dp[0] = m->LSS_Tdam(pka.recoilE());
        dp[1] = m->NRT(dp[0]);
        t(Event::NRT_LSS_damage,*i,dp);
        // NRT damage based on material with true Tdam
        dp[0] = pka.Tdam();
        dp[1] = m->NRT(dp[0]);
        t(Event::NRT_damage,*i,dp);
        break;
    default:
        break;
    }
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




