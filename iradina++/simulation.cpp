#include "simulation.h"
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
    dedx1(s.dedx1),
    de_strag_(s.de_strag_),
    max_fp_(s.max_fp_),
    max_impact_par_(s.max_impact_par_),
    scattering_matrix_(s.scattering_matrix_),
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
        sqrtfp_const[i] = std::sqrt(par_.flight_path_const / materials[i]->atomicDistance());

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

    /*
     * create max_fp_ tables for all ion - material
     * combinations, # =
     * (all target atoms + projectile ) x (all target materials)
     * For each combi, get a corteo dedx table
     */
    max_fp_ = ArrayNDf(natoms,nmat,nerg);
    float dEmin = 0.01f; /// @TODO: dEmin should be user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            float* p = &max_fp_(z1,im,0);
            const float* q = &dedx_(z1,im,0);
            for(dedx_index ie; ie!=ie.end(); ie++)
                p[ie] = dEmin*(*ie)/(*q++);
        }
    }

    dedx1 = ArrayNDf(nmat,nerg);
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
            float Nl0 = mat->atomicDensity()*mat->atomicDistance(); // at/nm^2
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
            case Corteo4bit:
                scattering_matrix_(z1,z2) = new xs_lab_zbl_corteo4bit;
                break;
            case Corteo6bit:
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
     * create arrays of sig(E) - "total" cross-section vs E
     * for each ion in each materials
     * # of combinations =
     * (all target atoms + projectile ) x (all materials)
     */
    max_impact_par_ = ArrayNDf(natoms, nmat, dedx_index::size);
    float Tmin = 1e-6f; // TODO: this should be user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            const material* m = materials[im];
            float* p = &max_impact_par_(z1,im,0);
            for(const atom* a : m->atoms())
            {
                int z2 = a->id();
                float x=a->X();
                for(dedx_index ie; ie!=ie.end(); ie++) {
                    float d = scattering_matrix_(z1,z2)->impactPar(*ie, Tmin);
                    p[ie] += x*d*d;
                }
            }
            for(dedx_index ie; ie!=ie.end(); ie++) p[ie] = std::sqrt(p[ie]);
        }
    }

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
                            const vector3& dir0, const float &mass_ratio, tally& t, pka_event *pka)
{
    // clone the projectile ion
    ion* j = q_.new_ion(*proj);

    // calc recoil direction from momentum conserv.
    float f1, f2;
    f1 = proj->erg() / recoil_erg;
    f2 = std::sqrt(f1);
    f1 = std::sqrt(f1+1);
    j->dir() = (f1*dir0 - f2*(proj->dir()))*mass_ratio;

    // adjust recoil atom type and energy
    j->myAtom() = target;
    // recoil energy is reduced by lattice binding energy
    j->erg() = recoil_erg - target->El();
    // add El to phonon energy
    float El = target->El();
    t(Event::Phonon,*j,&El);

    // add lattice energy recoil to pka Tdam
    if (pka) pka->Tdam() += target->El();

    // adjust recoil id and store in respective queue
    j->recoil_id()++;
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
            tion_(Event::NewRecoil,*j); // a PKA is also a recoil

            // calc LSS/NRT - QC type damage
            // based on material / average Ed, Z, M
//            const material* mat = target_->cell(pka.cellid());
//            float Tdam_LSS_ = mat->LSS_Tdam(pka.recoilE());
//            tally_.Tdam_LSS(pka.cellid()) += Tdam_LSS_;
//            tally_.Vnrt_LSS(pka.cellid()) += mat->NRT(Tdam_LSS_);

            const atom* z2 = j->myAtom();
            float dp[2];
            dp[0] = z2->LSS_Tdam(pka.recoilE());
            dp[1] = z2->NRT(dp[0]);
            tion_(Event::NRT_LSS_damage,*j,dp);

            // FullCascade or CascadesOnly
            if (par_.simulation_type != IonsOnly) {
                transport(j,tion_,&pka);
                // transport all secondary recoils
                ion* k;
                while ((k = q_.pop_recoil())!=nullptr) {
                    tion_(Event::NewRecoil,*k);
                    transport(k,tion_,&pka);
                }
                // store total Tdam & NRT vacancies
                dp[0] = pka.Tdam();
                dp[1] = z2->NRT(pka.Tdam());
                tion_(Event::NRT_damage,*j,dp);
            }
            pka_stream_.write(&pka);
        } // pka loop

        // add this ion's tally to total score
        tally_ += tion_;
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

/**
 * @brief Transport an ion through the target
 * @param i ion to transport
 * @return 0 if succesfull
 */
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
                ip > 2.f*mat->layerDistance()) {
                // TODO: register missed collision
                continue;
            }

            auto xs = scattering_matrix_(i->myAtom()->id(),z2->id());
            float T; // recoil energy
            float sintheta, costheta; // scattering angle in Lab sys
            assert(i->erg() > 0);
            xs->scatter(i->erg(), ip, T, sintheta, costheta);
            float nx, ny; // azimuthial dir
            rng.azimuth(nx,ny);

            vector3 dir0 = i->dir(); // store initial dir
            i->deflect(
                vector3(nx*sintheta, ny*sintheta, costheta)
                );
            i->erg() -= T;

            // TODO: special treatment of surface effects (sputtering etc.)

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
                    // Count replacement, energy goes to phonons {???}
                    t(Event::Replacement,*i);
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

/**
 * @brief Decide ion flight-path and impact parameter
 * @param i is the ion
 * @param m is the material the ion is travelling in
 * @param fp is the flight-path [nm]
 * @param ip is the impact parameter [nm]
 * @param sqrtfp is sqrt(fp/atomicDistance) - used for calculating straggling
 * @return
 */
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
        fp = m->atomicDistance();
        sqrtfp = 1.f;
        ip = m->meanImpactPar()*std::sqrt(rng.u01());
        break;
    case Constant:
        fp = par_.flight_path_const;
        sqrtfp = sqrtfp_const[m->id()];
        ip = m->meanImpactPar()/sqrtfp_const[m->id()]*std::sqrt(rng.u01());
        break;
    case SRIMlike:
        {
            float epsilon = i->erg() * m->meanF();
            float xsi = std::sqrt(epsilon * m->meanMinRedTransfer());
            float bmax = 1.f/(xsi + std::sqrt(xsi) + 0.125*std::pow(xsi,.1f));
            ipmax = bmax * m->meanA();
            fp = 1./(M_PI * m->atomicDensity() * ipmax * ipmax);
            sqrtfp = std::sqrt(fp/m->atomicDistance());
            ip = ipmax*std::sqrt(rng.u01());
        }
        break;
    case MendenhallWeller:
        ipmax = max_impact_par_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        if (ipmax < m->atomicDistance()) // TODO: check def of impact par
        {
            fp = 1.f/(M_PI*m->atomicDensity()*ipmax*ipmax);
            sqrtfp = std::sqrt(fp/m->atomicDistance());
            ip = std::sqrt(-std::log(rng.u01open()))*ipmax;
            if (ip > ipmax) ip = INFINITY;
        } else { // atomic spacing
            fp = m->atomicDistance();
            sqrtfp = 1.f;
            ip = m->atomicDistance()*std::sqrt(rng.u01());
        }
        break;
    case MyFFP:
        ipmax = max_impact_par_(i->myAtom()->id(),m->id(),dedx_index(i->erg()));
        if (ipmax < m->atomicDistance()) // TODO: check def of impact par
        {
            fp = 1.f/(M_PI*m->atomicDensity()*ipmax*ipmax);
            fp *= (-std::log(rng.u01open()));
            sqrtfp = std::sqrt(fp/m->atomicDistance());
            //ip = -std::log(urbg.u01open())*ipmax;
            //if (ip > ipmax) ip = INFINITY;
            ip = ipmax*std::sqrt(rng.u01lopen());
        } else { // atomic spacing
            fp = m->atomicDistance();
            sqrtfp = 1.f;
            ip = m->atomicDistance()*std::sqrt(rng.u01lopen());
        }
        break;
    default:
        fp = ip = 0.f;
    }
    return 0;
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




