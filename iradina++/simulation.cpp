#include "simulation.h"
#include "random_vars.h"
#include "dedx.h"
#include "event_stream.h"

#include <iostream>
#include <type_traits>

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

int simulation::run()
{
    tally tion, tpka;
    tion.init(target_->atoms().size(), target_->grid().ncells());
    tpka.init(target_->atoms().size(), target_->grid().ncells());

    pka.setNatoms(target_->atoms().size()-1);
    if (out_opts_.store_pka)
        pka_stream_.open(outFileName("pka").c_str(),pka.size());
    if (out_opts_.store_transmitted_ions)
        exit_stream_.open(outFileName("exit").c_str(),exit_ev.size());

    for(int k=0; k<par_.max_no_ions; k++) {

        // generate ion
        ion* i = q_.new_ion();
        tion(Event::NewSourceIon,*i);
        i->ion_id() = ++count_offset_;
        i->recoil_id() = 0;
        source_->source_ion(rng, *target_, *i);

        // transport the ion
        transport(i,tion);

        // transport all PKAs
        ion* j;
        while ((j = q_.pop_pka()) != nullptr) {
            // transport PKA
            pka.init(j);
            tpka(Event::NewRecoil,*j); // a PKA is also a recoil

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
            tpka(Event::NRT_LSS_damage,*j,dp);

            // FullCascade or CascadesOnly
            if (par_.simulation_type != IonsOnly) {
                transport(j,tpka,&pka);
                // transport all secondary recoils
                ion* k;
                while ((k = q_.pop_recoil())!=nullptr) {
                    tpka(Event::NewRecoil,*k);
                    transport(k,tpka,&pka);
                }
                // store total Tdam & NRT vacancies
                dp[0] = pka.Tdam();
                dp[1] = z2->NRT(pka.Tdam());
                tpka(Event::NRT_damage,*j,dp);
            }
            pka_stream_.write(&pka);

            tion += tpka;
            tpka.clear();
        } // pka loop

        // add this ion's tally to total score
        tally_ += tion;
        tion.clear();

        // update thread counter
        nion_thread_++;

    } // ion loop

    pka_stream_.close();
    exit_stream_.close();

    return 0;
}

float simulation::doDedx(const ion *i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl)
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
int simulation::transport(ion* i, tally &t, pka_event *pka)
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

            auto xs = scattering_matrix_[i->myAtom()->id()][z2->id()];
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
int simulation::flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp)
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
        ipmax = max_impact_par_[i->myAtom()->id()][m->id()][dedx_index(i->erg())];
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
        ipmax = max_impact_par_[i->myAtom()->id()][m->id()][dedx_index(i->erg())];
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




