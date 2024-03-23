#include "simulation.h"
#include "dedx.h"
#include "elements.h"
#include "out_file.h"
#include "event_stream.h"
#include "options.h"

#include <iostream>
#include <thread>
#include <chrono>

void calcStraggling(const float* dedx, const float* dedx1, int Z1, const float& M1,
                    int Z2, const float& Ns,
                    simulation::straggling_model_t model, float* strag);

simulation::simulation() :
    source_(new ion_beam),
    target_(new target),
    ref_count_(new int(0)),
    count_offset_(0),
    nion_thread_(0)
{
}

simulation::simulation(const parameters &p) :
    par_(p),
    source_(new ion_beam),
    target_(new target),
    ref_count_(new int(0)),
    count_offset_(0),
    nion_thread_(0)
{
}

simulation::simulation(const simulation &s) :
    par_(s.par_),
    out_opts_(s.out_opts_),
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
    rng()
{
    tally_.copy(s.tally_);
}

simulation::~simulation()
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

int simulation::init() {

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
    max_fp_ = Array3Df(natoms,nmat,nerg);
    float dEmin = 0.01f; // TODO: this should be user option
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int im=0; im<materials.size(); im++)
        {
            float* p = max_fp_[z1][im];
            const float* q = dedx_[z1][im];
            for(dedx_index ie; ie!=ie.end(); ie++)
                p[ie] = dEmin*(*ie)/(*q++);
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
    scattering_matrix_ = Array2D<abstract_xs_lab*>(natoms, natoms);
    for(int z1 = 0; z1<natoms; z1++)
    {
        for(int z2 = 1; z2<natoms; z2++)
        {
            switch (par_.scattering_calculation) {
            case Corteo4bit:
                scattering_matrix_[z1][z2] = new xs_lab_zbl_corteo4bit;
                break;
            case Corteo6bit:
                scattering_matrix_[z1][z2] = new xs_lab_zbl_corteo6bit;
                break;
            case ZBL_MAGICK:
                scattering_matrix_[z1][z2] = new xs_lab_zbl_magic;
                break;
            default:
                scattering_matrix_[z1][z2] = new xs_lab_zbl_corteo4bit;
                break;
            }

            scattering_matrix_[z1][z2]->init(atoms[z1]->Z(), atoms[z1]->M(),
                                             atoms[z2]->Z(), atoms[z2]->M());
        }
    }

    /*
     * create arrays of sig(E) - "total" cross-section vs E
     * for each ion in each materials
     * # of combinations =
     * (all target atoms + projectile ) x (all materials)
     */
    max_impact_par_ = Array3Df(natoms, nmat, dedx_index::size);
    float Tmin = 1e-6f; // TODO: this should be user option
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
    int ncells = target_->grid().ncells();
    tally_.init(natoms, ncells);

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
 * @param dedx (out) pointer to stopping power dE/dx table [eV/nm]
 * @param de_stragg (out) pointer to Straggling dE table [eV/nm^(1/2)]
 * @return 0 on succes
 */
int simulation::getDEtables(const atom* z1, const material* m,
                           const float* &dedx, const float* &de_stragg) const
{
    int ia = z1->id();
    int im = m->id();
    dedx = dedx_[ia][im];
    de_stragg = de_strag_[ia][im];
    return 0;
}

int simulation::exec(progress_callback cb, uint msInterval)
{
    using namespace std::chrono_literals;

    int nthreads = par_.threads;
    if (nthreads < 1) nthreads = 1;

    // TIMING
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_start);

    // create simulations
    std::vector< simulation* > sims(nthreads);
    sims[0] = this;
    thread_id_ = 0;
    for(int i=1; i<nthreads; i++) {
        sims[i] = new simulation(*this);
        sims[i]->thread_id_ = i;
    }

    // if no seeds given, generate random seeds
    std::vector<uint> myseeds = par_.seeds;
    if (myseeds.size() < nthreads) {
        myseeds.resize(nthreads);
        std::random_device rd;
        for(int i = 0; i<nthreads; i++)
            myseeds[i] = rd();

    }
    // set # ions and seed in each thread
    unsigned int maxIons = par_.max_no_ions;
    unsigned int N = maxIons;
    unsigned int Nth = N/nthreads;
    uint offset = 0;
    for(int i=0; i<nthreads; i++) {
        uint n = i==nthreads-1 ? N : Nth;
        sims[i]->setMaxIons(n);
        sims[i]->count_offset_ = offset;
        offset += n;
        N -= n;
        sims[i]->seed(myseeds[i]);
    }
    // create worker threads
    std::vector< std::thread* > threads;
    for(int i=0; i<nthreads; i++) {
        threads.push_back(
            new std::thread(&simulation::run, sims[i])
            );
    }
    // report progress if callback function is given
    if (cb) {
        std::vector<uint> progv(nthreads+1,0);
        unsigned int n = 0;
        while(n < maxIons) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msInterval));
            int i=0;
            for(; i<nthreads; i++) {
                progv[i] = sims[i]->ions_done();
                n += progv[i];
            }
            progv[i] = n;
            cb(progv);
        }
    }
    // wait for threads to finish...
    for(int i=0; i<nthreads; i++) threads[i]->join();

    // consolidate results
    for(int i=1; i<nthreads; i++)
        sims[0]->addTally(sims[i]->getTally());

    if (out_opts_.store_pka) {
        std::string h5fname = outFileName("pka");
        h5fname += ".h5";
        std::vector<event_stream *> ev(nthreads);

        for(int i=0; i<nthreads; i++)
            ev[i] = &(sims[i]->pka_stream_);

        event_stream::merge(ev, h5fname.c_str(), "pka");
    }
    if (out_opts_.store_transmitted_ions) {
        std::string h5fname = outFileName("exit");
        h5fname += ".h5";
        std::vector<event_stream *> ev(nthreads);

        for(int i=0; i<nthreads; i++)
            ev[i] = &(sims[i]->exit_stream_);

        event_stream::merge(ev, h5fname.c_str(), "exit");
    }


    // delete threads
    for(int i=0; i<nthreads; i++) {
        delete threads[i];
    }
    // delete simulation clones
    for(int i=1; i<nthreads; i++) {
        delete sims[i];
    }
    // restore my max_no_ions
    par_.max_no_ions = maxIons;

    // CALC TIME/ion CLOCK_PROCESS_CPUTIME_ID
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_end); // POSIX
    double t_secs = 1. * (t_end.tv_sec - t_start.tv_sec) / nthreads;
    t_secs += 1.e-9 * (t_end.tv_nsec - t_start.tv_nsec) / nthreads;
    ips_ = tally_.Nions()/t_secs;

    return 0;
}

ion* simulation::new_recoil(const ion* proj, const atom *target, const float& recoil_erg,
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

int simulation::saveTallys()
{
    out_file of(this);
    std::string fname(out_opts_.OutputFileBaseName);
    fname += ".h5";
    if (of.open(fname.c_str())!=0) return -1;
    of.save();
    of.close();
    return 0;
}

std::string simulation::outFileName(const char* type)
{
    std::stringstream ss;
    ss << out_opts_.OutputFileBaseName << '.' << type;
    if (thread_id_) ss << thread_id_;
    return ss.str();
}

void simulation::getOptions(options& opt) const
{
    opt.Simulation = par_;
    opt.Output = out_opts_;
    opt.IonBeam = source_->getParameters();
    opt.Target = target_->getDescription();
    target_->getMaterialDescriptions(opt.materials_desc);
    opt.regions_desc = target_->regions();
}





