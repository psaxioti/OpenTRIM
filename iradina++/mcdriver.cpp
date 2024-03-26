#include "mcdriver.h"
#include "out_file.h"
#include "elements.h"

#include <iostream>
#include <thread>
#include <chrono>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cctype>

#define CHECK_EMPTY_VEC(MatDesc,V,MatName) \
if (MatDesc.V.size()==0) throw std::invalid_argument("Empty "#V" in " + MatName);

#define CHECK_VEC_SIZE(MatDesc,V,N,MatName) \
if (MatDesc.V.size()!=N) throw std::invalid_argument("In material descriptor " + MatName + " the # of values in "#V" is not equal to those of Z." );

#define ANY_ZEROorNEG(V) \
std::any_of(V.begin(),V.end(),[](float c){ return c<=0.f; })

#define CHECK_VEC_ZEROorNEG(MatDesc,V,MatName) \
    if (ANY_ZEROorNEG(MatDesc.V)) throw std::invalid_argument("In material descriptor " + MatName + " "#V" has zero or negative values." );


mcdriver::mcdriver() :
    s_(nullptr), ips_(0.)
{

}

mcdriver::~mcdriver()
{
    if (s_) delete s_;
}

int mcdriver::saveTallys()
{
    out_file of(s_);
    std::string fname(out_opts_.OutputFileBaseName);
    fname += ".h5";
    if (of.open(fname.c_str())!=0) return -1;
    options opt;
    getOptions(opt);
    of.save(opt);
    of.close();
    return 0;
}

std::string mcdriver::outFileName(const char* type, int thread_id)
{
    std::stringstream ss;
    ss << out_opts_.OutputFileBaseName << '.' << type;
    if (thread_id) ss << thread_id;
    return ss.str();
}

void mcdriver::getOptions(options& opt) const
{
    opt.Driver = par_;
    opt.Output = out_opts_;
    opt.Simulation = s_->getParameters();
    opt.IonBeam = s_->getSource().getParameters();
    opt.Target = s_->getTarget().getDescription();
    s_->getTarget().getMaterialDescriptions(opt.materials_desc);
    opt.regions_desc = s_->getTarget().regions();
}

void mcdriver::setOptions(const options& o)
{
    if (s_) {
        delete s_;
        s_ = nullptr;
    }
    s_ = o.createSimulation();
    par_ = o.Driver;
    out_opts_ = o.Output;
    s_->init();
}

int mcdriver::exec(progress_callback cb, uint msInterval)
{
    using namespace std::chrono_literals;

    int nthreads = par_.threads;
    if (nthreads < 1) nthreads = 1;

    // TIMING
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_start);

    // create simulation clones
    std::vector< mccore* > sims(nthreads);
    sims[0] = s_;
    for(int i=1; i<nthreads; i++) sims[i] = new mccore(*s_);


    // open streams
    for(int i=0; i<nthreads; i++) {
        if (out_opts_.store_pka)
            sims[i]->open_pka_stream(outFileName("pka",i).c_str());
        if (out_opts_.store_transmitted_ions)
            sims[i]->open_exit_stream(outFileName("exit",i).c_str());
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
    unsigned int N = par_.max_no_ions;
    unsigned int Nth = N/nthreads;
    uint offset = 0;
    for(int i=0; i<nthreads; i++) {
        uint n = i==nthreads-1 ? N : Nth;
        sims[i]->setMaxIons(n);
        sims[i]->setCountOffset(offset);
        offset += n;
        N -= n;
        sims[i]->seed(myseeds[i]);
    }

    // create & start worker threads
    std::vector< std::thread* > threads;
    for(int i=0; i<nthreads; i++) {
        threads.push_back(
            new std::thread(&mccore::run, sims[i])
            );
    }

    // report progress if callback function is given
    if (cb) {
        thread_ion_count_.assign(nthreads+1,0);
        unsigned int n = 0;
        while(n < par_.max_no_ions) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msInterval));
            int i=0;
            for(; i<nthreads; i++) {
                thread_ion_count_[i] = sims[i]->ions_done();
                n += thread_ion_count_[i];
            }
            thread_ion_count_[i] = n;
            cb(*this);
        }
    }

    // wait for threads to finish...
    for(int i=0; i<nthreads; i++) threads[i]->join();

    // consolidate results
    for(int i=1; i<nthreads; i++)
        sims[0]->merge(*(sims[i]));

    if (out_opts_.store_pka) {
        std::string h5fname = outFileName("pka", 0);
        h5fname += ".h5";
        sims[0]->pka_stream().saveH5(h5fname.c_str(), "pka");
    }
    if (out_opts_.store_transmitted_ions) {
        std::string h5fname = outFileName("exit", 0);
        h5fname += ".h5";
        sims[0]->exit_stream().saveH5(h5fname.c_str(), "exit");
    }
    for(int i=0; i<nthreads; i++) sims[i]->remove_stream_files();


    // delete threads
    for(int i=0; i<nthreads; i++) {
        delete threads[i];
    }
    // delete simulation clones
    for(int i=1; i<nthreads; i++) {
        delete sims[i];
    }

    // CALC TIME/ion CLOCK_PROCESS_CPUTIME_ID
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_end); // POSIX
    double t_secs = 1. * (t_end.tv_sec - t_start.tv_sec) / nthreads;
    t_secs += 1.e-9 * (t_end.tv_nsec - t_start.tv_nsec) / nthreads;
    ips_ = s_->getTally().Nions()/t_secs;

    return 0;
}

void mcdriver::def_progress_callback(const mcdriver& d)
{
    static int i = 0;

    auto v = d.thread_ion_count();

    if (i==0) {
        for(int k=0; k<v.size()-1; k++)
            std::cout << "#" << k+1 << '\t';
        std::cout << "Total" << std::endl;
        for(int k=0; k<v.size()-1; k++)
            std::cout << "-------" << '\t';
        std::cout << "-------" << std::endl;
    }

    for(int k=0; k<v.size(); k++)
        std::cout << v[k] << '\t';
    std::cout << std::endl;

    i++;
}


int options::validate()
{
    // Driver
    if (Driver.max_no_ions <= 0)
        throw std::invalid_argument("Simulation.max_no_ions must be larger than 0.");

    if (Driver.threads <=0)
        throw std::invalid_argument("Simulation.threads must be 1 or larger.");

    if (!Driver.seeds.empty()) {
        if (Driver.seeds.size()!=Driver.threads) {
            std::stringstream ss;
            ss << "Driver.threads=" << Driver.threads << " while the # of "
               << "Driver.seeds is " << Driver.seeds.size() << "." << std::endl
               << "Either enter a # of seeds equal to the # of threads or no seeds at all.";
            throw std::invalid_argument(ss.str());
        }
    }

    // Simulation
    if (Simulation.flight_path_type==mccore::Constant &&
        Simulation.flight_path_const<=0.f)
        throw std::invalid_argument("Simulation.flight_path_type is \"Constant\" but Simulation.flight_path_const is negative.");

    // Output
    const std::string& fname = Output.OutputFileBaseName;
    if (fname.empty())
        throw std::invalid_argument("Output.OutputFileBaseName is empty.");

    if (std::any_of(fname.begin(),
                    fname.end(),
                    [](unsigned char c){ return !std::isalnum(c); }))
    {
        std::string msg = "Output.OutputFileBaseName=\"";
        msg += fname;
        msg += "\" contains non alphanumeric characters.";
        throw std::invalid_argument(msg);
    }

    // Target
    if (Target.materials.empty())
        throw std::invalid_argument("Target.materials is empty.");
    if (Target.regions.empty())
        throw std::invalid_argument("Target.regions is empty.");
    const ivector3& cell_count = Target.cell_count;
    if (std::any_of(cell_count.begin(),
                    cell_count.end(),
                    [](int c){ return c<=0; }))
    {
        std::stringstream msg;
        msg << "Target.cell_count=[";
        msg << cell_count;
        msg << "] contains zero or negative values.";
        throw std::invalid_argument(msg.str());
    }
    const vector3& cell_size = Target.cell_size;
    if (std::any_of(cell_size.begin(),
                    cell_size.end(),
                    [](float c){ return c<=0.f; }))
    {
        std::stringstream msg;
        msg << "Target.cell_size=[";
        msg << cell_size;
        msg << "] contains zero or negative values.";
        throw std::invalid_argument(msg.str());
    }

    if (Target.materials.size()!=materials_desc.size())
        throw std::invalid_argument("The # of Target.materials is not equal to "
                                    "the # of material descriptors.");
    if (Target.regions.size()!=regions_desc.size())
        throw std::invalid_argument("The # of Target.regions is not equal to "
                                    "the # of region descriptors.");

    // Material descriptors
    int k = 0;
    for(auto m : materials_desc) {
        auto mname = Target.materials[k++];
        if (m.density <= 0.f) {
            std::stringstream msg;
            msg << "Zero or negative density in material";
            msg << mname;
            throw std::invalid_argument(msg.str());
        }
        CHECK_EMPTY_VEC(m,Z,mname)
        CHECK_EMPTY_VEC(m,M,mname)
        CHECK_EMPTY_VEC(m,X,mname)
        CHECK_EMPTY_VEC(m,Ed,mname)
        CHECK_EMPTY_VEC(m,El,mname)
        CHECK_EMPTY_VEC(m,Es,mname)
        CHECK_EMPTY_VEC(m,Er,mname)

        int natoms = m.Z.size();
        if (std::any_of(m.Z.begin(),
                        m.Z.end(),
                        [](int z){ return (z < 1) || (z > elements::max_atomic_num);}))

        throw std::invalid_argument("Invalid Z number in material" + mname);

        CHECK_VEC_SIZE(m,M,natoms,mname)
        CHECK_VEC_SIZE(m,X,natoms,mname)
        CHECK_VEC_SIZE(m,Ed,natoms,mname)
        CHECK_VEC_SIZE(m,El,natoms,mname)
        CHECK_VEC_SIZE(m,Es,natoms,mname)
        CHECK_VEC_SIZE(m,Er,natoms,mname)
        CHECK_VEC_ZEROorNEG(m,M,mname)
        CHECK_VEC_ZEROorNEG(m,X,mname)
        CHECK_VEC_ZEROorNEG(m,Ed,mname)
        CHECK_VEC_ZEROorNEG(m,El,mname)
        CHECK_VEC_ZEROorNEG(m,Es,mname)
        CHECK_VEC_ZEROorNEG(m,Er,mname)
    }
    return 0;
}

mccore* options::createSimulation() const
{
    mccore* S = new mccore(Simulation);
    //S->setOutputOptions(Output);
    S->getSource().setParameters(IonBeam);

    target& T = S->getTarget();

    for(int i=0; i<materials_desc.size(); i++) {
        material* m = T.addMaterial(Target.materials[i].c_str());
        const material::material_desc_t& md = materials_desc[i];
        m->setMassDensity(md.density);
        for(int i=0; i<md.Z.size(); i++)
            m->addAtom(
                {.Z=md.Z[i],.M=md.M[i],.Ed=md.Ed[i],
                 .El=md.El[i],.Es=md.Es[i],.Er=md.Er[i]},
                md.X[i]);
    }

    grid3D& G = T.grid();
    G.setX(0, Target.cell_count.x()*Target.cell_size.x(),
           Target.cell_count.x(), Target.periodic_bc.x());
    G.setY(0, Target.cell_count.y()*Target.cell_size.y(),
           Target.cell_count.y(), Target.periodic_bc.y());
    G.setZ(0, Target.cell_count.z()*Target.cell_size.z(),
           Target.cell_count.z(), Target.periodic_bc.z());

    const std::vector<material*>& imat = T.materials();
    for(const target::region& rd : regions_desc) {
        box3D box;
        box.min() = rd.min;
        box.max() = rd.max;
        int i = materialIdx(rd.material_id);
        T.fill(box,imat[i]);
    }

    return S;
}
