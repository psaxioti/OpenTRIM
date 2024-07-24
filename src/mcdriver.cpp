#include "mcdriver.h"
#include "elements.h"

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

std::string mcdriver::outFileName(const char* type, int thread_id)
{
    std::stringstream ss;
    ss << out_opts_.OutputFileBaseName << '.' << type;
    if (thread_id) ss << thread_id;
    return ss.str();
}

std::string mcdriver::outFileName() const
{
    std::string s(out_opts_.OutputFileBaseName);
    s += ".h5";
    return s;
}

void mcdriver::getOptions(options& opt) const
{
    opt.Driver = par_;
    opt.Output = out_opts_;
    opt.Simulation = s_->getParameters();
    opt.IonBeam = s_->getSource().getParameters();
    opt.Target = s_->getTarget().getDescription();
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

int mcdriver::exec(progress_callback cb, uint msInterval, void *callback_user_data)
{
    using namespace std::chrono_literals;

    int nthreads = par_.threads;
    if (nthreads < 1) nthreads = 1;

    // TIMING
    start_time_ = std::time(nullptr);
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
        thread_ion_count_.assign(nthreads,0);
        ion_count_ = 0;
        while(ion_count_ < par_.max_no_ions) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msInterval));
            int i=0;
            for(; i<nthreads; i++) {
                thread_ion_count_[i] = sims[i]->ions_done();
                ion_count_ += thread_ion_count_[i];
            }
            cb(*this,callback_user_data);
        }
    }

    // wait for threads to finish...
    for(int i=0; i<nthreads; i++) threads[i]->join();

    // consolidate results
    for(int i=1; i<nthreads; i++)
        sims[0]->merge(*(sims[i]));

    // CALC TIME/ion CLOCK_PROCESS_CPUTIME_ID
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_end); // POSIX
    double t_secs = 1. * (t_end.tv_sec - t_start.tv_sec) / nthreads;
    t_secs += 1.e-9 * (t_end.tv_nsec - t_start.tv_nsec) / nthreads;
    ips_ = s_->getTally().Nions()/t_secs;
    end_time_ = std::time(nullptr);

    // delete threads
    for(int i=0; i<nthreads; i++) {
        delete threads[i];
    }
    // delete simulation clones
    for(int i=1; i<nthreads; i++) {
        delete sims[i];
    }

    return 0;
}

int mcdriver::options::validate()
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
    if (Simulation.scattering_calculation==mccore::InvalidScatteringOption)
        throw std::invalid_argument("Invalid Simulation.scattering_calculation");

    if (Simulation.flight_path_type==mccore::InvalidPath)
        throw std::invalid_argument("Invalid Simulation.flight_path_type");
    if (Simulation.flight_path_type==mccore::Constant &&
        Simulation.flight_path_const<=0.f)
        throw std::invalid_argument("Simulation.flight_path_type is \"Constant\" but Simulation.flight_path_const is negative.");

    if (Simulation.nrt_calculation==mccore::NRT_InvalidOption)
        throw std::invalid_argument("Invalid Simulation.nrt_calculation");

    if (Simulation.screening_type==Screening::None)
        throw std::invalid_argument("Invalid Simulation.screening_type");

    if (Simulation.straggling_model==StragglingModel::Invalid)
        throw std::invalid_argument("Invalid Simulation.straggling_model");

    // Ion source
    if (IonBeam.ion_distribution==ion_beam::InvalidIonDistribution)
        throw std::invalid_argument("Invalid IonBeam.ion_distribution");


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

    // Check Material descriptors
    for(auto p : Target.materials) {
        auto md = p.second;
        auto mname = p.first;
        if (md.density <= 0.f) {
            std::stringstream msg;
            msg << "Zero or negative density in material";
            msg << mname;
            throw std::invalid_argument(msg.str());
        }
        CHECK_EMPTY_VEC(md,Z,mname)
        CHECK_EMPTY_VEC(md,M,mname)
        CHECK_EMPTY_VEC(md,X,mname)
        CHECK_EMPTY_VEC(md,Ed,mname)
        CHECK_EMPTY_VEC(md,El,mname)
        CHECK_EMPTY_VEC(md,Es,mname)
        CHECK_EMPTY_VEC(md,Er,mname)

        int natoms = md.Z.size();
        if (std::any_of(md.Z.begin(),
                        md.Z.end(),
                        [](int z){ return (z < 1) || (z > elements::max_atomic_num);}))
        {
            throw std::invalid_argument("Invalid Z number in material" + mname);
        }

        CHECK_VEC_SIZE(md,M,natoms,mname)
        CHECK_VEC_SIZE(md,X,natoms,mname)
        CHECK_VEC_SIZE(md,Ed,natoms,mname)
        CHECK_VEC_SIZE(md,El,natoms,mname)
        CHECK_VEC_SIZE(md,Es,natoms,mname)
        CHECK_VEC_SIZE(md,Er,natoms,mname)
        CHECK_VEC_ZEROorNEG(md,M,mname)
        CHECK_VEC_ZEROorNEG(md,X,mname)
        CHECK_VEC_ZEROorNEG(md,Ed,mname)
        CHECK_VEC_ZEROorNEG(md,El,mname)
        CHECK_VEC_ZEROorNEG(md,Es,mname)
        CHECK_VEC_ZEROorNEG(md,Er,mname)
    }

    // Check Region descriptors
    for(auto p : Target.regions) {
        auto rd = p.second;
        auto rname = p.first;

        // check valid material
        if (Target.materials.find(rd.material_id)==Target.materials.end()) {
            std::stringstream msg;
            msg << "Region " << rname << " has invalid material_id: ";
            msg << rd.material_id;
            throw std::invalid_argument(msg.str());
        }

        // check max > min
        for (int i=0; i<3; i++) {
            if (rd.min[i] >= rd.max[i]) {
                static char axis[] = { 'x', 'y', 'z' };
                std::stringstream msg;
                msg << "Region " << rname << " has invalid axis limits: ";
                msg << axis[i] << "_min(" << rd.min[i] << ") >= ";
                msg << axis[i] << "_max(" << rd.max[i] << ")";
                throw std::invalid_argument(msg.str());
            }
        }

        // check that the region is within the simulation volume
        box3D rbox, // region box
            sbox; // simulation box
        rbox.min() = rd.min;
        rbox.max() = rd.max;
        sbox.min() = vector3(0,0,0);
        sbox.max() = vector3(Target.cell_count(0)*Target.cell_size(0),
                             Target.cell_count(1)*Target.cell_size(1),
                             Target.cell_count(2)*Target.cell_size(2));
        rbox = sbox.intersection(rbox);
        if (rbox.isEmpty()) {
            std::stringstream msg;
            msg << "Region " << rname << " does not intersect ";
            msg << "the simulation volume.";
            throw std::invalid_argument(msg.str());
        }
    }

    return 0;
}

mccore* mcdriver::options::createSimulation() const
{
    mccore* S = new mccore(Simulation);
    //S->setOutputOptions(Output);
    S->getSource().setParameters(IonBeam);

    target& T = S->getTarget();

    std::unordered_map<std::string, material*> materials_map;

    for(auto p : Target.materials) {
        material* m = T.addMaterial(p.first.c_str());
        materials_map[p.first] = m;
        const material::material_desc_t& md = p.second;
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

    for(auto p : Target.regions) {
        const target::region& rd = p.second;
        box3D box;
        box.min() = rd.min;
        box.max() = rd.max;
        T.fill(box,materials_map[rd.material_id]);
    }

    return S;
}
