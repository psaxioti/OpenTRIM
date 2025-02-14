#include "mcdriver.h"

#include <cxxopts.hpp>

#include <iostream>
#include <fstream>
#include <iomanip>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// helper class for getting real-time info
// for the running simulation
class running_sim_info {

    // timing
    typedef std::chrono::high_resolution_clock hr_clock_t;
    typedef hr_clock_t::time_point time_point;
    // sim data
    time_point tstart_; // time when sim started
    size_t nstart_; // # of ions when started
    size_t ncurr_; // # of ions current
    size_t ntarget_; // max # of ions
    int progress_; // 0..100
    double elapsed_;
    double total_elapsed_;
    double etc_; // estimated time to completion (s)
    double ips_; // ions per second

public:
    // init : called before simulation starts
    void init(const mcdriver& D);
    // update : called during simulation run
    void update(const mcdriver& D);
    // print cli progress bar
    void print();
    // getters
    int progress() const { return progress_; }
    double elapsed() const { return elapsed_; }
    double ips() const { return ips_; }
    double etc() const { return etc_; }
    size_t nions() const { return ncurr_; }
};

running_sim_info info;

void progress_callback(const mcdriver& d, void* )
{
    info.update(d);
    info.print();
}

int main(int argc, char* argv[])
{
    std::string program_name(PROJECT_NAME);
    std::transform(program_name.begin(), program_name.end(), program_name.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    cxxopts::Options cli_options(program_name, PROJECT_DESCRIPTION);

    cli_options.add_options()
        ("n","Number of histories to run (overrides config input)", cxxopts::value<int>())
        ("j","Number of threads (overrides config input)", cxxopts::value<int>())
        ("s,seed","random generator seed (overrides config input)",cxxopts::value<int>())
        ("i,input","input HDF5 file name",cxxopts::value<std::string>())
        ("f","JSON config file",cxxopts::value<std::string>())
        ("o,output","output HDF5 file name (overrides config input)",cxxopts::value<std::string>())
        ("t,template","print a template JSON config to stdout")
        ("v,version","Display version information")
        ("h,help","Display short help message");

    int n(-1), j(-1), s(-1);
    std::string input_config_file, input_file, output_file;

    try {
        auto result = cli_options.parse(argc,argv);

        if (result.count("help")) {
            cout << cli_options.help() << endl;
            return 0;
        }
        if (result.count("version")) {
            cout << PROJECT_NAME << " version " << PROJECT_VERSION << endl;
            cout << "Build time: " << BUILD_TIME << endl;
            cout << "Compiler: " << COMPILER_ID << " v" << COMPILER_VERSION << " on " SYSTEM_ID << endl;
            return 0;
        }
        if (result.count("template")) {
            mcdriver::options opt;
            opt.printJSON(cout);
            return 0;
        }
        if (result.count("n")) n = result["n"].as<int>();
        if (result.count("j")) j = result["j"].as<int>();
        if (result.count("s")) s = result["s"].as<int>();
        if (result.count("f")) input_config_file = result["f"].as<std::string>();
        if (result.count("i")) input_file = result["i"].as<std::string>();
        if (result.count("o")) output_file = result["o"].as<std::string>();
    }
    catch(const cxxopts::exceptions::exception& e)
    {
        cerr << "Error parsing command line: " << e.what() << endl;
        return -1;
    }

    if (input_config_file.empty() && input_file.empty()) {
        cerr << "Error: No input file." << endl;
        return -1;
    }

    if (!input_config_file.empty() && !input_file.empty()) {
        cerr << "Warning: JSON config ignored (HDF5 input will be used)." << endl;
        input_config_file.clear();
    }

    mcdriver D;

    if (!input_config_file.empty()) {

        cout << "Parsing JSON config from " << input_config_file << endl;

        std::ifstream is(input_config_file);
        mcdriver::options opt;
        if (opt.parseJSON(is,true,&cerr)!=0) return -1;

        // cli overrides
        if (n>0) opt.Driver.max_no_ions = n;
        if (j>0) opt.Driver.threads = j;
        if (s>0) opt.Driver.seed = s;
        if (!output_file.empty()) opt.Output.outfilename = output_file;

        D.init(opt);

    } else {

        cout << "Loading simulation from " << input_file << endl;

        if (D.load(input_file, &cerr)!=0) return -1;

        // cli overrides
        mcdriver::parameters par = D.driverOptions();
        if (n>0) par.max_no_ions = n;
        if (j>0) par.threads = j;
        D.setDriverOptions(par);

        mcdriver::output_options opts = D.outputOptions();
        if (!output_file.empty()) {
            opts.outfilename = output_file;
            D.setOutputOptions(opts);
        }
    }

    cout << "Starting simulation '" << D.outputOptions().title << "'..." << endl << endl;

    info.init(D);
    info.print();

    D.exec(progress_callback,500);
    
    const mcdriver::run_data& rd = D.run_history().back();
    cout << endl << endl
         << "Completed " << rd.total_ion_count << " ion histories." << endl;
    cout << "Cpu time (s) = " << rd.cpu_time << endl;
    cout << "Ions/cpu-s = " << rd.ips << endl;
    cout << "Real time (s) = " << info.elapsed() << endl;
    cout << "Ions/real-s = " << info.ips() << endl;
    cout << "Storing results in " << D.outFileName() << " ...";
    cout.flush();
    D.save(D.outFileName(),&cerr);
    cout << " OK." << endl;

    return 0;
}

void running_sim_info::init(const mcdriver& d)
{
    tstart_ = hr_clock_t::now();
    nstart_ = d.getSim()->ion_count();
    ncurr_ = nstart_;
    ntarget_ = d.driverOptions().max_no_ions;
    progress_ = int(100.0*ncurr_/ntarget_);
    elapsed_ = 0.;
    etc_ = 0;
    ips_ = 0.;
}

void running_sim_info::update(const mcdriver& d)
{
    time_point t = hr_clock_t::now();
    ncurr_ = d.getSim()->ion_count();
    // floating-point duration: no duration_cast needed
    const std::chrono::duration<double> fp_sec = t - tstart_;
    elapsed_ = fp_sec.count();
    ips_ = (ncurr_ - nstart_)/elapsed_;
    etc_ = ips_ > 0 ? (ntarget_ - ncurr_)/ips_ : 0;
    progress_ = int(100.0*ncurr_/ntarget_);
}

void running_sim_info::print()
{
    cout << '\r';
    cout << '[';

    int n = progress_ >> 1; // bar steps at 2%/step

    // add n #
    for (int j = 0; j < n; ++j) cout << '#';

    // add spaces
    for (int j = 0; j < 50-n; ++j) cout << ' ';

    // add closing bracket & trailing percentage characters
    cout << ']';

    char buff[8];
    sprintf(buff,"%3d%%",progress_);
    cout << buff;

    cout.flush();
}





