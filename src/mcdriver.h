#ifndef MCDRIVER_H
#define MCDRIVER_H

#include "mccore.h"

#include <ctime>

/**
 * \defgroup Driver Driver
 *
 * \brief The core of the Monte-Carlo BCA ion transport simulation
 *
 * @{
 *
 * @ingroup MC
 *
 * @}
 *
 *
 */

/**
 * @brief The mcdriver class facilitates the setup and running of a simulation.
 *
 * @ingroup Driver
 */
class mcdriver
{
public:
    /**
     * @brief Driver parameters/options
     */
    struct parameters {
        /// Ions to run
        unsigned int max_no_ions{100};
        /// Number of threads to use
        int threads{1};
        /// Seed for the random number generator
        std::vector<unsigned int> seeds;
    };


    struct output_options {
        /// Simulation title
        std::string title{"Ion Simulation"};
        std::string OutputFileBaseName{"out"};
        int storage_interval{1000};
        int store_transmitted_ions{0};
        int store_range_3d{0};
        int store_ion_paths{0};
        int store_path_limit{100};
        int store_recoil_cascades{0};
        int store_path_limit_recoils{4};
        int store_pka{0};
        int store_dedx{1};
    };

    typedef void (*progress_callback)(const mcdriver& v, void* p);

protected:

    // timing
    double ips_; // ions/s
    std::time_t start_time_, end_time_;

    parameters par_;
    output_options out_opts_;

    mccore* s_;

    std::string outFileName(const char* type, int thread_id);

    std::vector<uint> thread_ion_count_;
    uint ion_count_;



public:
    mcdriver();
    ~mcdriver();

    void getOptions(options& opt) const;
    void setOptions(const options& o);
    std::string outFileName() const;

    const output_options& outputOptions(const output_options& opts) const
    { return out_opts_; }
    void setOutputOptions(const output_options& opts)
    { out_opts_ = opts; }

    const mccore* getSim() { return s_; }

    double ips() const { return ips_; }
    int nThreads() const { return par_.threads; }

    int save();

    const std::vector<uint>& thread_ion_count() const { return thread_ion_count_; }
    uint ion_count() const { return ion_count_; }

    int exec(progress_callback cb = nullptr, uint msInterval = 1000, void* callback_user_data = 0);

};


struct options
{
    mcdriver::parameters Driver;
    mcdriver::output_options Output;
    mccore::parameters Simulation;
    ion_beam::parameters IonBeam;
    target::target_desc_t Target;

    int parseJSON(std::istream& js);
    void printJSON(std::ostream& os) const;

    int validate();
    mccore* createSimulation() const;

};

#endif // MCDRIVER_H
