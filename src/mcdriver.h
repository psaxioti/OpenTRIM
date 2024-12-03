#ifndef MCDRIVER_H
#define MCDRIVER_H

#include "mccore.h"

#include <ctime>
#include <thread>

/**
 * \defgroup Driver Driver
 *
 * \brief Classes for setting up and running a simulation
 *
 * The \ref mcdriver and its sub-classes can be used to perform
 * the following tasks:
 * 
 * - Parse the configuration options from JSON 
 * - Validate the configuration
 * - Create the \ref mccore object, generate the geometry and load all options
 * - Run the simulation
 * - Save the results
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
 * A typical usage scenario would be:
 *
 * @code{.cpp}
 * #include "mcdriver.h"
 * #include <iostream>
 * 
 * mcdriver::options opt;
 * opt.parseJSON(std::cin);
 * mcdriver d;
 * d.setOptions(opt);
 * d.exec();
 * d.save();
 * @endcode
 *
 * @ingroup Driver
 */
class mcdriver
{
public:
    /// mcdriver parameters for running the simulation
    struct parameters {
        /// Ions to run
        unsigned int max_no_ions{100};
        /// Number of threads to use
        int threads{1};
        /// Seed for the random number generator
        unsigned int seed{123456789};
    };

    /// mcdriver output options
    struct output_options {
        /// Simulation title
        std::string title{"Ion Simulation"};
        /// Base name for the output file
        std::string OutputFileBaseName{"out"};
        /// Interval in sec to store the output @todo
        int storage_interval{1000};
        /// Store ion exit events
        int store_transmitted_ions{0};
        /// @todo
        int store_range_3d{0};
        /// @todo
        int store_ion_paths{0};
        /// @todo
        int store_path_limit{100};
        /// @todo
        int store_recoil_cascades{0};
        /// @todo
        int store_path_limit_recoils{4};
        /// Store the pka events
        int store_pka{0};
        /// Store electron energy loss data
        int store_dedx{1};
    };

    /// Typedef for a function to be called during simulation execution
    typedef void (*progress_callback)(const mcdriver& v, void* p);

    /**
     * @brief mcdriver::options is a helper class for parsing and validating all simulation options
     */
    struct options
    {
        mcdriver::parameters Driver;
        mcdriver::output_options Output;
        mccore::parameters Simulation;
        mccore::transport_options Transport;
        ion_beam::parameters IonBeam;
        target::target_desc_t Target;

        /**
         * @brief Parse simulation options from JSON formatted input
         *
         * For a full list of available options and details on JSON
         * formatting see \ref json_config.
         *
         * The function first parses the whole JSON string. On formatting errors
         * the function stops and prints an error message to stderr.
         *
         * After that validate() is called to check the given options. Errors
         * are again reported to stderr.
         *
         * @param js a JSON formatted input stream
         * @param doValidation if true the function calls validate()
         * @return 0 if succesfull, negative value otherwise
         */
        int parseJSON(std::istream& js, bool doValidation = true, std::ostream *os = nullptr);

        /// Pretty print JSON formattet options to a stream
        void printJSON(std::ostream& os) const;

        /// Return options as a JSON string
        std::string toJSON() const;

        /**
         * @brief Validate the simulation options
         *
         * A number of checks are performed
         * including
         * - correct parameter range
         * - allowed parameter combinations
         * - target definition (geometry, materials, regions)
         *
         * On error, a std::invalid_argument exception is thrown.
         * exception::what() returns a relevant error message.
         *
         * @return
         */
        int validate();

        /**
         * @brief Create a simulation object from the given options
         * @return a pointer to a mccore object
         */
        mccore* createSimulation() const;

    };

protected:

    // timing
    double ips_; // ions/s
    std::time_t start_time_, end_time_;
    struct timespec t_start, t_end;

    // driver parameters
    parameters par_;
    output_options out_opts_;

    // the simulation object
    mccore* s_;

    // The following 2 are populated when the
    // simulation runs
    // 1. threads
    std::vector< std::thread > thread_pool_;
    // 2. simulation execution clones
    std::vector< mccore* > sim_clones_;

public:

    mcdriver();
    ~mcdriver();


    void getOptions(options& opt) const;
    void setOptions(const options& o);
    std::string outFileName() const;

    const output_options& outputOptions() const
    { return out_opts_; }
    void setOutputOptions(const output_options& opts)
    { out_opts_ = opts; }
    const parameters& driverOptions() const
    { return par_; }
    void setDriverOptions(const parameters& opts)
    { par_ = opts; }

    const mccore* getSim() const { return s_; }

    /// Signal a running simulation to abort
    void abort();
    /// Wait for a running simulation to finish
    void wait();
    /// Reset the simulation object
    void reset();

    /// total ions/s
    double ips() const { return ips_; }
    /// number of threads
    int nThreads() const { return par_.threads; }

    /**
     * @brief Saves all data and results in a HDF5 file
     *
     * For details on the structure of the output file see \ref out_file.
     *
     * @param h5filename the output file name
     * @param os optional stream pointer to write any error messages
     * @return 0 if succesful
     */
    int save(const std::string& h5filename, std::ostream* os = nullptr);

    /**
     * @brief Load a simulation from a HDF5 file
     * @param h5filename the name of the file
     * @param os optional stream pointer to write any error messages
     * @return 0 if succesfull
     */
    int load(const std::string& h5filename, std::ostream* os = nullptr);

    int exec(progress_callback cb = nullptr, size_t msInterval = 1000, void* callback_user_data = 0);
};




#endif // MCDRIVER_H
