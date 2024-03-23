#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "random_vars.h"
#include "ion_beam.h"
#include "ion.h"
#include "xs.h"
#include "tally.h"
#include "event_stream.h"
#include "dedx.h"

#include <atomic>

class out_file;
class options;

/**
 * \defgroup MC Simulation components
 *
 * \brief The main parts of a simulation.
 *
 * @{
 *
 * @}
 */

/**
 * @brief The simulation class is used for preparing and running an ion transport simulation.
 *
 * @ingroup MC
 */
class simulation
{
public:

    /**
     * @brief Type of simulation
     */
    enum simulation_type_t {
        FullCascade = 0,    /**< Full Damage Cascade, follow recoils */
        Invalid1 = 1,       // Compatibility. Iradina does not have simulation_type=1 & 2
        Invalid2 = 2,
        IonsOnly = 3,       /**< Ions only. Recoils are not followed. Damage estimate by NRT */
        Cascades = 4,       /**< Cascades generated and followed */
        InvalidSimulationType = -1
    };

    /**
     * @brief Detail of NRT implementation in multielement materials
     */
    enum nrt_calculation_t {
        NRT_element = 0, /**< NRT calculated for using Ed of struck atom (similar to SRIM) */
        NRT_average = 1,  /**< NRT calculated with average Ed (J.-P. Crocombette 2019) */
        NRT_InvalidOption = -1
    };

    /**
     * @brief Determines how ion scattering is simulated
     */
    enum scattering_calculation_t {
        Corteo4bit = 0,     /**< Using 4bit corteo tabulated ZBL cross-section */
        Corteo6bit = 1,     /**< Using 6bit corteo tabulated ZBL cross-section */
        ZBL_MAGICK = 2,     /**< Using the MAGIC formula for the ZBL cross-section (as in SRIM) */
        InvalidScatteringOption = -1
    };

    /**
     * @brief Flight path selection algorithm
     *
     * Determines the algorithm to select the
     * free flight path \f$\ell\f$ between collisions.
     */
    typedef enum {
        AtomicSpacing = 1,  /**< Constant, equal to interatomic distance */
        Constant = 2,       /**< Constant, equal to user supplied value */
        SRIMlike = 3,        /**< Algorithm similar to SRIM */
        MendenhallWeller = 4, /**< Algorithm from Mendenhall-Weller NIMB2005*/
        MyFFP = 5,            /**< Other algorithm */
        InvalidPath = -1
    } flight_path_type_t;

    /**
     * @brief The model used to calculate electronic straggling
     */
    enum straggling_model_t {
        NoStraggling = 0,       /**< No straggling calculated */
        BohrStraggling,         /**< Bohr straggling model */
        ChuStraggling,          /**< Chu straggling model */
        YangStraggling,         /**< Yang straggling model */
        InvalidStraggling = -1
    };

    /**
     * @brief Simulation parameters/options
     */
    struct parameters {
        /// Simulation title
        std::string title;
        /// Ions to run
        unsigned int max_no_ions{100};
        /// Type of the simulation
        simulation_type_t simulation_type{FullCascade};
        /// Way to calculate NRT vacancies in multielement materials
        nrt_calculation_t nrt_calculation{NRT_element};
        /// Way to calculate nuclear scattering
        scattering_calculation_t scattering_calculation{Corteo4bit};
        /// Free flight path selection algorithm
        flight_path_type_t flight_path_type{AtomicSpacing};
        /// Model to use for calculating electronic straggling
        straggling_model_t straggling_model{YangStraggling};
        /// The constant flight path (for relevant algorithm) [nm]
        float flight_path_const{0.1};
        /// Minimum energy cutoff for ion transport [eV]
        float min_energy{1.f};
        /// Number of threads to use
        int threads{1};
        /// Seed for the random number generator
        std::vector<unsigned int> seeds;
    };

    struct output_options {
        std::string OutputFileBaseName{"iradina++"};
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

    typedef void (*progress_callback)(const std::vector<uint>& v);

protected:

    ion_queue q_;

    // simulation paramenters
    parameters par_;
    output_options out_opts_;

    // components
    ion_beam* source_;
    target* target_;
    tally tally_;
    random_vars rng;
    event_stream pka_stream_, exit_stream_;
    pka_event pka;
    exit_event exit_ev;

    // ref counter
    std::shared_ptr<int> ref_count_;
    int thread_id_;
    uint count_offset_;

    // helper variable for flight path calc
    std::vector<float> sqrtfp_const;

    // Scattering cross-sections
    Array2D< abstract_xs_lab* > scattering_matrix_;
    // Stopping & Straggling Tables
    Array3Df dedx_; // stopping data (atoms x materials x energy)
    Array2Df dedx1; // proton stopping (materials x energy)
    Array3Df de_strag_; // straggling data (atoms x materials x energy)
    Array3Df max_fp_; // max fp for 1% dEdx (atoms x materials x energy)
    //
    Array3Df max_impact_par_;

    // timing
    double ips_; // ions/s
    std::atomic_uint nion_thread_; // thread safe simulated ion counter

    friend class out_file;

public:
    simulation();
    explicit simulation(const parameters& p);
    simulation(const simulation& S);
    ~simulation();

    void setIonBeam(const ion_beam::parameters& p) {
        source_->setParameters(p);
    }

    const output_options& outputOptions(const output_options& opts) const
    { return out_opts_; }
    void setOutputOptions(const output_options& opts)
    { out_opts_ = opts; }
    void setOutputFileBaseName(const char* n)
    { out_opts_.OutputFileBaseName = n; }
    void store_pka(bool on = true)
    { out_opts_.store_pka = on; }

    void setMaxIons(unsigned int n) { par_.max_no_ions = n; }

    void getOptions(options& opt) const;

    const parameters& getParameters() const { return par_; }
    const std::string& title() const { return par_.title; }
    simulation_type_t simulationType() const { return par_.simulation_type; }
    flight_path_type_t flightPathType() const { return par_.flight_path_type; }
    straggling_model_t stragglingModel() const { return par_.straggling_model; }
    int nThreads() const { return par_.threads; }

    void setTitle(const char* t) { par_.title = t; }
    void setSimulationType(simulation_type_t v) { par_.simulation_type = v; }
    void setFlightPathType(flight_path_type_t v) { par_.flight_path_type = v; }
    void setStragglingModel(straggling_model_t v) { par_.straggling_model = v; }
    void setNThreads(int n) { par_.threads = n; }
    void setSeeds(const std::vector<uint>& s) { par_.seeds = s; }

    double ips() const { return ips_; }

    target& getTarget() { return *target_; }
    const target& getTarget() const { return *target_; }
    ion_beam& getSource() { return *source_; }
    const ion_beam& getSource() const { return *source_; }
    const tally& getTally() const { return tally_; }

    int saveTallys();
    int init();
    int exec(progress_callback cb = 0, uint msInterval = 1000);
    void seed(unsigned int s) { rng.seed(s); }

protected:

    void addTally(const tally& t) { tally_ += t; }

    ion *new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio, tally &t, pka_event* pka);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;

    std::string outFileName(const char* type);

    int transport(ion* i, tally& t, pka_event* ev = nullptr);
    int flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp);
    int run();
    uint ions_done() { return nion_thread_.exchange(0); }
    float doDedx(const ion* i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl);

};

// helper 1D interpolation function
//   CI : a corteo index
//   array : a sequence container that can be accessed by array[i]
// On input
//   x : the point where we require an interpolation
//   i : should be CI(x)
//   data : the interpolation table
template<class CI, class array>
float interp1d(float x, CI i, array data)
{
    if (i==i.rbegin()) return data[i]; // x out of range
    float y1 = data[i], x1 = *i++;
    float y2 = data[i], x2 = *i;
    return y1 + (y2-y1)*(x-x1)/(x2-x1);
}

// helper 1D interpolation function for dedx
//   CI : a corteo index
//   array : a sequence container that can be accessed by array[i]
// On input
//   x : the point where we require an interpolation
//   i : should be CI(x)
//   data : the interpolation table
inline float interp_dedx(float E, const float* data)
{
    if (E <= dedx_index::minVal) return data[0];
    if (E >= dedx_index::maxVal) return data[dedx_index::size - 1];
    dedx_index i(E);
    float y1 = data[i], x1 = *i++;
    float y2 = data[i], x2 = *i;
    return y1 + (y2-y1)*(E-x1)/(x2-x1);
}


#endif // SIMULATION_H
