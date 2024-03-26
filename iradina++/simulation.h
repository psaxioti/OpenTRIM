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
 * \defgroup Core Core of the Monte-Carlo ion transport simulation
 *
 * \brief The core of the Monte-Carlo ion transport simulation
 *
 * @{
 *
 * @}
 */

/**
 * @brief The mccore class defines the core Monte-Carlo ion transport simulation algorithms.
 *
 * @ingroup MC
 */
class mccore
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
    };



protected:

    ion_queue q_;

    // simulation paramenters
    parameters par_;

    // components
    ion_beam* source_;
    target* target_;
    tally tally_, dtally_, tion_;
    random_vars rng;
    event_stream pka_stream_, exit_stream_;
    pka_event pka;
    exit_event exit_ev;

    // ref counter
    std::shared_ptr<int> ref_count_;

    // ion counters
    uint max_no_ions_;
    uint count_offset_;
    std::atomic_uint nion_thread_; // thread safe simulated ion counter

    // helper variable for flight path calc
    std::vector<float> sqrtfp_const;

    // Scattering cross-sections
    ArrayND< abstract_xs_lab* > scattering_matrix_;
    // Stopping & Straggling Tables
    ArrayNDf dedx_; // stopping data (atoms x materials x energy)
    ArrayNDf dedx1; // proton stopping (materials x energy)
    ArrayNDf de_strag_; // straggling data (atoms x materials x energy)
    ArrayNDf max_fp_; // max fp for 1% dEdx (atoms x materials x energy)
    //
    ArrayNDf max_impact_par_;

public:
    mccore();
    explicit mccore(const parameters& p);
    mccore(const mccore& S);
    ~mccore();

    void setMaxIons(unsigned int n) { max_no_ions_ = n; }
    void setCountOffset(uint n) { count_offset_ = n; }
    uint ions_done() { return nion_thread_.exchange(0); }

    const parameters& getParameters() const { return par_; }

    target& getTarget() { return *target_; }
    const target& getTarget() const { return *target_; }
    ion_beam& getSource() { return *source_; }
    const ion_beam& getSource() const { return *source_; }
    const tally& getTally() const { return tally_; }

    int open_pka_stream(const char* fname) {
        return pka_stream_.open(fname, pka.size());
    }
    const event_stream& pka_stream() { return pka_stream_; }
    int open_exit_stream(const char* fname) {
        return exit_stream_.open(fname, exit_ev.size());
    }
    const event_stream& exit_stream() { return exit_stream_; }
    void remove_stream_files();

    ArrayNDf dedx() const { return dedx_; }
    ArrayNDf de_strag() const { return de_strag_; }
    ArrayNDf max_impact_par() const { return max_impact_par_; }
    ArrayNDf max_fp() const { return max_fp_; }

    int init();

    void seed(unsigned int s) { rng.seed(s); }

    void merge(const mccore& other);

    int run();

protected:



    ion *new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio, tally &t, pka_event* pka);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;



    int transport(ion* i, tally& t, pka_event* ev = nullptr);
    int flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp);


    float doDedx(const ion* i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl);

    static float interp_dedx(float E, const float* data);

};


// helper 1D interpolation function for dedx
//   CI : a corteo index
//   array : a sequence container that can be accessed by array[i]
// On input
//   x : the point where we require an interpolation
//   i : should be CI(x)
//   data : the interpolation table
inline float mccore::interp_dedx(float E, const float* data)
{
    if (E <= dedx_index::minVal) return data[0];
    if (E >= dedx_index::maxVal) return data[dedx_index::size - 1];
    dedx_index i(E);
    float y1 = data[i], x1 = *i++;
    float y2 = data[i], x2 = *i;
    return y1 + (y2-y1)*(E-x1)/(x2-x1);
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
inline int mccore::getDEtables(const atom* z1, const material* m,
                            const float* &dedx, const float* &de_stragg) const
{
    int ia = z1->id();
    int im = m->id();
    dedx = &dedx_(ia,im,0);
    de_stragg = &de_strag_(ia,im,0);
    return 0;
}




#endif // SIMULATION_H
