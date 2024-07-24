#ifndef _MCCORE_H_
#define _MCCORE_H_

#include "target.h"
#include "random_vars.h"
#include "ion_beam.h"
#include "ion.h"
#include "xs.h"
#include "tally.h"
#include "event_stream.h"
#include "dedx.h"

#include <atomic>

/**
 * \defgroup MC libiradinapp shared library
 *
 * \brief A library containing the main ion transport simulation code.
 * 
 * The building blocks for the Monte-Carlo ion transport simulation are: \ref Geometry, 
 * \ref Ions "Ion transport", \ref RNG and \ref target "Target definition".
 * 
 * All these are combined into the simulation \ref Core.
 * 
 * The \ref Driver part facilitates the process of loading the configuration,
 * setting up the simulation, executing on multiple threads and adding up and storing the results.
 *
 * 
 */

/**
 * \defgroup Core Core
 *
 * \brief The core of the Monte-Carlo BCA ion transport simulation
 *
 * @{
 *
 *
 * @ingroup MC
 *
 * @}
 *
 *
 */

/**
 * @brief The mccore class defines the core Monte-Carlo BCA ion transport simulation algorithms.
 *
 * @ingroup Core
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
     * @brief Determines how ion scattering is simulated
     */
    enum scattering_calculation_t {
        Corteo4bitTable = 0,     /**< Using 4bit corteo-tabulated scattering integrals */
        ZBL_MAGICK = 1,     /**< Using the MAGIC formula for the scattering integrals (as in SRIM) */
        InvalidScatteringOption = -1
    };

    /**
     * @brief Determines electronic energy loss calculation
     */
    enum eloss_calculation_t {
        EnergyLossOff = 0, /**< Electronic energy loss disabled */
        EnergyLoss = 1,  /**< Only electronic energy loss is calculated */
        EnergyLossAndStraggling = 2, /**< Both energy loss & straggling are calculated */
        InvalidEnergyLoss = -1
    };

    /**
     * @brief Flight path selection algorithm
     *
     * Determines the algorithm to select the
     * free flight path \f$\ell\f$ between collisions.
     */
    enum flight_path_type_t {
        AtomicSpacing = 1,  /**< Constant, equal to interatomic distance */
        Constant = 2,       /**< Constant, equal to user supplied value */
        MendenhallWeller = 3, /**< Algorithm from Mendenhall-Weller NIMB2005*/
        IPP = 4,            /**< IPP algorithm */
        InvalidPath = -1
    };

    /**
     * @brief Detail of NRT implementation in multielement materials
     */
    enum nrt_calculation_t {
        NRT_element = 0, /**< NRT calculated using Ed of struck atom (similar to SRIM) */
        NRT_average = 1,  /**< NRT calculated with average Ed (J.-P. Crocombette 2019) */
        NRT_InvalidOption = -1
    };

    /**
     * @brief Simulation parameters/options
     */
    struct parameters {        
        /// Type of the simulation
        simulation_type_t simulation_type{FullCascade};
        /// screeninig
        Screening screening_type{Screening::ZBL};
        /// Way to calculate nuclear scattering
        scattering_calculation_t scattering_calculation{Corteo4bitTable};
        /// Free flight path selection algorithm
        flight_path_type_t flight_path_type{AtomicSpacing};
        /// Electronic energy loss calculation option
        eloss_calculation_t eloss_calculation{EnergyLoss};
        /// Model to use for calculating electronic straggling (if calculated)
        StragglingModel straggling_model{StragglingModel::Bohr};
        /// Way to calculate NRT vacancies in multielement materials
        nrt_calculation_t nrt_calculation{NRT_element};
        /// The constant flight path (for algorithm flight_path_type_t=Constant) [nm]
        float flight_path_const{0.1};
        /// Minimum energy cutoff for ion transport [eV]
        float min_energy{1.f};
        /// Minimum recoil energy
        float min_recoil_energy{1.f};
        /// Scattering at shorter distances than interatomic distance
        bool allow_sub_ml_scattering{false};
        /// Max mfp
        float max_mfp{std::numeric_limits<float>::max()};
        /// Max dE/E per mfp, dE = dEdx*mfp
        float max_rel_eloss{0.05f};
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

    // helper variables for flight path calc
    std::vector<float> sqrtfp_const;
    std::vector<float> ip0;

    // Scattering cross-sections
    ArrayND< abstract_xs_lab* > scattering_matrix_;

    // Electronic Stopping & Straggling Tables
    ArrayND< dedx_interp* > dedx_; // stopping data (atoms x materials)
    ArrayND< straggling_interp* > de_strag_; // straggling data (atoms x materials)

    // flight path selection par
    ArrayNDf mfp_, ipmax_, fp_max_, Tcutoff_;

public:
    mccore();
    explicit mccore(const parameters& p);
    mccore(const mccore& S);
    ~mccore();

    void setMaxIons(unsigned int n) { max_no_ions_ = n; }
    void setCountOffset(uint n) { count_offset_ = n; }
    uint ions_done() { return nion_thread_.exchange(0); }
    size_t ion_que_size() const { return q_.size(); }

    /// Returns the core simulation parameters
    const parameters& getParameters() const { return par_; }

    /// Return a reference to the \ref target object
    target& getTarget() { return *target_; }
    /// Return a const reference to the target object
    const target& getTarget() const { return *target_; }
    /// Return a reference to the source (ion_beam) object
    ion_beam& getSource() { return *source_; }
    /// Return a const reference to the source (ion_beam) object
    const ion_beam& getSource() const { return *source_; }
    /// Return a const reference to the tally object
    const tally& getTally() const { return tally_; }
    /// Return a const reference to the tally object
    const tally& getTallyVar() const { return dtally_; }
    /// Open file stream to store \ref pka_event data
    int open_pka_stream(const char* fname) {
        return pka_stream_.open(fname, pka);
    }
    /// Return reference to the pka stream
    const event_stream& pka_stream() { return pka_stream_; }
    /// Open file stream to store \ref exit_event data
    int open_exit_stream(const char* fname) {
        return exit_stream_.open(fname, exit_ev);
    }
    /// Return reference to the exit stream
    const event_stream& exit_stream() { return exit_stream_; }
    /// Removes stream files from the filesystem
    void remove_stream_files();

    // dedx interpolators
    ArrayND<dedx_interp*> dedx() const { return dedx_; }
    ArrayND<straggling_interp*> de_strag() const { return de_strag_; }

    // Tables of flight path selection parameters 
    ArrayNDf mfp() const { return mfp_; }
    ArrayNDf ipmax() const { return ipmax_; }
    ArrayNDf fpmax() const { return fp_max_; }
    ArrayNDf Tcutoff() const { return Tcutoff_; }

    /**
     * @brief Initialize internal variables of the mccore object
     *
     * The init() function performs the following:
     *  - calls the init() function of the underlying target and ion_beam objects
     *  - creates electronic stopping and straggling interpolation tables for all atom/material combinations
     *  - creates an array of lab frame scattering cross-section objects (\ref xs_lab) for all ion combinations
     *  - creates tables of mean free path and max impact parameter for flight path selection
     *  - allocates tally memory
     *
     * @return 0 if succesfull
     */
    int init();

    /**
     * @brief Seed the random number generator
     * @param s the seed value
     */
    void seed(unsigned int s) { rng.seed(s); }

    /**
     * @brief Merge the results from another simulation object
     *
     * The tallies and streamed events from @p other are
     * added to this object's tallies and event streams
     *
     * @param other a simulation object
     */
    void merge(const mccore& other);

    /**
     * @brief Run the specified number of ion histories
     *
     * This function runs the actual Monte-Carlo loop n times, repeating
     * the following steps:
     *  - Call the source object to generate new ion
     *  - Call transport() to propagate the ion through target
     *  - Call transport() for all recoils generated by the ion and all their secondary recoils
     *
     * @param n the number of histories to run
     * @return 0 if succesfull
     *
     * @todo Convert run() to run(uint n)
     */
    int run();

protected:

    /**
     * @brief flightPath selects the ion's flight path and impact parameter
     *
     * The selection is performed according to the specified algorithm by
     * mccore::flight_path_type_t.
     * 
     * If mccore::flight_path_type_t is equal to \ref AtomicSpacing or \ref Constant
     * then the flight path \f$\ell\f$ is precalculated and equal to either the material's atomic radius
     * \f$ R_{at} = \left(\frac{3}{4\pi}N\right)^{1/3} \f$ or 
     * to \ref parameters::flight_path_const, respectively.
     * 
     * In both of these cases the impact parameter \f$p\f$ is calculated as
     * $$
     * p = p_{max}\sqrt{u}, \quad p_{max} = (\pi \ell N)^{-1/2}
     * $$
     * where \f$u \in (0,1)\f$ is a random number.
     * 
     * The \ref MendenhallWeller algorithm is as follows:
     * - Get precomputed tables based on ion energy: \f$\ell = \ell(E)\f$, \f$p_{max} = p_{max}(E)\f$
     * - If \f$p_{max}<(\pi R_{at} N)^{-1/2}\f$
     * $$
     * p = p_{max}\sqrt{-\log{u}}
     * $$
     * and reject collision if \f$p>p_{max}\f$.
     * - otherwise \f$\ell = R_{at}\f$ and:
     * $$
     * p = p_{max}\sqrt{u}, \quad p_{max} = (\pi \ell N)^{-1/2}
     * $$
     * 
     * Finally, the \ref IPP algorithm does the following:
     * - Get precomputed tables based on ion energy: \f$\ell_0 = \ell_0(E)\f$, \f$p_{max} = p_{max}(E)\f$
     * - \f$\ell = -\ell_0 \log{u_1}\f$
     * - \f$p = p_{max}\sqrt{u_2}\f$
     * - Reject collision if \f$\ell > \ell_{max}\f$
     * 
     * 
     * @param[in] i pointer to the moving ion object
     * @param[in] m pointer to the material the ion is in
     * @param[out] fp  the selected flight path [nm]
     * @param[out] ip the selected impact parameter [nm]
     * @param[out]  sqrtfp the square root of the selected flight path
     * @return true is the ion should collide at the end of its flight path
     * 
     * @sa \ref flightpath
     */
    bool flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp,
                    std::array<const float *, 3>& fp_par_tbl);

    /**
     * @brief Transport an ion through the target
     *
     * This function propagates an ion by repeating the following steps in a loop:
     *   - Call flightPath() to select flight path and impact parameter
     *   - Call doDedx() to subtract electronic stopping
     *   - Select collision partner and call xs_lab::scatter
     *   - Generate a recoil if conditions allow by calling new_recoil()
     *
     * The loop is repeated until either the ion exits the target or
     * its energy becomes less than parameters::min_energy.
     *
     * Events (recoil generation, ion exit, etc) and energy changes (ionization, phonons, etc)
     *  are registered in the \ref tally object \p t.
     *
     * If \p ev is not null then this ion is part of a recoil cascade
     * and relevant information (e.g., damage energy, vacancies, etc.) if collected in
     * the \ref pka_event object.
     *
     * @param i pointer to the ion object
     * @param t reference to the tally object
     * @param ev pointer to a pka_event object
     * @return 0 if succesfull
     */
    int transport(ion* i, tally& t, pka_event* ev = nullptr);

    /**
     * @brief Generate a new recoil ion
     *
     * This function creates a new recoil ion of atomic species \p target
     * with energy \p recoil_erg, it then calculates its direction of motion
     * based on the initial direction of the scattered ion, \p dir0.
     *
     * The \ref tally object \p t is needed to register the lattice binding energy
     * (if the recoiling atom specifies a non-zero) that is released to phonons.
     *
     * If \p pka is not null then this recoil is part of a cascade
     * and the lattice binding energy is added to the damage energy of
     * the \ref pka_event object.
     *
     * The newly created recoil is put on the respective \ref ion_queue in order to
     * be transported later.
     *
     * @param proj the scattered projectile ion
     * @param target pointer to the atomic species of the target atom
     * @param recoil_erg recoil energy of the target atom
     * @param dir0 initial projectile direction
     * @param mass_ratio mass ratio projectile/target
     * @param t reference to a tally object
     * @param pka pointer to a pka event object
     * @return a pointer to the newly created moving ion
     */
    ion *new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& sqrt_mass_ratio);

    /**
     * @brief Calculate electronic energy loss and straggling for the moving ion
     *
     * The specific stopping & straggling values are obtained by
     * interpolation from the relevant tables by calling interp_dedx()
     *
     * Tables follow the corteo scheme, see \ref dedx_index.
     *
     * @param i pointer to the moving ion
     * @param m pointer to the material the ion is moving in
     * @param fp flight path [nm]
     * @param sqrtfp sqrt of the flight path (used for straggling)
     * @param stopping_tbl pointer to the relevant stopping table [eV/nm]
     * @param straggling_tbl pointer to the relevant straggling table [eV/sqrt{nm}]
     * @return the energy loss [eV]
     */
    float calcDedx(const ion* i, const material* m,
                   float fp, float sqrtfp,
                   const dedx_interp* stopping_tbl,
                   const straggling_interp* straggling_tbl);

    /// Get pointers to electronic loss and straggling tables for a specific ion/material combination
    int getDEtables(const atom* z1, const material* m,
                    const dedx_interp *&dedx,
                    const straggling_interp *&de_stragg) const;

    int getMFPtables(const atom* z1, const material* m,
                    float& fp, float& sqrt_fp,
                    std::array<const float *, 3>& fp_par_tbl) const;


    void pka_end(const ion* i, tally &t, const pka_event &pka);

};

inline int mccore::getDEtables(const atom* z1, const material* m,
                               const dedx_interp *&dedx,
                               const straggling_interp *&de_stragg) const
{
    assert(z1);
    assert(m);
    int ia = z1->id();
    int im = m->id();
    dedx = dedx_(ia,im);
    de_stragg = de_strag_(ia,im);
    return 0;
}

inline int mccore::getMFPtables(const atom *z1, const material *m,
                                float &fp, float &sqrt_fp,
                                std::array<const float *, 3>& fp_par_tbl) const
{
    assert(z1);
    assert(m);

    const float* & ipmax_tbl = fp_par_tbl[0];
    const float* & mfp_tbl   = fp_par_tbl[1];
    const float* & fpmax_tbl = fp_par_tbl[2];

    switch (par_.flight_path_type)
    {
    case AtomicSpacing:
        fp = m->atomicRadius();
        sqrt_fp = 1.f;
        ipmax_tbl = &(ip0[m->id()]);
        break;
    case Constant:
        fp = par_.flight_path_const;
        sqrt_fp = sqrtfp_const[m->id()];
        ipmax_tbl = &(ip0[m->id()]);
        break;
    case MendenhallWeller:
    case IPP:
        ipmax_tbl = &(ipmax_(z1->id(),m->id(),0));
        mfp_tbl = &(mfp_(z1->id(),m->id(),0));
        fpmax_tbl = &(fp_max_(z1->id(),m->id(),0));
        break;
    default:
        fp = sqrt_fp = 0.f;
        assert(false); // should never reach here
    }

    return 0;
}

#endif // SIMULATION_H
