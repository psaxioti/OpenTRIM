#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "random_vars.h"
#include "ion_beam.h"
#include "ion_queues.h"
#include "ion.h"
#include "xs.h"
#include "tally.h"
#include "event_stream.h"

#include <atomic>

class out_file;
class options;

/**
 * \defgroup MC Monte-Carlo simulation components
 * @{
 *
 * @}
 */

/**
 * @brief The simulation_base class forms the basis of all simulation classes
 *
 * @ingroup MC
 */
class simulation_base
{
public:

    /**
     * @brief Type of simulation
     */
    typedef enum {
        FullCascade = 0,    /**< Full Damage Cascade, follow recoils */
        Invalid1 = 1,       // Compatibility. Iradina does not have simulation_type=1 & 2
        Invalid2 = 2,
        IonsOnly = 3,       /**< Ions only. Recoils are not followed. Damage estimate by NRT */
        Cascades = 4,       /**< Cascades generated and followed */
        InvalidSimulationType = -1
    } simulation_type_t;

    /**
     * @brief Detail of NRT implementation
     */
    typedef enum {
        NRT_element = 0, /**< NRT calculated for using Ed of struck atom (similar to SRIM) */
        NRT_average = 1,  /**< NRT calculated with average Ed (J.-P. Crocombette 2019) */
        NRT_InvalidOption = -1
    } nrt_calculation_t;

    /**
     * @brief Determines how ion scattering is simulated
     */
    typedef enum {
        Corteo4bit = 0,
        Corteo6bit = 1,
        ZBL_MAGICK = 2,
        InvalidScatteringOption = -1
    } scattering_calculation_t;

    /**
     * @brief Flight path selection algorithm
     *
     * Determines the algorithm to select the
     * free flight path \f$\ell\f$ between collisions
     */
    typedef enum {
        AtomicSpacing = 1,  /**< Constant, equal to interatomic distance */
        Constant = 2,       /**< Constant, equal to user supplied value */
        SRIMlike = 3,        /**< Algorithm similar to SRIM */
        MendenhallWeller = 4,
        MyFFP = 5,
        InvalidPath = -1
    } flight_path_type_t;

    /**
     * @brief The model used to calculate ion straggling
     */
    typedef enum {
        NoStraggling = 0,
        BohrStraggling,
        ChuStraggling,
        YangStraggling,
        InvalidStraggling = -1
    } straggling_model_t;

    /**
     * @brief Type of simulation event
     */
    typedef enum {
        NewSourceIon = 0,
        NewRecoil,
        Scattering,
        IonExit,
        IonStop
    } simulation_event_t;

    struct parameters {
        std::string title;
        unsigned int max_no_ions{100};
        simulation_type_t simulation_type{FullCascade};
        nrt_calculation_t nrt_calculation{NRT_element};
        scattering_calculation_t scattering_calculation{Corteo4bit};
        flight_path_type_t flight_path_type{AtomicSpacing};
        straggling_model_t straggling_model{YangStraggling};
        float flight_path_const{0.1};
        float min_energy{1.f};
        int threads{1};
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

    ion_queues< ion > q_;

    // simulation paramenters
    parameters par_;
    output_options out_opts_;

    // components
    ion_beam* source_;
    target* target_;
    tally tally_;
    event_stream pka_stream_, exit_stream_;
    pka_event pka;
    exit_event exit_ev;

    // ref counter
    std::shared_ptr<int> ref_count_;
    int thread_id_;
    uint count_offset_;

    // helper variable for flight path calc
    std::vector<float> sqrtfp_const;

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

    virtual ~simulation_base();

    static simulation_base* fromParameters(const parameters &par);

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

    // simulation setup
    material* addMaterial(const char* name)
    { return target_->addMaterial(name); }

    //void setProjectile(int Z, float M, float E0);
    grid3D& grid() { return target_->grid(); }
    const grid3D& grid() const { return target_->grid(); }
    void fill(const box3D& box, const material* m) {
        target_->fill(box,m);
    }

    const target* getTarget() const { return target_; }
    const ion_beam* getSource() const { return source_; }
    const tally& getTally() const { return tally_; }
    void addTally(const tally& t) { tally_ += t; }


    int saveTallys();
    virtual int init();
    int exec(progress_callback cb = 0, uint msInterval = 1000);
    virtual void seed(unsigned int s) = 0;

protected:

    void new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio, pka_event* pka);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;

    // protected constructor
    // cannot instantiate simulation base objects
    simulation_base();
    simulation_base(const parameters& p);
    simulation_base(const simulation_base& s);

    std::string outFileName(const char* type);

    virtual int transport(ion* i, pka_event* ev = nullptr) = 0;
    virtual int flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp) = 0;
    virtual int run() = 0;
    virtual simulation_base* clone() const = 0;
    uint ions_done() { return nion_thread_.exchange(0); }

};

template<class _XScm>
class simulation : public simulation_base
{
public:

    typedef _XScm reducedXScm;
    typedef xs_lab<_XScm> scatteringXSlab;

private:

    typedef simulation< _XScm > _Myt;

    random_vars rng;
    Array2D<scatteringXSlab*> scattering_matrix_;

public:
    explicit simulation(const char* t = 0);
    explicit simulation(const parameters& p);
    simulation(const _Myt& S);
    ~simulation();

    virtual int init() override;

protected:
    virtual int transport(ion* i, pka_event* ev = nullptr) override;
    virtual int flightPath(const ion* i, const material* m, float& fp, float& ip, float& sqrtfp) override;
    void doDedx(ion* i, const material* m, float fp, float sqrtfp, const float* stopping_tbl, const float* straggling_tbl);
    virtual int run() override;
    virtual simulation_base* clone() const override { return new _Myt(*this); }
    virtual void seed(unsigned int s) override { rng.seed(s); }
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


typedef simulation<xs_zbl_magic> SimZBLMagic;
typedef simulation<xs_corteo4bit> SimCorteo4bit;
typedef simulation<xs_corteo6bit> SimCorteo6bit;


#endif // SIMULATION_H
