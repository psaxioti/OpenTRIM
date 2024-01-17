#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "random_vars.h"
#include "ion_beam.h"
#include "ion_queues.h"
#include "ion.h"
#include "xs.h"
#include "tally.h"



class out_file;

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
        IonsOnly = 3,       /**< Ions only. Recoils are not followed, No damage profiles stored (--- not deeply TESTED!) */
        KP1 = 4,            /**< KP quick calculation of damage, mono-elemental formula similar to SRIM (added by J.-P. Crocombette) */
        KP2 = 5             /**< KP quick calculation of damage, material averaging, more physical as 4 but different compared to SRIM (added by J.-P. Crocombette) */
    } simulation_type_t;

    /**
     * @brief Determines how ion scattering is simulated
     */
    typedef enum {
        Corteo4bit = 0,
        Corteo6bit = 1,
        ZBL_MAGICK = 2
    } scattering_calculation_t;

    /**
     * @brief Flight path selection algorithm
     *
     * Determines the algorithm to select the
     * free flight path \f$\ell\f$ between collisions
     */
    typedef enum {
        Poisson = 0,    /**< Poisson distributed \f$P(\ell) = e^{-\ell/ell_0}\f$*/
        AtomicSpacing,  /**< Constant, equal to interatomic distance */
        Constant,       /**< Constant, equal to user supplied value */
        SRIMlike        /**< Algorithm similar to SRIM */
    } flight_path_type_t;

    /**
     * @brief The model used to calculate ion straggling
     */
    typedef enum {
        NoStraggling = 0,
        BohrStraggling,
        ChuStraggling,
        YangStraggling
    } straggling_model_t;

    /**
     * @brief Type of random variable sampling
     */
    typedef enum {
        Sampled = 0,    /**< Sampling directly from distributions */
        Tabulated = 1   /**< Sampling from tabulated values */
    } random_var_t;

    /**
     * @brief The random number generator used
     */
    typedef enum {
        MinStd = 0,  /**< Minimum standard - fast but low statistical quality */
        MersenneTwister = 1   /**< slower but better statistics */
    } random_generator_t;

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
        unsigned int max_no_ions;
        simulation_type_t simulation_type;
        scattering_calculation_t scattering_calculation;
        flight_path_type_t flight_path_type;
        straggling_model_t straggling_model;
        float flight_path_const;
        float min_energy;
        random_var_t random_var_type;
        random_generator_t random_generator_type;
        int threads;
        parameters(); // set defaults
    };

protected:

    ion_queues< ion > q_;

    // simulation paramenters
    parameters par_;
    // components
    std::shared_ptr<ion_beam> source_;
    std::shared_ptr<target> target_;
    tally tally_;

    // helper variable for flight path calc
    std::vector<float> sqrtfp_const;

    // Stopping & Straggling Tables
    Array3Df dedx_; // stopping data (atoms x materials x energy)
    Array2Df dedx1; // proton stopping (materials x energy)
    Array3Df de_strag_; // straggling data (atoms x materials x energy)

    // timing
    double ms_per_ion_;

    friend class out_file;

public:

    virtual ~simulation_base();

    static simulation_base* fromParameters(const parameters &par);

    void setIonBeam(const ion_beam::parameters& p) {
        source_->setParameters(p);
    }

    void setMaxIons(unsigned int n) { par_.max_no_ions = n; }

    const parameters& getParameters() const { return par_; }
    const std::string& title() const { return par_.title; }
    simulation_type_t simulationType() const { return par_.simulation_type; }
    flight_path_type_t flightPathType() const { return par_.flight_path_type; }
    straggling_model_t stragglingModel() const { return par_.straggling_model; }
    random_var_t randomVarType() const { return par_.random_var_type; }
    random_generator_t rngType() const { return par_.random_generator_type; }

    double ms_per_ion() const { return ms_per_ion_; }

    // simulation setup
    material* addMaterial(const char* name, const float& density) {
        return target_->addMaterial(name, density);
    }
    //void setProjectile(int Z, float M, float E0);
    grid3D& grid() { return target_->grid(); }
    const grid3D& grid() const { return target_->grid(); }
    void fill(const box3D& box, const material* m) {
        target_->fill(box,m);
    }

    const target* getTarget() const { return target_.get(); }
    const ion_beam* getSource() const { return source_.get(); }
    const tally& getTally() const { return tally_; }


    int saveTallys();

    virtual int init();
    virtual int run() = 0;
    int run(int nthreads);

protected:

    ion* new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio);

    void tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed);

    float LSS_Tdam(int Z, float M, float T);
    float NRT(float Ed, float Tdam);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;

    // protected constructor
    // cannot instantiate simulation base objects
    simulation_base(const parameters& p);

    virtual simulation_base* clone();

    virtual int transport(ion* i) = 0;
    virtual float flightPath(const ion* i, const material* m, float& sqrtfp, float &pmax) = 0;
    virtual float impactPar(const ion* i, const material* m, const float &sqrtfp, const float &pmax) = 0;

};

template<class _XScm, class _RNG_E>
class simulation : public simulation_base
{
public:

    typedef _XScm reducedXScm;
    typedef XSlab<_XScm> scatteringXSlab;

    typedef _RNG_E rng_engine;
    typedef URBG_< _RNG_E > URBG;

private:

    typedef simulation< _XScm, _RNG_E > _Myt;

    URBG urbg;
    random_vars_base* rnd;
    Array2D< std::shared_ptr<scatteringXSlab> > scattering_matrix_;

public:
    simulation(const simulation_base::parameters& p);
    ~simulation();

    virtual int init() override;

    virtual int run() override;

protected:
    virtual int transport(ion* i) override;
    virtual float flightPath(const ion* i, const material* m, float& sqrtfp, float &pmax) override;
    virtual float impactPar(const ion* i, const material* m, const float &sqrtfp, const float &pmax) override;

    random_vars_base* createRandomVars();

    virtual simulation_base* clone() override;

};

typedef simulation<XS_zbl_magic,  std::mt19937> SimZBLMagic_MT;
typedef simulation<XS_corteo4bit, std::mt19937> SimCorteo4bit_MT;
// typedef simulation<XS_corteo6bit, std::mt19937> SimCorteo6bit_MT;

typedef simulation<XS_zbl_magic,  std::minstd_rand> SimZBLMagic_MSRAND;
typedef simulation<XS_corteo4bit, std::minstd_rand> SimCorteo4bit_MSRAND;
// typedef simulation<XS_corteo6bit, std::minstd_rand> SimCorteo6bit_MSRAND;

#endif // SIMULATION_H
