#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "random_vars.h"
#include "ion_beam.h"
#include "ion.h"
#include "xs.h"

#include <queue>

class out_file;

class simulation_base
{
public:

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

    typedef enum {
        NoStraggling = 0,
        BohrStraggling,
        ChuStraggling,
        YangStraggling
    } straggling_model_t;

    /**
     * @brief Type of simulation
     */
    typedef enum {
        FullCascade = 0,    /**< Full Damage Cascade, follow recoils */
        IonsOnly = 3,       /**< Ions only. Recoils are not followed, No damage profiles stored (--- not deeply TESTED!) */
        KP1 = 4,            /**< KP quick calculation of damage, mono-elemental formula similar to SRIM (added by J.-P. Crocombette) */
        KP2 = 5             /**< KP quick calculation of damage, material averaging, more physical as 4 but different compared to SRIM (added by J.-P. Crocombette) */
    } simulation_type_t;

    typedef enum {
        Sampled = 0,
        Tabulated
    } random_var_t;

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

protected:
    // FIFO ion buffers
    typedef std::queue<ion*> ion_queue_t;
    ion_queue_t ion_buffer_; // buffer of allocated ion objects
    ion_queue_t recoil_queue_; // queue of generated recoils
    ion_queue_t pka_queue_; // queue of generated PKAs

    // return finished ion to buffer
    void free_ion(ion* i) { ion_buffer_.push(i); }
    // pop an ion from the respective queue
    static ion* pop_one_(ion_queue_t& Q)
    { ion* i = Q.front(); Q.pop(); return i; }
    // ...
    ion* pop_pka() { return pop_one_(pka_queue_); }
    // ...
    ion* pop_recoil() { return pop_one_(recoil_queue_); }
    // push an ion. if recoil_id>0 goes to recoil queue
    void push(ion* i)
    { i->recoil_id > 1 ? recoil_queue_.push(i) : pka_queue_.push(i); }

    std::string name_;

    target target_;
    inventory inventory_;
    ion_beam source_;

    simulation_type_t simulation_type;
    flight_path_type_t flight_path_type;
    straggling_model_t straggling_model;
    random_var_t random_var_type;

    float flight_path_const;
    std::vector<float> sqrtfp_const;

    float energy_cutoff_;

    // statistics
    double ms_per_ion_;
    unsigned int Nions_;
    unsigned int Npkas_;
    unsigned int Nrecoils_;

    Array3Df dedx_; // stopping data (atoms x materials x energy)
    Array2Df dedx1; // proton stopping (materials x energy)
    Array3Df de_strag_; // straggling data (atoms x materials x energy)

    /*
     * Tallys
     */
    Array2Dui InterstitialTally;
    Array2Dui ReplacementTally;
    Array2Dui VacancyTally;
    Array2Dd KPTally;
    Array2Dd IonizationEnergyTally;
    Array2Dd PhononEnergyTally;

    out_file* out_file_;
    friend class out_file;

    ion* new_ion(const ion* parent = nullptr);
    ion* new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio);

    void tallyIonizationEnergy(const ion* i, const float& v);
    void tallyPhononEnergy(const ion* i, const float& v);
    void tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed);

    float LSS_Tdam(int Z, float M, float T);
    float NRT(float Ed, float Tdam);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;




public:
    simulation_base(const char* name);
    ~simulation_base();

    const std::string& name() { return name_; }
    void setName(const char* n) { name_ = n; }

    simulation_type_t simulationTypel() const { return simulation_type; }
    void setSimulationType(simulation_type_t m) { simulation_type = m; }

    flight_path_type_t flightPathType() const { return flight_path_type; }
    void setFlightPathType(flight_path_type_t m) { flight_path_type = m; }

    straggling_model_t stragglingModel() const { return straggling_model; }
    void setStragglingModel(straggling_model_t m) { straggling_model = m; }

    random_var_t randomVarType() const { return random_var_type; }
    void setRandomVarType(random_var_t m) { random_var_type = m; }

    unsigned int ion_histories() const { return Nions_; }
    double ms_per_ion() const { return ms_per_ion_; }
    unsigned int pkas() const { return Npkas_; }
    unsigned int recoils() const { return Nrecoils_; }

    material* addMaterial(const char* name, const float& density) {
        return inventory_.addMaterial(name, density);
    }
    void setProjectile(int Z, float M, float E0);
    grid3D& grid() { return target_.grid(); }
    const grid3D& grid() const { return target_.grid(); }
    void fill(const box3D& box, const material* m) {
        target_.fill(box,m);
    }

    const target& getTarget() const { return target_; }
    const inventory& getInventory() const { return inventory_; }
    const ion_beam& getSource() const { return source_; }

    virtual int init();

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

    URBG urbg;
    random_vars_base* rnd;

    Array2D<scatteringXSlab*> scattering_matrix_;

public:
    simulation(const char* name);
    ~simulation();

    virtual int init() override;

    int run(int count, const char* outfname);

    int transport(ion* i);

    float flightPath(const ion* i, const material* m, float& sqrtfp, float &pmax);
    float impactPar(const ion* i, const material* m, const float &sqrtfp, const float &pmax);

};

typedef simulation<XS_zbl_magic,  std::mt19937> SimZBLMagic_MT;
typedef simulation<XS_corteo4bit, std::mt19937> SimCorteo4bit_MT;
typedef simulation<XS_corteo6bit, std::mt19937> SimCorteo6bit_MT;

typedef simulation<XS_zbl_magic,  std::minstd_rand> SimZBLMagic_MSRAND;
typedef simulation<XS_corteo4bit, std::minstd_rand> SimCorteo4bit_MSRAND;
typedef simulation<XS_corteo6bit, std::minstd_rand> SimCorteo6bit_MSRAND;

#endif // SIMULATION_H
