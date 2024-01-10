#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "urbg.h"
#include "ion_beam.h"
#include "ion.h"
#include "xs.h"

#include <queue>

struct random_vars;
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

    typedef std::queue<ion*> ion_queue_t;
    ion_queue_t ion_buffer_;
    ion_queue_t ion_queue_;

    std::string name_;

    URBG urbg;
    target target_;
    inventory inventory_;
    ion_beam source_;
    random_vars* rnd;

    flight_path_type_t flight_path_type;
    simulation_type_t simulation_type;
    straggling_model_t straggling_model;

    float flight_path_const;
    std::vector<float> sqrtfp_const;

    float energy_cutoff_;

    double ms_per_ion_;

    unsigned int ion_histories_;

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

public:
    simulation_base(const char* name);
    ~simulation_base();

    const std::string& name() { return name_; }
    void setName(const char* n) { name_ = n; }

    void setStragglingModel(straggling_model_t m) { straggling_model = m; }

    unsigned int ion_histories() const { return ion_histories_; }
    double ms_per_ion() const { return ms_per_ion_; }

    material* addMaterial(const char* name, const float& density) {
        return inventory_.addMaterial(name, density);
    }
    void setProjectile(int Z, float M, float E0);

    int getDEtables(const atom* z1, const material* m,
                    const float *&dedx, const float *&de_stragg) const;

    grid3D& grid() { return target_.grid(); }
    const grid3D& grid() const { return target_.grid(); }
    void fill(const box3D& box, const material* m) {
        target_.fill(box,m);
    }

    const target& getTarget() const { return target_; }
    const inventory& getInventory() const { return inventory_; }
    const ion_beam& getSource() const { return source_; }

    virtual int init();

    void tallyIonizationEnergy(const ion* i, const float& v);
    void tallyPhononEnergy(const ion* i, const float& v);
    void tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed);

    ion* new_ion(const ion* parent = nullptr);
    ion* new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio);
    void free_ion(ion* i);
    ion* pop_ion();
    void push_ion(ion* i);

    float flightPath(const ion* i, const material* m, float& sqrtfp, float &pmax);
    float impactPar(const ion* i, const material* m, const float &sqrtfp, const float &pmax);

    float LSS_Tdam(int Z, float M, float T);
    float NRT(float Ed, float Tdam);
};

template<class _XScm>
class simulation : public simulation_base
{
public:

    typedef _XScm reducedXScm;
    typedef XSlab<_XScm> scatteringXSlab;

private:

    Array2D<scatteringXSlab*> scattering_matrix_;

public:
    simulation(const char* name);
    ~simulation();

    virtual int init() override;

    int run(int count, const char* outfname);



    int transport(ion* i);

};

typedef simulation<XS_zbl_magic>  SimZBLMagic;
typedef simulation<XS_corteo4bit> SimCorteo4bit;
typedef simulation<XS_corteo6bit> SimCorteo6bit;

#endif // SIMULATION_H
