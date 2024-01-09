#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "target.h"
#include "urbg.h"
#include "ion_beam.h"
#include "ion.h"
#include "xs.h"

#include <queue>

struct random_vars;
struct reducedXS_base;
class out_file;

class simulation
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

private:

    typedef std::queue<ion*> ion_queue_t;
    ion_queue_t ion_buffer_;
    ion_queue_t ion_queue_;

    std::string name_;
    URBG urbg;
    target target_;
    inventory inventory_;
    ion_beam source_;
    reducedXS xs_;
    random_vars* rnd;

    flight_path_type_t flight_path_type;
    simulation_type_t simulation_type;
    StragglingModel straggling_model;

    float flight_path_const;
    std::vector<float> sqrtfp_const;

    float energy_cutoff_;

    double ms_per_ion_;

    unsigned int ion_histories_;

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
    simulation(const char* name, const reducedXS& x = reducedXS());
    ~simulation();

    const std::string& name() { return name_; }
    void setName(const char* n) { name_ = n; }
    const reducedXS& xs() { return xs_; }
    void setXS(const reducedXS& x) { xs_ = x; }

    void setStragglingModel(StragglingModel m) { straggling_model = m; }

    unsigned int ion_histories() const { return ion_histories_; }
    double ms_per_ion() const { return ms_per_ion_; }

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

    int init();

    void tallyIonizationEnergy(const ion* i, const float& v);
    void tallyPhononEnergy(const ion* i, const float& v);
    void tallyKP(const ion* i, const float& Tdam, const float& nv, const float& Ed);

    int run(int count, const char* outfname);

    ion* new_ion(const ion* parent = nullptr);
    ion* new_recoil(const ion* proj, const atom* target, const float& recoil_erg,
                    const vector3& dir0, const float& mass_ratio);
    void free_ion(ion* i);
    ion* pop_ion();
    void push_ion(ion* i);

    int transport(ion* i);
    float flightPath(const ion* i, const material* m, float& sqrtfp, float &pmax);
    float impactPar(const ion* i, const material* m, const float &sqrtfp, const float &pmax);

    float LSS_Tdam(int Z, float M, float T);
    float NRT(float Ed, float Tdam);
};

#endif // SIMULATION_H
