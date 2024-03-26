#ifndef MCDRIVER_H
#define MCDRIVER_H

#include "simulation.h"

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

    typedef void (*progress_callback)(const mcdriver& v);

    static void def_progress_callback(const mcdriver& d);

protected:

    // timing
    double ips_; // ions/s

    parameters par_;
    output_options out_opts_;

    simulation* s_;

    std::string outFileName(const char* type, int thread_id);

    std::vector<uint> thread_ion_count_;



public:
    mcdriver();
    ~mcdriver();

    void getOptions(options& opt) const;
    void setOptions(const options& o);

    const output_options& outputOptions(const output_options& opts) const
    { return out_opts_; }
    void setOutputOptions(const output_options& opts)
    { out_opts_ = opts; }

    const simulation* getSim() { return s_; }

    double ips() const { return ips_; }
    int nThreads() const { return par_.threads; }

    int saveTallys();

    const std::vector<uint>& thread_ion_count() const { return thread_ion_count_; }

    int exec(progress_callback cb = def_progress_callback, uint msInterval = 1000);

};


struct options
{
    mcdriver::parameters Driver;
    mcdriver::output_options Output;
    simulation::parameters Simulation;
    ion_beam::parameters IonBeam;
    target::target_desc_t Target;
    std::vector< material::material_desc_t > materials_desc;
    std::vector< target::region > regions_desc;

    int parseJSON(std::istream& js);
    void printJSON(std::ostream& os) const;

    int validate();
    simulation* createSimulation() const;

    int regionIdx(const std::string& name) const {
        for(int i=0; i<Target.regions.size(); i++) {
            if (Target.regions[i] == name) {
                return i;
            }
        }
        return -1;
    }
    int materialIdx(const std::string& name) const {
        for(int i=0; i<Target.materials.size(); i++) {
            if (Target.materials[i] == name) {
                return i;
            }
        }
        return -1;
    }

};

#endif // MCDRIVER_H
