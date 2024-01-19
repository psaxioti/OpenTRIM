#ifndef SETTINGS_H
#define SETTINGS_H

#include "simulation.h"
#include "target.h"

struct settings
{
    struct material_desc {
        std::string name;
        float density;
        bool isMassDensity;
        // std::vector< atom::parameters > atoms;
        std::vector<int> Z;
        std::vector<float> M;
        std::vector<float> X;
        std::vector<float> Ed;
        std::vector<float> El;
        std::vector<float> Es;
        std::vector<float> Er;
        //material_desc& operator=(const material_desc& m);
        //material_desc(const material_desc& m);
        //material_desc() {}
    };

    struct region_desc {
        std::string name;
        int material_id;
        std::vector<float> extX;
        std::vector<float> extY;
        std::vector<float> extZ;
        //region_desc& operator=(const region_desc& m);
    };

    simulation_base::parameters psim_;
    ion_beam::parameters psrc_;
    std::vector< material_desc > materials;
    std::vector< region_desc > regions;
    ivector3 cell_count;
    ivector3 periodic_bc;
    vector3 cell_size;


    settings();

    int parse(const char* fname, bool verbose = false);
    const char* ini_template() const;
    simulation_base* createSimulation() const;

//    settings& operator=(const settings& s) {
//        psim_ = s.psim_; psrc_ = s.psrc_;
//        materials = s.materials; regions = s.regions;
//        cell_count = s.cell_count;
//        cell_size = s.cell_size;
//        periodic_bc = s.periodic_bc;
//        return *this;
//    }
};

struct controller
{
    struct worker {
        settings s_;
        simulation_base* S_;
        uint max_no_ions;
        uint seed;
        int run();
        //worker() {}
//        worker(const worker& w) :
//            s_(w.s_), max_no_ions(w.max_no_ions),
//            seed(w.seed)
//        {}
    };
    double ips;

    int run_simulation(const settings& s);
};

#endif // SETTINGS_H
