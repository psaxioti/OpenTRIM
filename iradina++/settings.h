#ifndef SETTINGS_H
#define SETTINGS_H

#include "simulation.h"

struct settings
{
    struct material_desc {
        std::string name;
        float density;
        bool isMassDensity;
        std::vector<int> Z;
        std::vector<float> M;
        std::vector<float> X;
        std::vector<float> Ed;
        std::vector<float> El;
        std::vector<float> Es;
        std::vector<float> Er;
    };

    struct region_desc {
        std::string name;
        int material_id;
        std::vector<float> extX;
        std::vector<float> extY;
        std::vector<float> extZ;
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

};

struct controller
{
    double ips;
    int run_simulation(const settings& s);
};

#endif // SETTINGS_H
