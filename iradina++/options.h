#ifndef OPTIONS_H
#define OPTIONS_H

#include "simulation.h"

#include <map>

struct options
{
    struct material_desc {
        std::string name;
        float density{1.f};
        bool isMassDensity{true};
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
        std::string material_id;
        std::vector<float> extX;
        std::vector<float> extY;
        std::vector<float> extZ;
    };

    struct target_desc_t {
        std::vector< std::string > materials;
        std::vector< std::string > regions;
        ivector3 cell_count{1, 1, 1};
        ivector3 periodic_bc{0, 1, 1};
        vector3 cell_size{100.f, 100.f, 100.f};
    };

    simulation_base::parameters sim_par;
    simulation_base::output_options out_opt;
    ion_beam::parameters src_par;
    target_desc_t target_desc;
    std::vector< material_desc > materials;
    std::vector< region_desc > regions;
    std::map<std::string, int> mat2idx;

    options();

    int fromJSON(std::istream& js);

    int parse(std::istream& is, bool verbose = false);
    void print(std::ostream& os);
    const char* ini_template() const;
    int validate();
    simulation_base* createSimulation() const;

};

#endif // OPTIONS_H
