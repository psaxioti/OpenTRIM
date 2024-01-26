#ifndef SETTINGS_H
#define SETTINGS_H

#include "simulation.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct settings;

void to_json(json& j, const settings& p);
void from_json(const json& j, settings& p);

struct settings
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
        int material_id{0};
        std::vector<float> extX;
        std::vector<float> extY;
        std::vector<float> extZ;
    };

    struct target_desc {
        std::vector< std::string > materials;
        std::vector< std::string > regions;
        ivector3 cell_count{1, 1, 1};
        ivector3 periodic_bc{0, 1, 1};
        vector3 cell_size{100.f, 100.f, 100.f};
    };

    simulation_base::parameters psim_;
    simulation_base::output_options out_opt_;
    ion_beam::parameters psrc_;
    target_desc target_desc_;
    std::vector< material_desc > materials;
    std::vector< region_desc > regions;

    ivector3 cell_count{1, 1, 1};
    ivector3 periodic_bc{0, 1, 1};
    vector3 cell_size{100.f, 100.f, 100.f};



    settings();

    static int test();

    int parse(json j);
    int parse(std::istream& is, bool verbose = false);
    void print(std::ostream& os);
    const char* ini_template() const;
    simulation_base* createSimulation() const;

};

#endif // SETTINGS_H
