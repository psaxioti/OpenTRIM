#ifndef OPTIONS_H
#define OPTIONS_H

#include "simulation.h"

#include <map>

struct options
{
    simulation_base::parameters Simulation;
    simulation_base::output_options Output;
    ion_beam::parameters IonBeam;
    target::target_desc_t Target;
    std::vector< material::material_desc_t > materials_desc;
    std::vector< target::region_desc_t > regions_desc;

    int parseJSON(std::istream& js);
    void printJSON(std::ostream& os);

    int parseINI(std::istream& is, bool verbose = false);
    const char* ini_template() const;

    int validate();
    simulation_base* createSimulation() const;

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

#endif // OPTIONS_H
