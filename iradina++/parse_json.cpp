#include "options.h"

#include <fstream>
#include <iostream>

#define JSON_DIAGNOSTICS 1
#include <nlohmann/json.hpp>

using std::cout;
using std::cerr;
using std::endl;

// Define a special json type for:
//   - keeping the order of elements
//   - use float for real numbers (to avoid e.g. 0.100001, etc)
using ojson = nlohmann::basic_json<nlohmann::ordered_map,
                         std::vector,
                         std::string,
                         bool,
                         std::int64_t,
                         std::uint64_t,
                         float>;


// json serialization of Eigen::AlignedVector3<T> types
namespace nlohmann {

template <typename T>
struct adl_serializer< Eigen::AlignedVector3<T> > {
    static void to_json(ojson& j, const Eigen::AlignedVector3<T>& V3) {
        j = {V3.x(), V3.y(), V3.z()};
    }

    static void from_json(const ojson& j, Eigen::AlignedVector3<T>& V3) {
        std::array<T,3> v = j.template get< std::array<T,3> >();
        V3 = Eigen::AlignedVector3<T>( v[0], v[1], v[2] );
    }
};

} // end namespace nlohmann

// define my macro for ojson
#define MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(Type, ...)  \
    inline void to_json(ojson& nlohmann_json_j, const Type& nlohmann_json_t) { NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_TO, __VA_ARGS__)) } \
    inline void from_json(const ojson& nlohmann_json_j, Type& nlohmann_json_t) { const Type nlohmann_json_default_obj{}; NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_FROM_WITH_DEFAULT, __VA_ARGS__)) }

// serialization of options struct
void to_json(ojson& j, const options& p);
void from_json(const ojson& j, options& p);

int options::parseJSON(std::istream& js)
{
    try {
        ojson j = ojson::parse(js,nullptr,true,true);
        *this = j.template get<options>();
        validate();
    }
    catch (const ojson::exception& e) {
        cerr << "Error reading json input:" << endl;
        cerr << e.what() << std::endl;
        return -1;
    }
    catch (const std::invalid_argument& e) {
        cerr << "Invalid option:" << endl;
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(target::region_desc_t,
                                   material_id,
                                   min, max)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(material::material_desc_t,
                                   density, isMassDensity,
                                   Z, M, X, Ed, El, Es, Er)

NLOHMANN_JSON_SERIALIZE_ENUM(ion_beam::ion_distribution_t, {
                                {ion_beam::InvalidIonDistribution, nullptr},
                                {ion_beam::SurfaceRandom, "SurfaceRandom"},
                                {ion_beam::SurfaceCentered, "SurfaceCentered"},
                                {ion_beam::FixedPos, "FixedPos"},
                                {ion_beam::VolumeCentered, "VolumeCentered"},
                                {ion_beam::VolumeRandom, "VolumeRandom"},
                                {ion_beam::SurfaceRandom, 0},
                                {ion_beam::SurfaceCentered, 1},
                                {ion_beam::FixedPos, 2},
                                {ion_beam::VolumeCentered, 3},
                                {ion_beam::VolumeRandom, 4}
                            })

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(ion_beam::parameters,
                                   ion_distribution,
                                   ionZ, ionM, ionE0,
                                   dir, pos)


NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::simulation_type_t, {
                                {simulation_base::InvalidSimulationType, nullptr},
                                {simulation_base::FullCascade, "FullCascade"},
                                {simulation_base::IonsOnly, "IonsOnly"},
                                {simulation_base::FullCascade, 0},
                                {simulation_base::FullCascade, 1},
                                {simulation_base::FullCascade, 2},
                                {simulation_base::IonsOnly, 3}
                            })

NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::nrt_calculation_t, {
                                {simulation_base::NRT_InvalidOption, nullptr},
                                {simulation_base::NRT_element, "NRT_element"},
                                {simulation_base::NRT_average, "NRT_average"},
                                {simulation_base::NRT_element, 0},
                                {simulation_base::NRT_average, 1}
                            })


NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::flight_path_type_t, {
                                {simulation_base::InvalidPath, nullptr},
                                {simulation_base::AtomicSpacing, "AtomicSpacing"},
                                {simulation_base::Constant, "Constant"},
                                {simulation_base::SRIMlike, "SRIMlike"},
                                {simulation_base::MendenhallWeller, "MendenhallWeller"},
                                {simulation_base::MyFFP, "MyFFP"},
                                {simulation_base::AtomicSpacing, 1},
                                {simulation_base::Constant, 2},
                                {simulation_base::SRIMlike, 3},
                                {simulation_base::MendenhallWeller, 4},
                                {simulation_base::MyFFP, 5}
                            })   

NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::scattering_calculation_t, {
                                {simulation_base::InvalidScatteringOption, nullptr},
                                {simulation_base::Corteo4bit, "Corteo4bit"},
                                {simulation_base::Corteo6bit, "Corteo6bit"},
                                {simulation_base::ZBL_MAGICK, "ZBL_MAGICK"},
                                {simulation_base::Corteo4bit, 0},
                                {simulation_base::Corteo6bit, 1},
                                {simulation_base::ZBL_MAGICK, 2}
                            })    

NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::straggling_model_t, {
                                {simulation_base::InvalidStraggling, nullptr},
                                {simulation_base::NoStraggling, "NoStraggling"},
                                {simulation_base::BohrStraggling, "BohrStraggling"},
                                {simulation_base::ChuStraggling, "ChuStraggling"},
                                {simulation_base::YangStraggling, "YangStraggling"},
                                {simulation_base::NoStraggling, 0},
                                {simulation_base::BohrStraggling, 1},
                                {simulation_base::ChuStraggling, 2},
                                {simulation_base::YangStraggling, 3}
                            })                                                                                        

NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::random_generator_t, {
                                {simulation_base::InvalidRandomGenerator, nullptr},
                                {simulation_base::MersenneTwister, "MersenneTwister"},
                                {simulation_base::MinStd, "MinStd"},
                                {simulation_base::Xoshiro128p, "Xoshiro128p"},
                                {simulation_base::MersenneTwister, 0},
                                {simulation_base::MinStd, 1},
                                {simulation_base::Xoshiro128p, 2}
                            })

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::parameters,
                                   title, max_no_ions, simulation_type,
                                   nrt_calculation, scattering_calculation,
                                   flight_path_type, straggling_model,
                                   flight_path_const, min_energy,
                                   random_generator_type,
                                   threads, seeds)


MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::output_options,
                                                OutputFileBaseName,
                                                storage_interval,
                                                store_transmitted_ions,
                                                store_range_3d,
                                                store_ion_paths,
                                                store_path_limit,
                                                store_recoil_cascades,
                                                store_path_limit_recoils,
                                                store_pka,
                                                store_dedx
                                                )

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(target::target_desc_t,
                                                materials,
                                                regions,
                                                cell_count,
                                                periodic_bc,
                                                cell_size)

void to_json(ojson& j, const options& p)
{
    j["Simulation"] = p.Simulation;
    j["IonBeam"] = p.IonBeam;
    j["Target"] = p.Target;
    j["Output"] = p.Output;
    for(int i=0; i<p.materials_desc.size(); i++) {
        const material::material_desc_t& md = p.materials_desc[i];
        const std::string& name = p.Target.materials[i];
        j[name]=md;
    }
    for(int i=0; i<p.regions_desc.size(); i++) {
        const target::region_desc_t& rd = p.regions_desc[i];
        const std::string& name = p.Target.regions[i];
        j[name]=rd;
    }
}
void from_json(const ojson& j, options& p)
{
    p = options();

    if (j.contains("Simulation"))
        p.Simulation = j["Simulation"];
    if (j.contains("Output"))
        p.Output = j["Output"];
    if (j.contains(("IonBeam")))
        p.IonBeam = j["IonBeam"];
    if (!j.contains("Target")) {
        throw std::invalid_argument("Required section \"Target\" not found.");
    }
    p.Target = j["Target"];

    if (p.Target.materials.empty()) {
        throw std::invalid_argument("No materials specified in \"Target\".");
    }
    if (p.Target.regions.empty()) {
        throw std::invalid_argument("No regions specified in \"Target\".");
    }

    for(const std::string& m : p.Target.materials)
    {
        if (j.contains(m)) {
            auto jm = j[m];
            material::material_desc_t md = jm.template get< material::material_desc_t >();
            p.materials_desc.push_back(md);
        } else {
            std::string msg;
            msg = "Definition of material \"";
            msg += m;
            msg += "\" not found";
            throw std::invalid_argument(msg);
        }
    }
    for(const std::string& r : p.Target.regions)
    {
        if (j.contains(r)) {
            auto jr = j[r];
            target::region_desc_t rd = jr.template get< target::region_desc_t >();
            if (p.materialIdx(rd.material_id)<0) {
                std::string msg;
                msg = "In region \"";
                msg += r;
                msg += "\" the specified material_id=\"";
                msg += rd.material_id;
                msg += "\" is not among the target materials.";
                throw std::invalid_argument(msg);
            }
            p.regions_desc.push_back(rd);
        } else {
            std::string msg;
            msg = "Definition of region \"";
            msg += r;
            msg += "\" not found";
            throw std::invalid_argument(msg);
        }
    }
}

void options::printJSON(std::ostream& os)
{
    ojson j(*this);
    os << j.dump(4) << endl;
}

