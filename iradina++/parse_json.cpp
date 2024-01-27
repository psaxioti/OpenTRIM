#include "settings.h"

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

// using ordered_json = nlohmann::ordered_json;
using std::cout;
using std::cerr;
using std::endl;

using ojson = nlohmann::basic_json<nlohmann::ordered_map,
                         std::vector,
                         std::string,
                         bool,
                         std::int64_t,
                         std::uint64_t,
                         float>;



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

} // end namespace

#define MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(Type, ...)  \
    inline void to_json(ojson& nlohmann_json_j, const Type& nlohmann_json_t) { NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_TO, __VA_ARGS__)) } \
    inline void from_json(const ojson& nlohmann_json_j, Type& nlohmann_json_t) { const Type nlohmann_json_default_obj{}; NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_FROM_WITH_DEFAULT, __VA_ARGS__)) }


void to_json(ojson& j, const settings& p);
void from_json(const ojson& j, settings& p);

int settings::fromJSON(std::istream& js)
{
    try {
        ojson j = ojson::parse(js,nullptr,true,true);
        *this = j.template get<settings>();
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

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::region_desc,
                                   name,
                                   material_id,
                                   extX, extY, extZ)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::material_desc,
                                   name,
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
                                {simulation_base::Poisson, "Poisson"},
                                {simulation_base::AtomicSpacing, "AtomicSpacing"},
                                {simulation_base::Constant, "Constant"},
                                {simulation_base::SRIMlike, "SRIMlike"},
                                {simulation_base::Poisson, 0},
                                {simulation_base::AtomicSpacing, 1},
                                {simulation_base::Constant, 2},
                                {simulation_base::SRIMlike, 3}
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

NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::random_var_t, {
                                {simulation_base::InvalidRandomVar, nullptr},
                                {simulation_base::Sampled, "Sampled"},
                                {simulation_base::Tabulated, "Tabulated"},
                                {simulation_base::Sampled, 0},
                                {simulation_base::Tabulated, 1}
                            })


NLOHMANN_JSON_SERIALIZE_ENUM(simulation_base::random_generator_t, {
                                {simulation_base::InvalidRandomGenerator, nullptr},
                                {simulation_base::MersenneTwister, "MersenneTwister"},
                                {simulation_base::MinStd, "MinStd"},
                                {simulation_base::MersenneTwister, 0},
                                {simulation_base::MinStd, 1}
                            })

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::parameters,
                                   title, max_no_ions, simulation_type,
                                   nrt_calculation, scattering_calculation,
                                   flight_path_type, straggling_model,
                                   flight_path_const, min_energy,
                                   random_var_type, random_generator_type,
                                   threads)


MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::output_options,
                                                outFileBaseName,
                                                storage_interval,
                                                store_transmitted_ions,
                                                store_range_3d,
                                                store_ion_paths,
                                                store_path_limit,
                                                store_recoil_cascades,
                                                store_path_limit_recoils,
                                                store_pka
                                                )

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::target_desc,
                                                materials,
                                                regions,
                                                cell_count,
                                                periodic_bc,
                                                cell_size)

void to_json(ojson& j, const settings& p)
{
    j["Simulation"] = p.psim_;
    j["IonBeam"] = p.psrc_;
    j["Target"] = p.target_desc_;
    j["Output"] = p.out_opt_;
    for(const settings::material_desc& md : p.materials) j[md.name] = md;
    for(const settings::region_desc& rd : p.regions)  j[rd.name] = ojson(rd);
}
void from_json(const ojson& j, settings& p)
{
    if (j.contains("Simulation"))
        p.psim_ = j["Simulation"];
    if (j.contains("Output"))
        p.out_opt_ = j["Output"];
    if (j.contains(("IonBeam")))
        p.psrc_ = j["IonBeam"];
    if (!j.contains("Target")) {
        throw std::invalid_argument("Required section \"Target\" not found.");
    }
    p.target_desc_ = j["Target"];

    if (p.target_desc_.materials.empty()) {
        throw std::invalid_argument("No materials specified in \"Target\".");
    }
    if (p.target_desc_.regions.empty()) {
        throw std::invalid_argument("No regions specified in \"Target\".");
    }

    for(const std::string& m : p.target_desc_.materials)
    {
        if (j.contains(m)) {
            auto jm = j[m];
            settings::material_desc md = jm.template get< settings::material_desc >();
            p.materials.push_back(md);
        } else {
            std::string msg;
            msg = "Definition of material \"";
            msg += m;
            msg += "\" not found";
            throw std::invalid_argument(msg);
        }
    }
    for(const std::string& r : p.target_desc_.regions)
    {
        if (j.contains(r)) {
            auto jr = j[r];
            settings::region_desc rd = jr.template get< settings::region_desc >();
            p.regions.push_back(rd);
        } else {
            std::string msg;
            msg = "Definition of region \"";
            msg += r;
            msg += "\" not found";
            throw std::invalid_argument(msg);
        }
    }
}

void settings::print(std::ostream& os)
{
    ojson j(*this);
    os << j.dump(4) << endl;
}

