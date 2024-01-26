#include "settings.h"

#include <fstream>

namespace nlohmann {

template <typename T>
struct adl_serializer< Eigen::AlignedVector3<T> > {
    static void to_json(json& j, const Eigen::AlignedVector3<T>& V3) {
        j = {V3.x(), V3.y(), V3.z()};
    }

    static void from_json(const json& j, Eigen::AlignedVector3<T>& V3) {
        std::array<T,3> v = j.template get< std::array<T,3> >();
        V3 = Eigen::AlignedVector3<T>( v[0], v[1], v[2] );
    }
};

} // end namespace

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::region_desc,
                                   name,
                                   material_id,
                                   extX, extY, extZ)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::material_desc,
                                   name,
                                   density, isMassDensity,
                                   Z, M, X, Ed, El, Es, Er);

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

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(ion_beam::parameters,
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

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::parameters,
                                   title, max_no_ions, simulation_type,
                                   nrt_calculation, scattering_calculation,
                                   flight_path_type, straggling_model,
                                   flight_path_const, min_energy,
                                   random_var_type, random_generator_type,
                                   threads)


NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(simulation_base::output_options,
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

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(settings::target_desc,
                                                materials,
                                                regions,
                                                cell_count,
                                                periodic_bc,
                                                cell_size)

void to_json(json& j, const settings& p)
{

}
void from_json(const json& j, settings& p)
{
    p.psim_ = j["Simulation"];
    p.out_opt_ = j["Output"];
    p.psrc_ = j["IonBeam"];
    p.target_desc_ = j["Target"];

    for(const std::string& m : p.target_desc_.materials)
    {
        auto jm = j[m];
        settings::material_desc md = jm.template get< settings::material_desc >();
        p.materials.push_back(md);
    }
    for(const std::string& r : p.target_desc_.regions)
    {
        auto jr = j[r];
        settings::region_desc rd = jr.template get< settings::region_desc >();
        p.regions.push_back(rd);
    }
}

int settings::parse(json j)
{

    return 0;
}

int settings::test() {

    std::ifstream f("FeFC.json");
    json j = json::parse(f);


    settings p =
        j.template get<settings>();

    return 0;
}
