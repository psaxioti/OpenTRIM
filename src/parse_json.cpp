#include "mcdriver.h"

#include <fstream>
#include <iostream>

#include "json_defs_p.h"

#include "periodic_table.h"

using std::cout;
using std::cerr;
using std::endl;

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



// serialization of options struct
void to_json(ojson& j, const mcdriver::options& p);
void from_json(const ojson& j, mcdriver::options& p);

int mcdriver::options::parseJSON(std::istream& js, bool doValidation, std::ostream* os)
{
    if (os) {
        try {
            ojson j = ojson::parse(js,nullptr,true,true);
            *this = j.template get<mcdriver::options>();
            if (doValidation) validate();
        }
        catch (const ojson::exception& e) {
            *os << "Error reading json input:" << endl;
            *os << e.what() << std::endl;
            return -1;
        }
        catch (const std::invalid_argument& e) {
            *os << "Invalid option:" << endl;
            *os << e.what() << endl;
            return -1;
        }
    } else {
        ojson j = ojson::parse(js,nullptr,true,true);
        *this = j.template get<mcdriver::options>();
        if (doValidation) validate();
    }
    return 0;
}

/* serialization of element_t struct
 *
 * This is needed so that an empty element is
 * always invalid input (without needing to call
 * validate() )
 *
 * This brakes the logic that json parsing checks
 * only for json grammar
 *
 * The options.set function needs to get a full element_t definition
 * in json. Otherwise it will throw error if you try to change
 * e.g. first the symbol and then the atomic_number
 *
 * TODO
 * Must be changed in the future to a more clear
 * logic
 * e.g. throw an invalid_input exception during parsing
 * 
 */

void to_json(ojson& j, const element_t& p)
{
    j = ojson{
        {"symbol", p.symbol},
        {"atomic_number", p.atomic_number},
        {"atomic_mass", p.atomic_mass}
    };
}

void check_element_def(const ojson& j, element_t& p)
{
    if (p.symbol.empty() && p.atomic_number<=0) {
        std::stringstream msg;
        msg << "Undefined element. ";
        msg << "Specify either \"symbol\" or \"atomic_number\"";
        throw ojson::exception(ojson::other_error::create(1000, msg.str(), &j));
    }

    if (p.symbol.empty()) {
        if (p.atomic_number>dedx_max_Z) {
            std::stringstream msg;
            msg << "Element with atomic_number Z=" << p.atomic_number;
            msg << " beyong the maximum possible Z=" << dedx_max_Z;
            throw ojson::exception(ojson::other_error::create(1001, msg.str(), &j));
        }
        p.symbol = periodic_table::at(p.atomic_number).symbol;
    } else {
        auto& el = periodic_table::at(p.symbol);
        if (!el.is_valid()) {
            std::stringstream msg;
            msg << "Invalid element symbol=" << p.symbol;
            throw ojson::exception(ojson::other_error::create(1002, msg.str(), &j));
        }
        if (el.Z>dedx_max_Z) {
            std::stringstream msg;
            msg << "Element " << p.symbol << "(Z=" << el.Z << ")";
            msg << "is beyong the maximum possible Z=" << dedx_max_Z;
            throw ojson::exception(ojson::other_error::create(1002, msg.str(), &j));
        }
        if (p.atomic_number>0 && el.Z!=p.atomic_number) {
            std::stringstream msg;
            msg << "Incompatible element symbol "<< p.symbol << "(Z=" << el.Z << ")";
            msg << " and specified atomic_number Z=" << p.atomic_number;
            throw ojson::exception(ojson::other_error::create(1002, msg.str(), &j));
        }
        p.atomic_number=el.Z;
    }

    if (p.atomic_mass==0.0f) p.atomic_mass = periodic_table::at(p.atomic_number).mass;
}

void from_json(const ojson& nlohmann_json_j, element_t& nlohmann_json_t)
{
    const element_t nlohmann_json_default_obj{};
    NLOHMANN_JSON_FROM_WITH_DEFAULT(symbol);
    NLOHMANN_JSON_FROM_WITH_DEFAULT(atomic_number);
    NLOHMANN_JSON_FROM_WITH_DEFAULT(atomic_mass);
    check_element_def(nlohmann_json_j, nlohmann_json_t);
}

// enum serialization definitions

NLOHMANN_JSON_SERIALIZE_ENUM(
    ion_beam::distribution_t, {
        {ion_beam::InvalidDistribution, nullptr},
        {ion_beam::SingleValue, "SingleValue"},
        {ion_beam::Uniform, "Uniform"},
        {ion_beam::Gaussian, "Gaussian"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    ion_beam::geometry_t, {
        {ion_beam::InvalidGeometry, nullptr},
        {ion_beam::Surface, "Surface"},
        {ion_beam::Volume, "Volume"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    mccore::simulation_type_t, {
        {mccore::InvalidSimulationType, nullptr},
        {mccore::FullCascade, "FullCascade"},
        {mccore::IonsOnly, "IonsOnly"},
        {mccore::CascadesOnly, "CascadesOnly"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    mccore::scattering_calculation_t, {
        {mccore::InvalidScatteringOption, nullptr},
        {mccore::Corteo4bitTable, "Corteo4bitTable"},
        {mccore::ZBL_MAGICK, "ZBL_MAGICK"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    dedx_calc::eloss_calculation_t,
    {
        {dedx_calc::InvalidEnergyLoss, nullptr},
        {dedx_calc::EnergyLossOff, "EnergyLossOff"},
        {dedx_calc::EnergyLoss, "EnergyLoss"},
        {dedx_calc::EnergyLossAndStraggling, "EnergyLossAndStraggling"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    mccore::nrt_calculation_t, {
        {mccore::NRT_InvalidOption, nullptr},
        {mccore::NRT_element, "NRT_element"},
        {mccore::NRT_average, "NRT_average"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    flight_path_calc::flight_path_type_t, {
        {flight_path_calc::InvalidPath, nullptr},
        {flight_path_calc::AtomicSpacing, "AtomicSpacing"},
        {flight_path_calc::Constant, "Constant"},
        {flight_path_calc::MendenhallWeller, "MendenhallWeller"},
        {flight_path_calc::IPP, "IPP"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    Screening, {
        {Screening::Invalid, nullptr},
        {Screening::None, "None"},
        {Screening::LenzJensen, "LenzJensen"},
        {Screening::KrC, "KrC"},
        {Screening::Moliere, "Moliere"},
        {Screening::ZBL, "ZBL"}
    })

NLOHMANN_JSON_SERIALIZE_ENUM(
    StragglingModel, {
        {StragglingModel::Invalid, nullptr},
        {StragglingModel::Bohr, "BohrStraggling"},
        {StragglingModel::Chu,  "ChuStraggling"},
        {StragglingModel::Yang, "YangStraggling"}
    })

// option struct serialization

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    target::region,
    id,
    material_id,
    origin, size)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    atom::parameters,
    element,
    X, Ed, El, Es, Er)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    material::material_desc_t,
    id,
    density,
    composition)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    ion_beam::energy_distribution_t,
    type, center, fwhm)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    ion_beam::spatial_distribution_t,
    geometry, type, center, fwhm)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    ion_beam::angular_distribution_t,
    type, center, fwhm)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    ion_beam::parameters,
    ion,
    energy_distribution,
    spatial_distribution,
    angular_distribution)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    mccore::parameters,
    simulation_type, screening_type,
    scattering_calculation,
    eloss_calculation, straggling_model,
    nrt_calculation)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    mccore::transport_options,
    flight_path_type,
    flight_path_const, min_energy, min_recoil_energy,
    allow_sub_ml_scattering, max_mfp, max_rel_eloss)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    mcdriver::parameters,
    max_no_ions,
    threads, seed)

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    mcdriver::output_options,
    title,
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

MY_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    target::target_desc_t,
    origin,
    size,
    cell_count,
    periodic_bc,
    materials,
    regions
    )


void to_json(ojson& j, const mcdriver::options& p)
{
    j["Simulation"] = p.Simulation;
    j["Transport"] = p.Transport;
    j["IonBeam"] = p.IonBeam;
    j["Target"] = p.Target;
    j["Output"] = p.Output;
    j["Driver"] = p.Driver;
}
void from_json(const ojson& j, mcdriver::options& p)
{
    p = mcdriver::options();

    if (j.contains("Simulation"))
        p.Simulation = j["Simulation"];
    if (j.contains("Transport"))
        p.Transport = j["Transport"];
    if (j.contains("Output"))
        p.Output = j["Output"];
    if (j.contains(("IonBeam")))
        p.IonBeam = j["IonBeam"];
    if (j.contains(("Driver")))
        p.Driver = j["Driver"];
    if (!j.contains("Target")) {
        throw std::invalid_argument("Required section \"Target\" not found.");
    }
    p.Target = j["Target"];
}

void mcdriver::options::printJSON(std::ostream& os) const
{
    ojson j(*this);
    os << j.dump(4) << endl;
}

std::string mcdriver::options::toJSON() const
{
    std::ostringstream os;
    printJSON(os);
    return os.str();
}

bool mcdriver::options::set(const std::string &path, const std::string &json_str, std::ostream* os)
{
    if (os) {
    try {
        set_impl_(path,json_str);
    }
    catch (const ojson::exception& e)
    {
        *os << "Failed setting option " << path << " to " << json_str << endl
            << "Error: " << e.what() << endl
            << "exception id: " << e.id << endl;
        return false;
    }
    catch (const std::invalid_argument& e)
    {
        *os << "Failed setting option " << path << " to " << json_str << endl
            << "Error: " << e.what() << endl;
        return false;
    }
    } else set_impl_(path,json_str);
    return true;
}

bool mcdriver::options::get(const std::string &path, std::string &json_str, std::ostream* os) const
{
    if (os) {
        try {
            get_impl_(path,json_str);
        }
        catch (const ojson::exception& e)
        {
            *os << "Failed getting option " <<  path << endl
                << "error: " << e.what() << endl
                << "exception id: " << e.id << endl;
            return false;
        }
    } else get_impl_(path,json_str);
    return true;
}

void mcdriver::options::set_impl_(const std::string &path, const std::string &json_str)
{
    ojson j(*this);
    ojson::json_pointer ptr(path.c_str());
    ojson v = ojson::parse(json_str);
    j.at(ptr) = v;
    *this = j;
    validate(true);
}

void mcdriver::options::get_impl_(const std::string &path, std::string &json_str) const
{
    ojson j(*this);
    ojson::json_pointer ptr(path.c_str());
    ojson::const_reference vref = j.at(ptr);
    std::ostringstream ss;
    ss << vref.dump();
    json_str = ss.str();
}

