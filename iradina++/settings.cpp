#include "settings.h"

#include <iostream>
#include <inicpp.h>

#include <chrono>
#include <thread>

using std::cout;
using std::cerr;
using std::endl;

typedef ini::IniFileCaseInsensitive inifile_t;
typedef ini::IniSectionCaseInsensitive inisection_t;

bool ini_verbose_flag_;

settings::settings() :
    cell_count({1,1,1}),
    periodic_bc({0,1,1}),
    cell_size({100.f, 100.f, 100.f})
{

}

/**
 * @brief Read an option from an ini-file section
 *
 * Read the value of an option of type T and return
 * it in the variable v.
 *
 * If the type T cannot be read from the option field
 * a std::invalid_argument error is thrown
 *
 * If the required flag is set and the option is
 * not present in the file, a std::invalid_argument error
 * is thrown
 *
 * @param s is the section of the ini file
 * @param name is the name of the option
 * @param v is the returned value for the option
 * @param required if set true, the option must be present
 */
template<typename T>
void read_option(const inisection_t& s,
                 const std::string& name,
                 T& v,
                 bool required = false)
{
    const auto& opt = s.find(name);
    if (opt!=s.end()) {
        try { v = opt->second.as<T>(); }
        catch (std::invalid_argument& s) {
            std::string msg("Error reading ");
            msg += name;
            msg += " : ";
            msg += s.what();
            throw std::invalid_argument(msg);
        }
        if (ini_verbose_flag_) cout << name << " = " << v << endl;
    } else {
        if (required) {
            std::stringstream msg;
            msg << "Error : required option " << name;
            msg << " is missing";
            throw std::invalid_argument(msg.str());
        }
    }
}

/**
 * @brief Reads an enum-typed option from the ini file
 *
 * The function reads actually an int which is then cast to
 * the enum of type T.
 *
 * TODO: maybe the option could be given as a string which will
 * be then converted to enum
 *
 * @param s is the ini section
 * @param name is the name of the option
 * @param v is the returned value
 * @param vmin the starting value of the v range
 * @param vmax the ending value of the v range
 * @param required if set shows that the option is required
 */
template<typename T>
void read_enum_option(const inisection_t& s,
                      const std::string& name,
                      T& v, int vmin, int vmax,
                      bool required = false)
{
    const auto& opt = s.find(name);
    if (opt!=s.end()) {
        try {
            int i = opt->second.as<int>();
            if (i<vmin || i>vmax) {
                std::stringstream msg;
                msg << "Error reading " << name << " : ";
                msg << "Invalid value: " << i;
                msg << "Valid input: " << vmin << "..." << vmax;
                throw std::invalid_argument(msg.str());
            }
            v = static_cast<T>(i);
        }
        catch (std::invalid_argument& s) {
            std::string msg("Error reading ");
            msg += name;
            msg += " : ";
            msg += s.what();
            throw std::invalid_argument(msg);
        }
        if (ini_verbose_flag_) cout << name << " = " << v << endl;
    } else {
        if (required) {
            std::stringstream msg;
            msg << "Error : required option " << name;
            msg << " is missing";
            throw std::invalid_argument(msg.str());
        }
    }
}

namespace ini
{
/**
 *  Conversion functor to parse a list of values from an ini field
 *
 *  The values of type T are put into a std::vector<T>
 *
 *  The list can be space- or comma-separated
 */
template<typename T>
struct Convert<std::vector<T>>
{
    /** Decodes a std::vector from a string. */
    void decode(const std::string &ivalue, std::vector<T> &result)
    {
        result.clear();

        // replace ',' by ' '
        std::string value(ivalue);
        int i;
        while ((i=value.find(','))!=std::string::npos) value.at(i) = ' ';

        std::istringstream is(value);

        while(is.good()) {
            T v;
            is >> v;
            if (is.fail()) {
                throw std::invalid_argument("Failed to parse array");
            }
            result.push_back(v);
            is >> std::ws;
            // if there is comma get it
            char c = is.peek();
            if (c==',') is.get();
        }
    }

    /** Encodes a std::vector to a string. */
    void encode(const std::vector<T> &value, std::string &result)
    {
        // variable to store the encoded element value
        std::string encoded;
        // string stream to build the result stream
        std::stringstream ss;
        for(size_t i = 0; i < value.size(); ++i)
        {
            // use the conversion functor for the type contained in
            // the vector, so the vector can use any type that
            // is compatible with inifile-cp
            Convert<T> conv;
            conv.encode(value[i], encoded);
            ss << encoded;

            // if this is not the last element add a comma as separator
            if(i != value.size() - 1)
                ss << ',';
        }
        // store the created string in the result
        result = ss.str();
    }
};
}

/**
 * @brief Read an option that contains an of values
 *
 * The values of type T are returned as a std::vector<T>.
 *
 * If required is set to true and the option is not present
 * a std::invalid_argiment is raised.
 *
 * If sz is set to a positive number, the list
 * must have exactly sz values. Otherwise a std::invalid_argument
 * is raised.
 *
 * @param s
 * @param name
 * @param v
 * @param required
 * @param sz
 */
template<typename T>
void read_option_array(const inisection_t& s,
                       const std::string& name,
                       std::vector<T>& v, bool required = false,
                       int sz = -1)
{
    const auto& opt = s.find(name);
    if (opt!=s.end()) {
        try { v = opt->second.as< std::vector<T> >(); }
        catch (std::invalid_argument& s) {
            std::string msg("Error reading ");
            msg += name;
            msg += " : ";
            msg += s.what();
            throw std::invalid_argument(msg);
        }
        if (sz > 0 && v.size()!=sz) {
            std::stringstream msg;
            msg << "Error reading " << name << " : ";
            msg << "Must be a " << sz << "-element vector";
            throw std::invalid_argument(msg.str());
        }
        if (ini_verbose_flag_) {
            cout << name << "[" << v.size() << "] = ";
            for(int i=0; i< std::min(3UL,v.size()); i++)
                cout << v[i] << ' ';
            if (v.size() > 3) cout << "...";
            cout << endl;
        }
    } else {
        if (required) {
            std::stringstream msg;
            msg << "Error : required option " << name;
            msg << " is missing";
            throw std::invalid_argument(msg.str());
        }
    }
}

template<typename T>
void read_option_v3(const inisection_t& s,
                    const std::string& name,
                    Eigen::AlignedVector3<T>& v3,
                    bool required = false)
{
    const auto& opt = s.find(name);
    if (opt!=s.end()) {
        std::vector<T> temp;
        read_option_array(s,name,temp,required,3);
        v3[0] = temp[0];
        v3[1] = temp[1];
        v3[2] = temp[2];
    }
}

int parse_material(const inisection_t& s, settings::material_desc& md) {
    int count = 0;
    md.density = 0;
    try {
        read_option(s,"ElementCount",count,true);
        if (count < 1) {
            cerr << "Invalid ElementCount in material " << md.name << endl;
            return -1;
        }
        read_option(s,"Density",md.density,true);
        if (md.density <= 0.f) {
            cerr << "Invalid Density" << endl;
            return -1;
        }
        read_option_array(s,"ElementsZ",            md.Z, true, count);
        read_option_array(s,"ElementsM",            md.M, true, count);
        read_option_array(s,"ElementsConc",         md.X, true, count);
        read_option_array(s,"ElementsDispEnergy",   md.Ed,true, count);
        read_option_array(s,"ElementsLattEnergy",   md.El,true, count);
        read_option_array(s,"ElementsSurfEnergy",   md.Es,true, count);
        read_option_array(s,"ElementsReplEnergy",   md.Er,true, count);

    } catch (std::invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }

    return 0;
}

int parse_region(const inisection_t& s,
                 const std::vector<std::string>& material_ids,
                 settings::region_desc& rd)
{
    try {
        std::string mname;
        read_option(s,"material",mname,true);
        auto i = std::find(material_ids.begin(), material_ids.end(),mname);
        if (i == material_ids.end()) {
            cerr << "In region " << rd.name << " the specified material " << mname
                 << " is not in the [Target] materials list " << endl;
            return -1;
        }
        rd.material_id = i - material_ids.begin();

        read_option_array(s,"extent_x", rd.extX, true, 2);
        read_option_array(s,"extent_y", rd.extY, true, 2);
        read_option_array(s,"extent_z", rd.extZ, true, 2);

    } catch (std::invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}


int settings::parse(const char* fname, bool verbose)
{
    ini_verbose_flag_ = verbose;

    if (verbose) cout << "Parsing config file " << fname << endl;

    // load the config file
    ini::IniFileCaseInsensitive fconfig;

    try {
        fconfig.load(fname);
    } catch (std::exception& e) {
        cerr << "Error reading config: " << e.what() << endl;
        return -1;
    }

    // Simulation
     {
        const auto& opts = fconfig.find("Simulation");
        if (opts != fconfig.end()) {
            if (verbose) cout << "[Simulation]" << endl;
            try {
                read_option(opts->second,"title", psim_.title, true);
                read_option(opts->second,"max_no_ions", psim_.max_no_ions, true);
                read_enum_option(opts->second,"simulation_type", psim_.simulation_type,
                                 simulation_base::FullCascade, simulation_base::KP2,
                                 true);

                // iradina does not define simulation type 1 & 2
                if (psim_.simulation_type == simulation_base::Invalid1 ||
                    psim_.simulation_type == simulation_base::Invalid2)
                {
                    std::stringstream msg;
                    msg << "simulation_type: Invalid option " << psim_.simulation_type;
                    throw std::invalid_argument(msg.str());
                }

                // KP2 not yet implementd
                if (psim_.simulation_type == simulation_base::KP2)
                {
                    std::stringstream msg;
                    msg << "simulation_type=" << psim_.simulation_type;
                    msg << " not yet implemented.";
                    throw std::invalid_argument(msg.str());
                }

                read_enum_option(opts->second,"scattering_calculation", psim_.scattering_calculation,
                                 simulation_base::Corteo4bit, simulation_base::ZBL_MAGICK);
                read_enum_option(opts->second,"flight_path_type", psim_.flight_path_type,
                                 simulation_base::Poisson, simulation_base::SRIMlike);
                read_enum_option(opts->second,"straggling_model", psim_.straggling_model,
                                 simulation_base::NoStraggling, simulation_base::YangStraggling);

                read_option(opts->second,"flight_path_const", psim_.flight_path_const);
                read_option(opts->second,"min_energy", psim_.min_energy);

                read_enum_option(opts->second,"random_var_type", psim_.random_var_type,
                                 simulation_base::Sampled, simulation_base::Tabulated);
                read_enum_option(opts->second,"random_generator_type", psim_.random_generator_type,
                                 simulation_base::MinStd, simulation_base::MersenneTwister);
                read_option(opts->second,"threads", psim_.threads);
                read_option_array(opts->second,"seeds", psim_.seeds);


            } catch (std::invalid_argument& e) {
                cerr << e.what() << endl;
                return -1;
            }
        } else {
            cerr << "Required section [Simulation] not found" << endl;
            return -1;
        }
    }


    // ion beam
    {
        const auto& opts = fconfig.find("IonBeam");
        if (opts != fconfig.end()) {
            cout << "[IonBeam]" << endl;
            try {
                read_enum_option(opts->second,"ion_distribution", psrc_.ion_distribution,
                                 ion_beam::SurfaceRandom, ion_beam::VolumeRandom);
                read_option(opts->second,"ionZ",psrc_.ionZ,true);
                read_option(opts->second,"ionM",psrc_.ionM,true);
                read_option(opts->second,"ionE0",psrc_.ionE0,true);

                read_option_v3(opts->second,"ion_dir",psrc_.dir);
                psrc_.dir.normalize();
                read_option_v3(opts->second,"ion_pos",psrc_.pos);

            } catch (std::invalid_argument& e) {
                cerr << e.what() << endl;
                return -1;
            }
        } else {
            cerr << "Error : Required section [IonBeam] not found" << endl;
            return -1;
        }
    }

    // target
    std::vector<std::string> material_ids, region_ids;
    {
        const auto& opts = fconfig.find("Target");
        if (opts != fconfig.end()) {
            cout << "[Target]" << endl;
            try {
                read_option_v3(opts->second,"cell_count",cell_count,true);
                read_option_v3(opts->second,"cell_size",cell_size,true);
                read_option_v3(opts->second,"periodic_boundary",periodic_bc,true);
                read_option_array(opts->second,"materials",material_ids,true);
                read_option_array(opts->second,"regions",region_ids,true);
            } catch (std::invalid_argument& e) {
                cerr << e.what() << endl;
                return -1;
            }
        } else {
            cerr << "Required section [Target] not found" << endl;
            return -1;
        }
        if (material_ids.empty()) {
            cerr << "Error : No materials specified." << endl;
            return -1;
        }
        if (region_ids.empty()) {
            cerr << "Error : No regions specified." << endl;
            return -1;
        }

    }

    // check that we have readable materials and regions
    for(std::string name : material_ids) {
        if (verbose) cout << "[" << name << "]" << endl;
        const auto& opts = fconfig.find(name);
        if (opts == fconfig.end()) {
            cerr << "Error : material section [" << name << "] not found." << endl;
            return -1;
        }
        material_desc md;
        md.name = name;
        if (parse_material(opts->second, md)==0) materials.push_back(md);
        else {
            cerr << "Error reading material " << name << endl;
            return -1;
        }

    }
    for(std::string name : region_ids) {
        if (verbose) cout << "[" << name << "]" << endl;
        const auto& opts = fconfig.find(name);
        if (opts == fconfig.end()) {
            cerr << "Error : region section [" << name << "] not found." << endl;
            return -1;
        }
        region_desc rd;
        rd.name = name;
        if (parse_region(opts->second, material_ids, rd)==0) regions.push_back(rd);
        else {
            cerr << "Error reading region " << name << endl;
            return -1;
        }
    }
    return 0;
}

static const char* ini_template_ =
    "# Template configuration file for iradina++\n"
    "# \n"
    "# The input is divided in sections.\n"
    "# All sections are required\n"
    "#\n"
    "\n"
    "[Simulation]                        # Required Section\n"
    "title = A short description         # Required. Free string\n"
    "max_no_ions=10000                   # Required.\n"
    "simulation_type=0                   # Required. Values 0,1,..,5\n"
    "scattering_calculation=0            # Optional. Values 0(default),1,2\n"
    "flight_length_type=0                # Optional. Values 0(def),...,3\n"
    "flight_length_constant=0.0          # Optional. [nm] Set only if flight_length_type=2\n"
    "straggling_model=3                  # Optional. Values 0,1,2,3(default)\n"
    "min_energy=5.0                      # Optional. Defaults to 5. [eV]\n"
    "random_var_type = 1                 # Optional. Values 0(default), 1\n"
    "random_generator = 0                # Optional. Values 0(default), 1\n"
    "\n"
    "[IonBeam]                           # Required section\n"
    "ionZ=26                             # Required\n"
    "ionM=53.476                         # Required\n"
    "ionE0=2000000.0                     # Required [eV]\n"
    "ion_dir=1, 1, 0                     # Optional. Default = 1 0 0\n"
    "ion_distribution=1                  # Optional. Values 0(default),1,2,3\n"
    "ion_pos=0 0 0                       # Optional. Def=0 0 0. [nm]. Only if required by distribution\n"
    "\n"
    "[Target]                            # Required section\n"
    "Materials = Fe                      # Required. A list of string identifiers for the materials\n"
    "regions = R1                        # Required. A list of string identifiers for the regions\n"
    "cell_count=100 1 1                  # Required. 3-element int array (x,y,z)\n"
    "cell_size=12.0 1200 1200            # Required. [nm]. 3-element float array (x,y,z)\n"
    "periodic_boundary=0 1 1             # Required. 3-element bool array (0=false or 1=true)\n"
    "\n"
    "[Fe]                                # Required. One section for each material\n"
    "                                    # The name corresponds to the one given in [Target]/Materials\n"
    "ElementCount=1                      # Required\n"
    "Density=7.874                       # Required. [g/cm^3]\n"
    "ElementsZ=26                        # Required. Array of ElementCount values\n"
    "ElementsM=53.476                    # Required. Array of ElementCount values\n"
    "ElementsConc=1                      # Required. Array of ElementCount values\n"
    "ElementsDispEnergy=40               # Required. Array of ElementCount values\n"
    "ElementsLattEnergy=3                # Required. Array of ElementCount values\n"
    "ElementsSurfEnergy=4.3              # Required. Array of ElementCount values\n"
    "ElementsReplEnergy=40               # Required. Array of ElementCount values\n"
    "\n"
    "\n"
    "[R1]                                # Required. One section for each region\n"
    "                                    # The name corresponds to the one given in [Target]/Regions\n"
    "material=Fe                         # Required. Identifier of material to fill the region\n"
    "extent_x=0 1200                     # Required. Region extent Xmin, Xmax\n"
    "extent_y=0 1200                     # Required. Region extent Ymin, Ymax\n"
    "extent_z=0 1200                     # Required. Region extent Zmin, Zmax\n"
    "\n";

const char* settings::ini_template() const
{
    return ini_template_;
}

simulation_base* settings::createSimulation() const
{
    simulation_base* S = simulation_base::fromParameters(psim_);
    S->setIonBeam(psrc_);

    for(const material_desc& md : materials) {
        material* m = S->addMaterial(md.name.c_str());
        m->setMassDensity(md.density);
        for(int i=0; i<md.Z.size(); i++)
            m->addAtom(
                {.Z=md.Z[i],.M=md.M[i],.Ed=md.Ed[i],
                 .El=md.El[i],.Es=md.Es[i],.Er=md.Er[i]},
                md.X[i]);
    }

    grid3D& G = S->grid();
    G.setX(0, cell_count.x()*cell_size.x(), cell_count.x());
    G.setY(0, cell_count.y()*cell_size.y(), cell_count.y());
    G.setZ(0, cell_count.z()*cell_size.z(), cell_count.z());

    const std::vector<material*>& imat = S->getTarget()->materials();
    for(const region_desc& rd : regions) {
        box3D box;
        box.min() = vector3(rd.extX[0],rd.extY[0],rd.extZ[0]);
        box.max() = vector3(rd.extX[1],rd.extY[1],rd.extZ[1]);
        S->fill(box,imat[rd.material_id]);
    }

    return S;
}



