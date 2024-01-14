#include "simulation.h"

#include <iostream>
#include <inicpp.h>

using std::cout;
using std::cerr;
using std::endl;

typedef ini::IniFileCaseInsensitive inifile_t;
typedef ini::IniSectionCaseInsensitive inisection_t;

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
        cout << name << " = " << v << endl;
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
        cout << name << " = " << v << endl;
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
        cout << name << "[" << v.size() << "] = ";
        for(int i=0; i< std::min(3UL,v.size()); i++)
            cout << v[i] << ' ';
        if (v.size() > 3) cout << "...";
        cout << endl;
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

struct material_data {
    std::string name;
    int count;
    float density;
    std::vector<int> Z;
    std::vector<float> M;
    std::vector<float> X;
    std::vector<float> Ed;
    std::vector<float> El;
    std::vector<float> Es;
    std::vector<float> Er;
};

int parse_material(const inisection_t& s, material_data& md) {
    md.count = 0;
    md.density = 0;
    try {
        read_option(s,"ElementCount",md.count,true);
        if (md.count < 1) {
            cerr << "Invalid ElementCount in " << endl;
            return -1;
        }
        read_option(s,"Density",md.density,true);
        if (md.density <= 0.f) {
            cerr << "Invalid Density" << endl;
            return -1;
        }
        read_option_array(s,"ElementsZ",            md.Z, true, md.count);
        read_option_array(s,"ElementsM",            md.M, true, md.count);
        read_option_array(s,"ElementsConc",         md.X, true, md.count);
        read_option_array(s,"ElementsDispEnergy",   md.Ed,true, md.count);
        read_option_array(s,"ElementsLattEnergy",   md.El,true, md.count);
        read_option_array(s,"ElementsSurfEnergy",   md.Es,true, md.count);
        read_option_array(s,"ElementsReplEnergy",   md.Er,true, md.count);

    } catch (std::invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}

struct region_data {
    int material_id;
    std::vector<float> extX;
    std::vector<float> extY;
    std::vector<float> extZ;
};

int parse_region(const inisection_t& s,
                 const std::vector<std::string>& materials,
                 region_data& rd)
{
    try {
        std::string mname;
        read_option(s,"material",mname,true);
        auto i = std::find(materials.begin(), materials.end(),mname);
        if (i == materials.end()) {
            cerr << "Material " << mname << " is not in the [Targey] materials list " << endl;
            return -1;
        }
        rd.material_id = i - materials.begin();

        read_option_array(s,"extent_x", rd.extX, true, 2);
        read_option_array(s,"extent_y", rd.extY, true, 2);
        read_option_array(s,"extent_z", rd.extZ, true, 2);

    } catch (std::invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}


int parse(std::string& fname, simulation_base* &S)
{
    cout << "Parsing config file " << fname << endl;

    // load the config file
    ini::IniFileCaseInsensitive fconfig;

    try {
        fconfig.load(fname);
    } catch (std::exception& e) {
        cerr << "Error reading config: " << e.what() << endl;
        return -1;
    }

    // Simulation
    simulation_base::parameters p;
    {
        const auto& opts = fconfig.find("Simulation");
        if (opts != fconfig.end()) {
            cout << "[Simulation]" << endl;
            try {
                read_option(opts->second,"title", p.title, true);
                read_option(opts->second,"max_no_ions", p.max_no_ions, true);
                read_enum_option(opts->second,"simulation_type", p.simulation_type,
                                 simulation_base::FullCascade, simulation_base::KP2,
                                 true);

                // iradina does not define simulation type 1 & 2
                if (p.simulation_type == simulation_base::Invalid1 ||
                    p.simulation_type == simulation_base::Invalid2)
                {
                    std::stringstream msg;
                    msg << "simulation_type: Invalid option " << p.simulation_type;
                    throw std::invalid_argument(msg.str());
                }

                // KP2 not yet implementd
                if (p.simulation_type == simulation_base::KP2)
                {
                    std::stringstream msg;
                    msg << "simulation_type=" << p.simulation_type;
                    msg << " not yet implemented.";
                    throw std::invalid_argument(msg.str());
                }

                read_enum_option(opts->second,"scattering_calculation", p.scattering_calculation,
                                 simulation_base::Corteo4bit, simulation_base::ZBL_MAGICK);
                read_enum_option(opts->second,"flight_path_type", p.flight_path_type,
                                 simulation_base::Poisson, simulation_base::SRIMlike);
                read_enum_option(opts->second,"straggling_model", p.straggling_model,
                                 simulation_base::NoStraggling, simulation_base::YangStraggling);

                read_option(opts->second,"flight_path_const", p.flight_path_const);
                read_option(opts->second,"min_energy", p.min_energy);

                read_enum_option(opts->second,"random_var_type", p.random_var_type,
                                 simulation_base::Sampled, simulation_base::Tabulated);
                read_enum_option(opts->second,"random_generator_type", p.random_generator_type,
                                 simulation_base::MinStd, simulation_base::MersenneTwister);


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
    ion_beam::parameters ibp;
    {
        const auto& opts = fconfig.find("IonBeam");
        if (opts != fconfig.end()) {
            cout << "[IonBeam]" << endl;
            try {
                read_enum_option(opts->second,"ion_distribution", ibp.ion_distribution,
                                 ion_beam::SurfaceRandom, ion_beam::VolumeRandom);
                read_option(opts->second,"ionZ",ibp.ionZ_,true);
                read_option(opts->second,"ionM",ibp.ionM_,true);
                read_option(opts->second,"ionE0",ibp.ionE0_,true);

                read_option_v3(opts->second,"ion_dir",ibp.dir_);
                ibp.dir_.normalize();
                read_option_v3(opts->second,"ion_pos",ibp.pos_);

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
    ivector3 cell_count = {1, 1, 1};
    ivector3 periodic_bc = {0, 1, 1};
    vector3 cell_size = {100.f, 100.f, 100.f};
    std::vector<std::string> materials, regions;
    {
        const auto& opts = fconfig.find("Target");
        if (opts != fconfig.end()) {
            cout << "[Target]" << endl;
            try {
                read_option_v3(opts->second,"cell_count",cell_count,true);
                read_option_v3(opts->second,"cell_size",cell_size,true);
                read_option_v3(opts->second,"periodic_boundary",periodic_bc,true);
                read_option_array(opts->second,"materials",materials,true);
                read_option_array(opts->second,"regions",regions,true);
            } catch (std::invalid_argument& e) {
                cerr << e.what() << endl;
                return -1;
            }
        } else {
            cerr << "Required section [IonBeam] not found" << endl;
            return -1;
        }
        if (materials.empty()) {
            cerr << "Error : No materials specified." << endl;
            return -1;
        }
        if (regions.empty()) {
            cerr << "Error : No regions specified." << endl;
            return -1;
        }

    }

    // check that we have readable materials and regions
    std::vector<material_data> mdat;
    for(std::string name : materials) {
        cout << "[" << name << "]" << endl;
        const auto& opts = fconfig.find(name);
        if (opts == fconfig.end()) {
            cerr << "Error : material section [" << name << "] not found." << endl;
            return -1;
        }
        material_data md;
        md.name = name;
        if (parse_material(opts->second, md)==0) mdat.push_back(md);
        else {
            cerr << "Error reading material " << name << endl;
            return -1;
        }

    }
    std::vector<region_data> rdat;
    for(std::string name : regions) {
        cout << "[" << name << "]" << endl;
        const auto& opts = fconfig.find(name);
        if (opts == fconfig.end()) {
            cerr << "Error : region section [" << name << "] not found." << endl;
            return -1;
        }
        region_data rd;
        if (parse_region(opts->second, materials, rd)==0) rdat.push_back(rd);
        else {
            cerr << "Error reading region " << name << endl;
            return -1;
        }
    }

    S = simulation_base::fromParameters(p);
    S->setIonBeam(ibp);

    for(const material_data& md : mdat) {
        material* m = S->addMaterial(md.name.c_str(),md.density);
        for(int i=0; i<md.count; i++)
            m->addAtom(md.Z[i],md.M[i],md.X[i],md.Ed[i],md.El[i],md.Es[i],md.Er[i]);
    }

    grid3D& G = S->grid();
    G.setX(0, cell_count.x()*cell_size.x(), cell_count.x());
    G.setY(0, cell_count.y()*cell_size.y(), cell_count.y());
    G.setZ(0, cell_count.z()*cell_size.z(), cell_count.z());

    const std::vector<material*>& imat = S->getInventory().materials();
    for(const region_data& rd : rdat) {
        box3D box;
        box.min() = vector3(rd.extX[0],rd.extY[0],rd.extZ[0]);
        box.max() = vector3(rd.extX[1],rd.extY[1],rd.extZ[1]);
        S->fill(box,imat[rd.material_id]);
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

void print_ini_template()
{
    cout << ini_template_;
}
