#include "options.h"
#include "elements.h"

#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cctype>

#define CHECK_EMPTY_VEC(MatDesc,V,MatName) \
    if (MatDesc.V.size()==0) throw std::invalid_argument("Empty "#V" in " + MatName);

#define CHECK_VEC_SIZE(MatDesc,V,N,MatName) \
    if (MatDesc.V.size()!=N) throw std::invalid_argument("In material descriptor " + MatName + " the # of values in "#V" is not equal to those of Z." );

#define ANY_ZEROorNEG(V) \
    std::any_of(V.begin(),V.end(),[](float c){ return c<=0.f; })

#define CHECK_VEC_ZEROorNEG(MatDesc,V,MatName) \
    if (ANY_ZEROorNEG(MatDesc.V)) throw std::invalid_argument("In material descriptor " + MatName + " "#V" has zero or negative values." );




int options::validate()
{
    // Simulation
    if (Simulation.max_no_ions <= 0)
        throw std::invalid_argument("Simulation.max_no_ions must be larger than 0.");

    if (Simulation.flight_path_type==simulation::Constant &&
        Simulation.flight_path_const<=0.f)
        throw std::invalid_argument("Simulation.flight_path_type is \"Constant\" but Simulation.flight_path_const is negative.");

    if (Simulation.threads <=0)
        throw std::invalid_argument("Simulation.threads must be 1 or larger.");

    if (!Simulation.seeds.empty()) {
        if (Simulation.seeds.size()!=Simulation.threads) {
            std::stringstream ss;
            ss << "Simulation.threads=" << Simulation.threads << " while the # of "
               << "Simulation.seeds is " << Simulation.seeds.size() << "." << std::endl
               << "Either enter a # of seeds equal to the # of threads or no seeds at all.";
            throw std::invalid_argument(ss.str());
        }
    }

    // Output
    const std::string& fname = Output.OutputFileBaseName;
    if (fname.empty())
        throw std::invalid_argument("Output.OutputFileBaseName is empty.");

    if (std::any_of(fname.begin(),
                    fname.end(),
                    [](unsigned char c){ return !std::isalnum(c); }))
    {
        std::string msg = "Output.OutputFileBaseName=\"";
        msg += fname;
        msg += "\" contains non alphanumeric characters.";
        throw std::invalid_argument(msg);
    }

    // Target
    if (Target.materials.empty())
        throw std::invalid_argument("Target.materials is empty.");
    if (Target.regions.empty())
        throw std::invalid_argument("Target.regions is empty.");
    const ivector3& cell_count = Target.cell_count;
    if (std::any_of(cell_count.begin(),
                    cell_count.end(),
                    [](int c){ return c<=0; }))
    {
        std::stringstream msg;
        msg << "Target.cell_count=[";
        msg << cell_count;
        msg << "] contains zero or negative values.";
        throw std::invalid_argument(msg.str());
    }
    const vector3& cell_size = Target.cell_size;
    if (std::any_of(cell_size.begin(),
                    cell_size.end(),
                    [](float c){ return c<=0.f; }))
    {
        std::stringstream msg;
        msg << "Target.cell_size=[";
        msg << cell_size;
        msg << "] contains zero or negative values.";
        throw std::invalid_argument(msg.str());
    }

    if (Target.materials.size()!=materials_desc.size())
        throw std::invalid_argument("The # of Target.materials is not equal to "
                                    "the # of material descriptors.");
    if (Target.regions.size()!=regions_desc.size())
        throw std::invalid_argument("The # of Target.regions is not equal to "
                                    "the # of region descriptors.");

    // Material descriptors
    int k = 0;
    for(auto m : materials_desc) {
        auto mname = Target.materials[k++];
        if (m.density <= 0.f) {
            std::stringstream msg;
            msg << "Zero or negative density in material";
            msg << mname;
            throw std::invalid_argument(msg.str());
        }
        CHECK_EMPTY_VEC(m,Z,mname)
        CHECK_EMPTY_VEC(m,M,mname)
        CHECK_EMPTY_VEC(m,X,mname)
        CHECK_EMPTY_VEC(m,Ed,mname)
        CHECK_EMPTY_VEC(m,El,mname)
        CHECK_EMPTY_VEC(m,Es,mname)
        CHECK_EMPTY_VEC(m,Er,mname)

        int natoms = m.Z.size();
        if (std::any_of(m.Z.begin(),
                    m.Z.end(),
                    [](int z){ return (z < 1) || (z > elements::max_atomic_num);}))

            throw std::invalid_argument("Invalid Z number in material" + mname);

        CHECK_VEC_SIZE(m,M,natoms,mname)
        CHECK_VEC_SIZE(m,X,natoms,mname)
        CHECK_VEC_SIZE(m,Ed,natoms,mname)
        CHECK_VEC_SIZE(m,El,natoms,mname)
        CHECK_VEC_SIZE(m,Es,natoms,mname)
        CHECK_VEC_SIZE(m,Er,natoms,mname)
        CHECK_VEC_ZEROorNEG(m,M,mname)
        CHECK_VEC_ZEROorNEG(m,X,mname)
        CHECK_VEC_ZEROorNEG(m,Ed,mname)
        CHECK_VEC_ZEROorNEG(m,El,mname)
        CHECK_VEC_ZEROorNEG(m,Es,mname)
        CHECK_VEC_ZEROorNEG(m,Er,mname)
    }


    return 0;
}

simulation* options::createSimulation() const
{
    simulation* S = new simulation(Simulation);
    S->setOutputOptions(Output);
    S->setIonBeam(IonBeam);

    target& T = S->getTarget();

    for(int i=0; i<materials_desc.size(); i++) {
        material* m = T.addMaterial(Target.materials[i].c_str());
        const material::material_desc_t& md = materials_desc[i];
        m->setMassDensity(md.density);
        for(int i=0; i<md.Z.size(); i++)
            m->addAtom(
                {.Z=md.Z[i],.M=md.M[i],.Ed=md.Ed[i],
                 .El=md.El[i],.Es=md.Es[i],.Er=md.Er[i]},
                md.X[i]);
    }

    grid3D& G = T.grid();
    G.setX(0, Target.cell_count.x()*Target.cell_size.x(),
           Target.cell_count.x(), Target.periodic_bc.x());
    G.setY(0, Target.cell_count.y()*Target.cell_size.y(),
           Target.cell_count.y(), Target.periodic_bc.y());
    G.setZ(0, Target.cell_count.z()*Target.cell_size.z(),
           Target.cell_count.z(), Target.periodic_bc.z());

    const std::vector<material*>& imat = T.materials();
    for(const target::region& rd : regions_desc) {
        box3D box;
        box.min() = rd.min;
        box.max() = rd.max;
        int i = materialIdx(rd.material_id);
        T.fill(box,imat[i]);
    }

    return S;
}


