#include "options.h"

#include <stdexcept>
#include <sstream>

options::options()
{}

int options::validate()
{
    // Simulation
    if (sim_par.validate()!=0 ||
        out_opt.validate()!=0) return -1;
    return 0;
}

simulation_base* options::createSimulation() const
{
    simulation_base* S = simulation_base::fromParameters(sim_par);
    S->setOutputOptions(out_opt);
    S->setIonBeam(src_par);

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
    G.setX(0, target_desc.cell_count.x()*target_desc.cell_size.x(), target_desc.cell_count.x());
    G.setY(0, target_desc.cell_count.y()*target_desc.cell_size.y(), target_desc.cell_count.y());
    G.setZ(0, target_desc.cell_count.z()*target_desc.cell_size.z(), target_desc.cell_count.z());

    const std::vector<material*>& imat = S->getTarget()->materials();
    for(const region_desc& rd : regions) {
        box3D box;
        box.min() = vector3(rd.extX[0],rd.extY[0],rd.extZ[0]);
        box.max() = vector3(rd.extX[1],rd.extY[1],rd.extZ[1]);
        int i = mat2idx.at(rd.material_id);
        S->fill(box,imat[i]);
    }

    return S;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << v[0];
    for(int i=1; i<v.size(); i++)
        os << ", " << v[i];
    return os;
}





