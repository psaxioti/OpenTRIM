#include "ion_beam.h"
#include "urbg.h"
#include "ion.h"
#include "target.h"

ion_beam::ion_beam() : source_type(SurfaceRandom),
    atom_(nullptr), E0_(1.e6f), counter(0)
{

}

void ion_beam::source_ion(URBG& g, const target& t, ion& i)
{
    // beam in +x direction
    i.dir = {1.f, 0.f, 0.f};

    const grid1D& x = t.grid().x();
    const grid1D& y = t.grid().y();
    const grid1D& z = t.grid().z();

    vector3 p(0.,0.,0.);

    // initial position
    switch (source_type) {
    case SurfaceRandom:
        p.x() = x.front();
        p.y() = y.front() + y.w*g.u01open();
        p.z() = z.front() + z.w*g.u01open();
        break;
    case SurfaceCentered:
        p.x() = x.front();
        p.y() = y.front() + y.w/2;
        p.z() = z.front() + z.w/2;
        break;
    case VolumeRandom:
        p.x() = x.front() + x.w*g.u01open();
        p.y() = y.front() + y.w*g.u01open();
        p.z() = z.front() + z.w*g.u01open();
        break;
    case VolumeCentered:
        p.x() = x.front() + x.w/2;
        p.y() = y.front() + y.w/2;
        p.z() = z.front() + z.w/2;
        break;
    }
    i.setPos(p);
    i.erg = E0_;
    i.atom_ = atom_;
    i.ion_id = ++counter;
    i.recoil_id = 0;
}
