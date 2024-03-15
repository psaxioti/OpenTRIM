#include "ion_beam.h"
#include "random_vars.h"
#include "ion.h"
#include "target.h"

//ion_beam::parameters::parameters() {
//    ion_distribution = ion_beam::SurfaceRandom;
//    ionZ_ = 1; // proton beam
//    ionM_ = 1; // ion mass
//    ionE0_ = 1e6f; // initial energy eV
//    dir_ = {1.f, 0.f, 0.f}; // initial direction
//    pos_ = {0.f, 0.f, 0.f}; // initial position
//}

ion_beam::ion_beam() :
    par_(), atom_(nullptr)
{

}

ion_beam::ion_beam(const ion_beam& i) :
    par_(i.par_), atom_(i.atom_)
{

}

template<class _U>
void ion_beam::source_ion(_U& g, const target& t, ion& i)
{
    i.setGrid(&(t.grid()));

    // beam in +x direction
    i.dir() = par_.dir;

    const grid1D& x = t.grid().x();
    const grid1D& y = t.grid().y();
    const grid1D& z = t.grid().z();

    vector3 p(0.,0.,0.);

    // initial position
    switch (par_.ion_distribution) {
    case SurfaceRandom:
        p.x() = x.front();
        p.y() = y.front() + y.w()*g.u01open();
        p.z() = z.front() + z.w()*g.u01open();
        break;
    case SurfaceCentered:
        p.x() = x.front();
        p.y() = y.front() + y.w()/2;
        p.z() = z.front() + z.w()/2;
        break;
    case FixedPos:
        p = par_.pos;
    case VolumeRandom:
        p.x() = x.front() + x.w()*g.u01open();
        p.y() = y.front() + y.w()*g.u01open();
        p.z() = z.front() + z.w()*g.u01open();
        break;
    case VolumeCentered:
        p.x() = x.front() + x.w()/2;
        p.y() = y.front() + y.w()/2;
        p.z() = z.front() + z.w()/2;
        break;
    default:
        break;
    }
    i.setPos(p);
    i.erg() = par_.ionE0;
    i.myAtom() = atom_;

}

template void ion_beam::source_ion< rng_mt >(rng_mt& g, const target& t, ion& i);
template void ion_beam::source_ion< rng_msrand >(rng_msrand& g, const target& t, ion& i);
template void ion_beam::source_ion< rng_xoshiro >(rng_xoshiro& g, const target& t, ion& i);
