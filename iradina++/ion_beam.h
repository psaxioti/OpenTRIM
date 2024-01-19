#ifndef ION_BEAM_H
#define ION_BEAM_H

#include "geometry.h"

class target;
class ion;
class atom;

class ion_beam
{
public:

    typedef enum {
        SurfaceRandom = 0,
        SurfaceCentered,
        FixedPos,
        VolumeCentered,
        VolumeRandom
    } ion_distribution_t;

    struct parameters {
        ion_distribution_t ion_distribution;
        int ionZ_; // atomic number
        float ionM_; // ion mass
        float ionE0_; // initial energy eV
        vector3 dir_; // initial direction
        vector3 pos_; // initial position
        parameters();
    };

protected:

    parameters par_;
    const atom* atom_; // atomic species

public:
    ion_beam();
    ion_beam(const ion_beam& i);

    void setParameters(const parameters& p) { par_ = p; }
    const parameters& getParameters() const { return par_; }

    int ionZ() const { return par_.ionZ_; }
    float ionM() const { return par_.ionM_; }
    float ionE0() const { return par_.ionE0_; }
    ion_distribution_t distributionType() const { return par_.ion_distribution; }
    vector3 ionDir() const { return par_.dir_; }
    vector3 ionPos() const { return par_.pos_; }

    const atom* projectile() const { return atom_; }
    void setProjectile(const atom* at) { atom_ = at; }

    template<class _U>
    void source_ion(_U& g, const target& t, ion& i);

};

#endif // ION_BEAM_H
