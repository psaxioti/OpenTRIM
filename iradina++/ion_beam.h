#ifndef ION_BEAM_H
#define ION_BEAM_H

class target;
class ion;
class atom;

class ion_beam
{
public:

    typedef enum {
        SurfaceRandom = 0,
        SurfaceCentered,
        VolumeRandom,
        VolumeCentered
    } source_type_t;

protected:

    source_type_t source_type;
    const atom* atom_; // atomic species
    float E0_; // initial energy eV
    unsigned int counter;

public:

    ion_beam();

    const atom* projectile() const { return atom_; }
    void setProjectile(const atom* at, float E0)
    { atom_ = at; E0_ = E0; }

    template<class _U>
    void source_ion(_U& g, const target& t, ion& i);

};

#endif // ION_BEAM_H
