#ifndef ION_BEAM_H
#define ION_BEAM_H

#include "geometry.h"

class target;
class ion;
class atom;
class random_vars;

/**
 * @brief The ion_beam class is used for ejecting ions into the simulation
 *
 * The ions are generated according to the spatial, direction and energy distributions
 * specified by the class \ref parameters.
 *
 * @ingroup Ions
 */
class ion_beam
{
public:

    /**
     * @brief Types of source ions spatial distribution
     */
    enum ion_distribution_t {
        SurfaceRandom = 0,          /**< Randomly distributed on the entrance surface */
        SurfaceCentered,            /**< At the center of the entrance surface */
        FixedPos,                   /**< At a fixed posistion */
        VolumeCentered,             /**< At the center of the target volume */
        VolumeRandom,               /**< Randomly distributed in the target volume */
        InvalidIonDistribution = -1
    };

    /**
     * @brief The parameters of the ion_beam class
     */
    struct parameters {
        /// Spatial distribution type
        ion_distribution_t ion_distribution{SurfaceRandom};
        /// Atomic number of generated ions. Default 1.
        int ionZ{1}; // atomic number
        /// Mass of generated ions. Default 1.0
        float ionM{1.f}; // ion mass
        /// Initial energy of generated ions [eV]. Default 1000.0 eV
        float ionE0{1000.f}; // initial energy eV
        /// Initial direction vector of generated ions. Default (1,0,0)
        vector3 dir{1,0,0}; // initial direction
        /// Initial position of generated ions [nm]. Used if ion_distribution_t==FixedPos. Default (0,0,0)
        vector3 pos{0,0,0}; // initial position
        //parameters();
    };

protected:

    parameters par_;
    const atom* atom_; // atomic species

public:
    /// Default constructor
    ion_beam();
    /// Copy constructor
    ion_beam(const ion_beam& i);

    /// Set ion_beam class parameter values
    void setParameters(const parameters& p) { par_ = p; }
    /// Returns the ion_beam parameters
    const parameters& getParameters() const { return par_; }

    /// Return the atomic number of generated ions
    int ionZ() const { return par_.ionZ; }
    /// Return the mass of generated ions
    float ionM() const { return par_.ionM; }
    /// Return the energy of generated ions [eV]
    float ionE0() const { return par_.ionE0; }
    /// Return the spatial distribution type of the ion beam
    ion_distribution_t distributionType() const { return par_.ion_distribution; }
    /// Return the direction of generated ions (if it is fixed)
    vector3 ionDir() const { return par_.dir; }
    /// Return the position of generated ions (if it is fixed)
    vector3 ionPos() const { return par_.pos; }

    /// Return the atom class of the projectile ions
    const atom* projectile() const { return atom_; }
    /// Set the atom class of the projectile ions
    void setProjectile(const atom* at) { atom_ = at; }

    /**
     * @brief source_ion generates an ion
     * @tparam _U a random number generator class
     * @param g the random number generator
     * @param t the simulation target
     * @param i the generated ion
     */
    void source_ion(random_vars& g, const target& t, ion& i);

};

#endif // ION_BEAM_H
