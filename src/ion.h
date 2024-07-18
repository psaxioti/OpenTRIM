#ifndef _ION_H_
#define _ION_H_

#include <cmath>
#include <queue>

#include "geometry.h"

class atom;

/**
 * \defgroup Ions Ion generation & transport
 *
 * @brief Classes for performing ion transport and injecting new ions into the simulation.
 *
 * @{
 *
 * @ingroup MC
 *
 * @}
 *
 *
 *
 */


/**
 * @brief Enum characterizing the type of boundary crossing for an ion
 * 
 * @ingroup Ions
 */
enum class BoundaryCrossing {
    None,       /**< No boundary crossing occured. */
    Internal,   /**< The ion crossed an internal cell boundary. The ion changes cell */ 
    External    /**< The ion crossed an external boundary of the simulation volume. The ion exits the simulation */ 
};

/**
 * @brief The ion class represents a moving ion in the simulation.
 *
 * The class stores the ion data, such as the position and direction
 * vectors, the energy, the atomic species,
 * the id of the geometric cell the ion is currently in, the id of the
 * parent source ion, and the recoil generation id of the current ion.
 *
 * Furthermore, here are defined 2 very important simulation functions:
 * \ref propagate(), which advances the position of the moving ion and
 * \ref deflect(), which changes the ion's direction after scattering.
 *
 * @ingroup Ions
 */
class ion
{
    vector3 pos_; // position = x,y,z in nm
    vector3 dir_; // direction cosines
    float erg_; // energy in eV
    ivector3 icell_;
    int cellid_, prev_cellid_;
    int ion_id_;
    int recoil_id_;
    const atom* atom_;
    const grid3D* grid_;

    // counters
    uint vac_,impl_,repl_,ncoll_;
    float path_,ioniz_,phonon_,recoil_;

public:

    /// Default constructor
    ion() :
        pos_(0.f,0.f,0.f),
        dir_(0.f,0.f,1.f),
        erg_(1.f),
        icell_(), cellid_(-1), prev_cellid_(-1),
        ion_id_(0), recoil_id_(0),
        atom_(nullptr), grid_(nullptr),
        vac_(0),impl_(0),repl_(0),ncoll_(0),
        path_(0),ioniz_(0),phonon_(0), recoil_(0)
    {}

    /// Returns the ion's position vector [nm]
    const vector3& pos() const { return pos_; }

    /// Returns the vector of the ion's direction cosines
    const vector3& dir() const { return dir_; }

    /// Returns the ion's kinetic energy
    float erg() const { return erg_; }

    /// Returns the index vector of the cell the ion is currently in
    const ivector3& icell() const { return icell_; }

    /// Returns the id of the cell the ion is currently in
    int cellid() const { return cellid_; }

    /// Returns the id of the cell the ion was previously in
    int prev_cellid() const { return prev_cellid_; }

    /// Returns the history id that the current ion belongs to
    int ion_id() const { return ion_id_; }

    /**
     * @brief Returns the recoil generation id of the current ion
     *
     * Source ions have a recoil_id equal to 0, PKAs have recoil_id==1, etc.
     *
     * @return The current ion's recoil_id
     */
    int recoil_id() const { return recoil_id_; }

    /// Returns a pointer to the \ref atom class describing the atomic species of the current ion
    const atom* myAtom() const { return atom_; }

    void de_phonon(float de) {
        erg_ -= de; phonon_ += de;
        assert(erg_>=0);
        assert(finite(erg_));
    }
    void de_ioniz(float de) {
        erg_ -= de; ioniz_ += de;
        assert(erg_>0);
        assert(finite(erg_));
    }
    void de_recoil(float de) {
        erg_ -= de;
        recoil_ += de;
        assert(erg_>=0);
        assert(finite(erg_));
    }
    float phonon() const { return phonon_; }
    float ioniz() const { return ioniz_; }
    float recoil() const { return recoil_; }
    float path() const { return path_; }
    uint ncoll() const { return ncoll_; }

    void add_coll() { ncoll_++; }
    void add_vac() { vac_++; }
    void add_repl() { repl_++; }
    void add_impl() { impl_++; }

    /// Return reference to the vector of direction cosines
    vector3& dir() { return dir_; }   

    // int& ion_id() { return ion_id_; }
    // int& recoil_id() { return recoil_id_; }
    // int& cellid() { return cellid_; }
    // const atom*& myAtom() { return atom_; }
    void setDir(const vector3 d) {
        dir_ = d;
        dir_.normalize();
        assert(dir_.allFinite());
        //float dn = std::abs(dir_.norm()-1.f);
        //assert(dn==0.f);
    }
    int setPos(const vector3& x);
    void setAtom(const atom* a) {
        atom_ = a;
    }
    void setErg(float e) {
        erg_ = e;
        assert(finite(erg_));
        assert(erg_>0);
    }
    void incRecoilId() {
        recoil_id_++;
    }
    void setId(int id) {
        ion_id_ = id;
    }
    void resetRecoilId() {
        recoil_id_ = 0;
    }
    void setGrid(const grid3D* g) {
        grid_ = g;
    }

    /**
     * @brief Deflect the ion after scattering.
     * 
     * Changes the direction after scattering at angles \f$ (\theta,\phi) \f$.
     * 
     * @param n the vector \f$ \mathbf{n} = (\cos\phi\,\sin\theta, \sin\phi\,\sin\theta,\cos\theta) \f$
     * @see  \ref deflect_vector()
     */
    void deflect(const vector3& n) { deflect_vector(dir_,n); }

    void reset_counters();


    BoundaryCrossing propagate(float& s);
};

/**
 * @brief The ion_queue class handles queues of ion objects waiting to be simulated.
 *
 * There are 2 queues, one for the primary recoils (PKAs) and one for
 * all higher recoil generations.
 *
 * When a simulated ion generates primary recoils, these are put on
 * the PKA queue.
 *
 * After the program finishes the ion history, it runs the PKA recoils from the
 * queue. These may produce new recoils, which are put on the secondary recoil queue.
 * The program proceeds to process all recoils in both queues until there is no new recoil.
 *
 * Then the program moves to the next ion.
 *
 * Objects of the ion class are allocated by ion_queue and provided to the simulation
 * by the \ref new_ion() functions.
 *
 * When the simulation finishes an ion history, it returns the ion object to the
 * ion_queue by calling \ref free_ion(). The object buffer can then be used for
 * another ion history.
 *
 * @ingroup Ions
 */
class ion_queue {

    // FIFO ion buffer
    typedef std::queue< ion* > ion_queue_t;

    ion_queue_t ion_buffer_; // buffer of allocated ion objects
    ion_queue_t recoil_queue_; // queue of generated recoils
    ion_queue_t pka_queue_; // queue of generated PKAs

    // pop an ion from the respective queue
    static ion* pop_one_(ion_queue_t& Q)
    {
        if (Q.empty()) return nullptr;
        ion* i = Q.front();
        Q.pop();
        return i;
    }

public:
    /// Returns a new ion object
    ion* new_ion() {
        ion* i;
        if (ion_buffer_.empty()) i = new ion;
        else {
            i = ion_buffer_.front();
            ion_buffer_.pop();
        }
        return i;
    }
    /// Returns a new ion object, initialized with data copied from p
    ion* new_ion(const ion& p) {
        ion* i;
        if (ion_buffer_.empty()) i = new ion(p);
        else {
            i = ion_buffer_.front();
            ion_buffer_.pop();
            *i = p;
        }
        return i;
    }

    /// Push an ion object to the PKA queue
    void push_pka(ion* i) { pka_queue_.push(i); }
    /// Push an ion object to the recoil queue
    void push_recoil(ion* i) { recoil_queue_.push(i); }
    /// Pop a PKA ion object from the queue. If the queue is empty, a nullptr is returned.
    ion* pop_pka() { return pop_one_(pka_queue_); }
    /// Pop a recoil ion object from the queue. If the queue is empty, a nullptr is returned.
    ion* pop_recoil() { return pop_one_(recoil_queue_); }

    /// Release a used ion object
    void free_ion(ion* i) { ion_buffer_.push(i); }

    /// Clear all allocated ion objects from memory
    void clear() {
        while (!ion_buffer_.empty()) {
            ion* i = ion_buffer_.front();
            ion_buffer_.pop();
            delete i;
        }
    }


};

#endif // ION_H
