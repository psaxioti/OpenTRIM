#ifndef TALLY_H
#define TALLY_H

#include "arrays.h"
#include "ion.h"

/**
 * \defgroup Tallies Tallies and event streams
 *
 * \brief The results of the Monte-Carlo simulation
 *
 * @ingroup MC
 *
 */

/**
 * @brief The types of Monte-Carlo events
 *
 * @ingroup Tallies
 *
 */
enum class Event : uint32_t {
    NewSourceIon = 1 << 0, /**< A new ion track is started. */
    NewRecoil = 1 << 1, /**< A new recoil track is started. */
    Scattering = 1 << 2, /**< An ion scattering occured. */
    IonExit = 1 << 3, /**< An ion exits the simulation volume. */
    IonStop = 1 << 4, /**< An ion stops inside the simulation volume. */
    BoundaryCrossing = 1 << 5, /**< An ion crosses an internal boundary. */
    Replacement = 1 << 6, /**< A replacement event occurs. */
    Vacancy = 1 << 7, /**< A vacancy is created. */
    CascadeComplete = 1 << 8, /**< A PKA cascade is complete. */
    NewFlightPath = 1 << 9,
    NEvent = 1 << 10
};

/**
 * @brief The tally class obtains and organizes Monte-Carlo scores in tables
 *
 * The \ref tally_t enum enumerates the standard score tables.
 *
 * All tables are 2-dimentional, \f$ [N_{at} \times N_c]\f$, where
 * \f$N_{at}\f$ is the number of atoms in the simulation (including the projectile ion)
 * and \f$N_c\f$ the number of cells, except the very first table (i=0)
 * which stores the total scores for histories, PKAs, etc.
 *
 * Tally scores are stored in double real number format.
 *
 * The element \f$A_k(i,j)\f$ of the \f$k\f$-th table stores the score
 * of the given quantity generated by an atom of id\f$=i\f$ in the
 * \f$j\f$-th cell. E.g., a vacancy of atom with id=2 generated in the cell
 * with index 10 will score +1 in \f$A_1(2,10)\f$.
 *
 * For some tables the projectile (\f$i=0\f$)
 * does not have a meaningfull contribution; e.g., vacancies are
 * only from one of the target atoms, not the projectile.
 * Thus, in this case the 0-th row of the respective table is all zeros.
 * Note that the id of an atom in the simulation
 * is not defined by the atomic species only.
 * E.g., we may have Fe projectiles (i=0) impinging
 * on a target containig Fe2O3 (Fe with id=1) and bulk Fe (Fe with id=3).
 *
 * The tables are divided in groups:
 * - Energy deposition (Ionization, phonons, etc)
 * - Defects (Vacancies, replacements, etc.)
 * - Damage Parameters (Damage energy, NRT vacancies, etc.)
 * - Ion statistics (Flight path, collisions)
 *
 * The grouping is preserved when the tables are saved in the \ref out_file "output archive".
 *
 * In memory, the tally scores add up the contributions from all ion histories.
 * When saved into the \ref out_file "HDF5 output file" all scores are divided by
 * the number of histories, thus, they are normalized "per ion".
 *
 * @ingroup Tallies
 * @see out_file
 *
 */
class tally
{

public:
    /// @brief The type of tally table
    enum tally_t {
        cT = 0, /**< Table of totals (for all the following quantities) */
        cV, /**< Vacancies  */
        cI, /**< Interstitials or Implanted ions  */
        cR, /**< Replacements  */
        cRecombinations, /**< Recombinations */
        cP, /**< PKAs  */
        cL, /**< Lost ions (they exit the simulation)  */
        eIoniz, /**< Ionization energy  */
        eLattice, /**< Lattice energy = thermal energy deposited to the lattice */
        eStored, /**< Energy stored in lattice defects (Frenkel Pairs) */
        eRecoil, /**< Energy transfered to recoil atoms */
        ePKA, /**< PKA recoil energy  */
        eLost, /**< Energy of lost ions  */
        dpTdam, /**< Damage energy  */
        dpTdam_LSS, /**< Damage energy predicted by LSS theory */
        dpVnrt, /**< NRT Vacancies based on Tdam */
        dpVnrt_LSS, /**< NRT Vacancies based on Tdam_LSS */
        isFlightPath, /**< Ion flight path */
        isCollision, /**< Ion collisions */
        tEnd
    };

    /// Number of standard tally tables
    static const int std_tallies = tEnd;
    /// Return the short name of the i-th tally table
    static const char *arrayName(int i);
    /// Return a description for the i-th tally table
    static const char *arrayDescription(int i);
    /// Return the group name for the i-th tally table
    static const char *arrayGroup(int i);

    /*
     * Totals: 1 (N-ions,V-vacancies,I-interst/implant,R-replace,P-pka, L-lost,D-displ) 7 x 1
     *
     * Cell counters: 5 (V,I,R,P,L) atoms x cells
     *
     * Energy : 3 (Ioniz, Phonon, PKA, Lost) atoms x cells
     *
     * Damage : 4 (Tdam, Tdam_LSS, Vnrt, Vnrt_LSS) atom x cells
     *
     */

protected:
    std::array<ArrayNDd, std_tallies> A;

    enum counter_t { H = 0, V = 1, I = 2, R = 3, Rec = 4, P = 5, L = 6 };

    uint32_t EventMask_{ static_cast<uint32_t>(Event::NEvent) - 1 };

public:
    std::vector<std::string> arrayNames() const;

    const ArrayNDd &vacancies() const { return A[cV]; }
    const ArrayNDd &implantations() const { return A[cI]; }
    const ArrayNDd &replacements() const { return A[cR]; }
    const ArrayNDd &recombinations() const { return A[cRecombinations]; }
    const ArrayNDd &pkas() const { return A[cP]; }
    const ArrayNDd &lost() const { return A[cL]; }

    const ArrayNDd &ionization() const { return A[eIoniz]; }
    const ArrayNDd &stored() const { return A[eStored]; }
    const ArrayNDd &lattice() const { return A[eLattice]; }
    const ArrayNDd &lostE() const { return A[eLost]; }
    const ArrayNDd &pkaE() const { return A[ePKA]; }

    const ArrayNDd &Tdam() const { return A[dpTdam]; }
    const ArrayNDd &Tdam_LSS() const { return A[dpTdam_LSS]; }
    const ArrayNDd &Vnrt() const { return A[dpVnrt]; }
    const ArrayNDd &Vnrt_LSS() const { return A[dpVnrt_LSS]; }

    /// Return a const reference to the i-th tally score table
    const ArrayNDd &at(int i) const { return A[i]; }
    /// Return a reference to the i-th tally score table
    ArrayNDd &at(int i) { return A[i]; }

    /// @brief Initialize tally buffers for given # of atoms and cells
    void init(int natoms, int ncells)
    {
        A[0] = ArrayNDd(std_tallies); // the "total sums" for all tallies
        for (int i = 1; i < std_tallies; i++)
            A[i] = ArrayNDd(natoms, ncells);
    }

    /// @brief Zero-out all tally scores
    void clear()
    {
        for (int i = 0; i < std_tallies; i++)
            A[i].clear();
    }

    /// Compute sums of all tallies
    void computeSums()
    {
        A[0][0] = 1; // 1 history
        for (size_t i = 1; i < std_tallies; i++)
            for (size_t j = 0; j < A[i].size(); j++)
                A[0][i] += A[i][j];
    }

    /// @brief Add the scores from another tally
    /// @param t another tally object
    /// @return this object
    tally &operator+=(const tally &t)
    {
        for (int i = 0; i < A.size(); i++)
            A[i] += t.A[i];
        return *this;
    }

    /// @brief Add the squared scores from another tally, i.e., x[i] += x'[i]*x'[i]
    /// @param t another tally object
    void addSquared(const tally &t)
    {
        for (int i = 0; i < A.size(); i++)
            A[i].addSquared(t.A[i]);
    }

    /// @brief Copy contents from another tally
    /// @param t another tally object
    void copy(const tally &t)
    {
        for (int i = 0; i < A.size(); i++)
            A[i] = t.A[i].copy();
    }

    /// @brief Copy contents to another tally
    /// @param t another tally object
    void copyTo(tally &t) const
    {
        for (int i = 0; i < A.size(); i++)
            A[i].copyTo(t.A[i]);
    }

    /// @brief Copy contents to another tally
    /// @param t another tally object
    tally clone() const
    {
        tally t;
        for (int i = 0; i < A.size(); i++)
            t.A[i] = A[i].copy();
        return t;
    }

    /// @brief Score an event
    /// @param ev the event type
    /// @param i the ion causing the event
    /// @param pv pointer to additional data, if available
    void operator()(Event ev, const ion &i, const void *pv = 0);

    bool debugCheck(int id, double E0);
    bool debugCheck(double E0);
    double totalErg(int id);
    double totalErg();
};

#endif // TALLY_H
