#ifndef TALLY_H
#define TALLY_H

#include "arrays.h"
#include "ion.h"
#include "target.h"

enum class Event : uint32_t {
    NewSourceIon = 1 << 0,
    NewRecoil    = 1 << 1,
    Scattering   = 1 << 2,
    IonExit      = 1 << 3,
    IonStop      = 1 << 4,
    BoundaryCrossing = 1 << 5,
    Ioniz        = 1 << 6,
    Phonon       = 1 << 7,
    Replacement  = 1 << 8,
    Vacancy      = 1 << 9,
    CascadeComplete = 1 << 10,
    NewFlightPath  = 1 << 11,
    NEvent         = 1 << 12
};

class tally {

public:

    static const int std_tallies = 16;
    static const char* arrayName(int i);
    static const char* arrayDescription(int i);
    static const char* arrayGroup(int i);

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

    std::array< ArrayNDd, std_tallies > A;

    enum counter_t {
        H = 0,
        V = 1,
        I = 2,
        R = 3,
        P = 4,
        L = 5,
        D = 6
    };

    enum tally_t {
        cT = 0,
        cV = 1,
        cI = 2,
        cR = 3,
        cP = 4,
        cL = 5,
        eIoniz = 6,
        ePhonon = 7,
        ePKA = 8,
        eLost = 9,
        dpTdam = 10,
        dpTdam_LSS = 11,
        dpVnrt = 12,
        dpVnrt_LSS = 13,
        isFlightPath = 14,
        isCollision = 15
    };

    uint32_t EventMask_ { static_cast<uint32_t>(Event::NEvent) - 1 };

public:

    const double& Nions() const { return A[cT](H); }
    const double& Npkas() const { return A[cT](P); }
    const double& Ndisp() const { return A[cT](D); }
    const double& Nrepl() const { return A[cT](R); }
    const double& Nimpl() const { return A[cT](I); }
    const double& Nvac()  const { return A[cT](V); }
    const double& Nlost() const { return A[cT](L); }

    const ArrayNDd& vacancies() const { return A[cV]; }
    const ArrayNDd& implantations() const { return A[cI]; }
    const ArrayNDd& replacements() const { return A[cR]; }
    const ArrayNDd& pkas() const { return A[cP]; }
    const ArrayNDd& lost() const { return A[cL]; }

    const ArrayNDd& ionization() const { return A[eIoniz]; }
    const ArrayNDd& phonons() const { return A[ePhonon]; }
    const ArrayNDd& lostE() const { return A[eLost]; }
    const ArrayNDd& pkaE() const { return A[ePKA]; }

    const ArrayNDd& Tdam() const { return A[dpTdam]; }
    const ArrayNDd& Tdam_LSS() const { return A[dpTdam_LSS]; }
    const ArrayNDd& Vnrt() const { return A[dpVnrt]; }
    const ArrayNDd& Vnrt_LSS() const { return A[dpVnrt_LSS]; }

    const ArrayNDd& at(int i) const { return A[i]; }

    int init(int natoms, int ncells) {
        A[0] = ArrayNDd(7);
        for(int i=1; i<std_tallies; i++)
            A[i] = ArrayNDd(natoms,ncells);
        return 0;
    }

    void clear() {
        for(int i=0; i<std_tallies; i++)
            A[i].clear();
    }

    tally& operator+=(const tally& t) {
        for(int i=0; i<A.size(); i++)
            A[i] += t.A[i];
        return *this;
    }

    void addSquared(const tally& t) {
        for(int i=0; i<A.size(); i++)
            A[i].addSquared(t.A[i]);
    }

    void copy(const tally& t) {
        for(int i=0; i<A.size(); i++)
            A[i] = t.A[i].copy();
    }

    inline void operator()(Event ev, const ion& i, const void* pv = 0)
    {
        int iid = i.myAtom()->id(), cid=i.cellid(), pid=i.prev_cellid();
        const float* p;
        const atom* a;
        switch (ev) {
        case Event::BoundaryCrossing:
            A[isCollision](iid,pid) += i.ncoll();
            A[isFlightPath](iid,pid) += i.path();
            A[ePhonon](iid,pid) += i.phonon();
            A[eIoniz](iid,pid) += i.ioniz();
            break;
        case Event::Replacement:
            // add a replacement
            A[cT](R)++;
            A[cR](iid,cid)++; // this atom, current cell
            // Remove a Vac
            A[cT](V)--;
            // for the replaced atom (id passed in pointer pv), in current cell
            a = reinterpret_cast<const atom*>(pv);
            A[cV](a->id(),cid)--;
            // if this was a recoil, add a V at the starting cell
            if (i.recoil_id()) {
                A[cT](D)++;
                A[cT](V)++;
                A[cV](iid,i.cellid0())++;
            } else A[cT](H)++;
            A[isCollision](iid,cid) += i.ncoll();
            A[isFlightPath](iid,cid) += i.path();
            A[eIoniz](iid,cid) += i.ioniz();
            A[ePhonon](iid,cid) += i.erg() + i.phonon();
            break;
        case Event::IonStop:
            A[cT](I)++;
            A[cI](iid,cid)++;
            if (i.recoil_id()) {
                A[cT](D)++;
                A[cT](V)++;
                A[cV](iid,i.cellid0())++;
            }  else A[cT](H)++;
            A[isCollision](iid,cid) += i.ncoll();
            A[isFlightPath](iid,cid) += i.path();
            A[eIoniz](iid,cid) += i.ioniz();
            A[ePhonon](iid,cid) += i.erg() + i.phonon();
            break;
        case Event::IonExit:
            A[cT](L)++;
            A[cL](iid,pid)++;
            if (i.recoil_id()) {
                A[cT](D)++;
                A[cT](V)++;
                A[cV](iid,i.cellid0())++;
            } else A[cT](H)++;
            A[isCollision](iid,pid) += i.ncoll();
            A[isFlightPath](iid,pid) += i.path();
            A[eIoniz](iid,pid) += i.ioniz();
            A[ePhonon](iid,pid) += i.phonon();
            A[eLost](iid,pid) += i.erg();
            break;
        case Event::CascadeComplete:
            A[cT](P) += 1;
            A[cP](iid,cid)++;
            p = reinterpret_cast<const float*>(pv);
            A[ePKA](iid,cid) += p[0];
            A[dpTdam_LSS](iid,cid) += p[1];
            A[dpVnrt_LSS](iid,cid) += p[2];
            A[dpTdam](iid,cid) += p[3];
            A[dpVnrt](iid,cid) += p[4];
            break;
        default:
            break;
        }

    }


};



#endif // TALLY_H
