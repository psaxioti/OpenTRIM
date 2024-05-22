#ifndef TALLY_H
#define TALLY_H

#include "arrays.h"
#include "ion.h"
#include "target.h"

enum class Event {
    NewSourceIon = 0,
    NewRecoil,
    Scattering,
    IonExit,
    IonStop,
    Ioniz,
    Phonon,
    Replacement,
    Vacancy,
    NRT_LSS_damage,
    NRT_damage,
    NewFP
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
        csFP = 14,
        csCollision = 15
    };

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

    inline void operator()(Event ev, const ion& i, const float* p = 0)
    {
        switch (ev) {
        case Event::NewSourceIon:
            A[cT](H)++;
            break;
        case Event::NewRecoil:
            A[cT](D)++;
            if (i.recoil_id()==1) { // PKA
                A[cT](P) += 1;
                A[cP](i.myAtom()->id(),i.cellid())++;
                A[ePKA](i.myAtom()->id(),i.cellid()) += i.erg();
            }
            break;
        case Event::Replacement:
            A[cT](R)++;
            A[cR](i.myAtom()->id(),i.cellid())++;
            A[ePhonon](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::Vacancy:
            A[cT](V)++;
            A[cV](i.myAtom()->id(),i.cellid())++;
            break;
        case Event::IonStop:
            A[cT](I)++;
            A[cI](i.myAtom()->id(),i.cellid())++;
            A[ePhonon](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::IonExit:
            A[cT](L)++;
            A[cL](i.myAtom()->id(),i.cellid())++;
            A[eLost](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::Ioniz:
            A[eIoniz](i.myAtom()->id(),i.cellid()) += p[0];
            break;
        case Event::Phonon:
            A[ePhonon](i.myAtom()->id(),i.cellid()) += p[0];
            break;
        case Event::NRT_LSS_damage:
            A[dpTdam_LSS](i.myAtom()->id(),i.cellid()) += p[0];
            A[dpVnrt_LSS](i.myAtom()->id(),i.cellid()) += p[1];
            break;
        case Event::NRT_damage:
            A[dpTdam](i.myAtom()->id(),i.cellid()) += p[0];
            A[dpVnrt](i.myAtom()->id(),i.cellid()) += p[1];
            break;
        case Event::NewFP:
            A[csFP](i.myAtom()->id(),i.cellid()) += p[0];
            break;
        case Event::Scattering:
            A[csCollision](i.myAtom()->id(),i.cellid())++;
            break;
        default:
            break;
        }

    }


};



#endif // TALLY_H
