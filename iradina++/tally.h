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
    NRT_damage
};

class tally {



    /*
     * Totals: 1 (N-ions,V-vacancies,I-interst/implant,R-replace,L-lost,D-displ,P-pka) 7 x 1
     *
     * Cell counters: 4 (V,I,R,E) atoms x cells
     *
     * Energy : 3 (Ioniz, Phonon, Lost) atoms x cells
     *
     * Damage : 4 (Tdam, Tdam_LSS, Vnrt, Vnrt_LSS) cells?? or atom x cells
     *
     */

    std::array< ArrayNDd, 1+4+3+4 > A;

    enum counter_t {
        N = 0,
        V = 1,
        I = 2,
        R = 3,
        L = 4,
        D = 5,
        P= 6
    };

    enum energy_t {
        Ioniz = 5,
        Phonon = 6,
        LostE = 7
    };

    enum dmg_par_t {
        dpTdam = 8,
        dpTdam_LSS = 9,
        dpVnrt = 10,
        dpVnrt_LSS = 11
    };

public:

    static const int std_tallies = 12;
    static const char* arrayName(int i);
    static const char* arrayDescription(int i);
    static const char* arrayGroup(int i);

    const double& Nions() const { return A[0](N); }
    const double& Npkas() const { return A[0](P); }
    const double& Ndisp() const { return A[0](D); }
    const double& Nrepl() const { return A[0](R); }
    const double& Nimpl() const { return A[0](I); }
    const double& Nvac()  const { return A[0](V); }
    const double& Nlost() const { return A[0](L); }

    const ArrayNDd& vacancies() const { return A[V]; }
    const ArrayNDd& implantations() const { return A[I]; }
    const ArrayNDd& replacements() const { return A[R]; }
    const ArrayNDd& lost() const { return A[L]; }

    const ArrayNDd& ionization() const { return A[Ioniz]; }
    const ArrayNDd& phonons() const { return A[Phonon]; }
    const ArrayNDd& lostE() const { return A[LostE]; }

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
            A[0](N)++;
            break;
        case Event::NewRecoil:
            A[0](D)++;
            A[0](P) += (i.recoil_id()==1);
            break;
        case Event::Replacement:
            A[0](R)++;
            A[R](i.myAtom()->id(),i.cellid())++;
            A[Phonon](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::Vacancy:
            A[0](V)++;
            A[V](i.myAtom()->id(),i.cellid())++;
            break;
        case Event::IonStop:
            A[0](I)++;
            A[I](i.myAtom()->id(),i.cellid())++;
            A[Phonon](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::IonExit:
            A[0](L)++;
            A[L](i.myAtom()->id(),i.cellid())++;
            A[LostE](i.myAtom()->id(),i.cellid()) += i.erg();
            break;
        case Event::Ioniz:
            A[Ioniz](i.myAtom()->id(),i.cellid()) += p[0];
            break;
        case Event::Phonon:
            A[Phonon](i.myAtom()->id(),i.cellid()) += p[0];
            break;
        case Event::NRT_LSS_damage:
            A[dpTdam_LSS](i.myAtom()->id(),i.cellid()) += p[0];
            A[dpVnrt_LSS](i.myAtom()->id(),i.cellid()) += p[1];
            break;
        case Event::NRT_damage:
            A[dpTdam](i.myAtom()->id(),i.cellid()) += p[0];
            A[dpVnrt](i.myAtom()->id(),i.cellid()) += p[1];
            break;
        case Event::Scattering:
        default:
            break;
        }

    }


};



#endif // TALLY_H
