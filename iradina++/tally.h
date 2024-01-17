#ifndef TALLY_H
#define TALLY_H

#include "arrays.h"

class tally {

    std::array< unsigned int, 6 > total_counters;
    std::array< Array2Dui, 3 > cell_counters;
    std::array< Array2Dd, 2 > energy_counters;

    typedef enum {
        Vac = 0,
        Impl,
        Repl,
        Rec,
        Pka,
        Ion
    } counter_t;

    typedef enum {
        Ioniz = 0,
        Phonon,
        Recoil
    } energy_t;



public:

    typedef unsigned int uint;

    uint Nions() const { return total_counters[Ion]; }
    uint& Nions() { return total_counters[Ion]; }
    uint Npkas() const { return total_counters[Pka]; }
    uint& Npkas() { return total_counters[Pka]; }
    uint Nrecoils() const { return total_counters[Rec]; }
    uint& Nrecoils() { return total_counters[Rec]; }
    uint Nrepl() const { return total_counters[Repl]; }
    uint& Nrepl() { return total_counters[Repl]; }
    uint Nimpl() const { return total_counters[Impl]; }
    uint& Nimpl() { return total_counters[Impl]; }
    uint Nvac() const { return total_counters[Vac]; }
    uint& Nvac() { return total_counters[Vac]; }

    uint& vacancies(int i, int j)
    { return cell_counters[Vac](i,j); }
    const uint& vacancies(int i, int j) const
    { return cell_counters[Vac](i,j); }
    Array2Dui vacancies() const
    { return cell_counters[Vac]; }

    uint& implantations(int i, int j)
    { return cell_counters[Impl](i,j); }
    const uint& implantations(int i, int j) const
    { return cell_counters[Impl](i,j); }
    Array2Dui implantations() const
    { return cell_counters[Impl]; }

    uint& replacements(int i, int j)
    { return cell_counters[Repl](i,j); }
    const uint& replacements(int i, int j) const
    { return cell_counters[Repl](i,j); }
    Array2Dui replacements() const
    { return cell_counters[Repl]; }

    double& ionization(int i, int j)
    { return energy_counters[Ioniz](i,j); }
    const double& ionization(int i, int j) const
    { return energy_counters[Ioniz](i,j); }
    Array2Dd ionization() const
    { return energy_counters[Vac]; }

    double& phonons(int i, int j)
    { return energy_counters[Phonon](i,j); }
    const double& phonons(int i, int j) const
    { return energy_counters[Phonon](i,j); }
    Array2Dd phonons() const
    { return energy_counters[Phonon]; }

    int init(int natoms, int ncells) {
        total_counters.fill(0);
        for(int i=0; i<cell_counters.size(); i++)
            cell_counters[i] = Array2Dui(natoms,ncells);
        for(int i=0; i<energy_counters.size(); i++)
            energy_counters[i] = Array2Dd(natoms,ncells);
        return 0;
    }

    tally& operator+=(const tally& t) {
        for(int i=0; i<total_counters.size(); i++)
            total_counters[i] += t.total_counters[i];
        for(int i=0; i<cell_counters.size(); i++)
            cell_counters[i] += t.cell_counters[i];
        for(int i=0; i<energy_counters.size(); i++)
            energy_counters[i] += t.energy_counters[i];
        return *this;
    }

    void copy(const tally& t) {
        for(int i=0; i<total_counters.size(); i++)
            total_counters[i] = t.total_counters[i];
        for(int i=0; i<cell_counters.size(); i++)
            cell_counters[i] = t.cell_counters[i].copy();
        for(int i=0; i<energy_counters.size(); i++)
            energy_counters[i] = t.energy_counters[i].copy();
    }


};



#endif // TALLY_H
