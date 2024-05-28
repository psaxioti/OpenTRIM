#include "tally.h"

const char* tally::arrayName(int i)
{
    static const char* names[] = {
        "Totals",
        "Vacancies",
        "Implantations",
        "Replacements",
        "PKAs",
        "Lost",
        "Ionization",
        "Phonons",
        "PKA",
        "Lost",
        "Tdam",
        "Tdam_LSS",
        "Vnrt",
        "Vnrt_LSS",
        "flight_path",
        "collisions",
        "X"
    };

    return (i < std_tallies) ? names[i] : names[std_tallies];
}

const char* tally::arrayDescription(int i)
{
    static const char* desc[] = {
        "Totals counts (histories, V, I, R, L, D, PKA)",
        "Vacancies",
        "Implantations",
        "Replacements",
        "PKAs",
        "Ions that exit the simulation volume",
        "Energy deposited to electron ionization [eV]",
        "Energy deposited to the lattice [eV]",
        "PKA recoil energy [eV]",
        "Energy lost due to ions exiting the simulation [eV]",
        "Damage energy  [eV]",
        "Damage energy estimated by the LSS approximation [eV]",
        "Vacancies per the NRT model using Tdam",
        "Vacancies per the NRT model using Tdam_LSS",
        "flight path [nm]",
        "ion collisions",
        "X"
    };

    return (i < std_tallies) ? desc[i] : desc[std_tallies];
}

const char* tally::arrayGroup(int i)
{
    static const char* desc[] = {
        "totals",
        "defects",
        "defects",
        "defects",
        "defects",
        "defects",
        "energy_deposition",
        "energy_deposition",
        "energy_deposition",
        "energy_deposition",
        "damage",
        "damage",
        "damage",
        "damage",
        "ion_stat",
        "ion_stat",
        "X"
    };

    return (i < std_tallies) ? desc[i] : desc[std_tallies];
}


