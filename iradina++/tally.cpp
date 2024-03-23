#include "tally.h"

#include "ion.h"
#include "target.h"

const char* tally::arrayName(int i)
{
    static const char* names[] = {
        "Totals",
        "Vacancies",
        "Implantations",
        "Replacements",
        "Lost",
        "IonizationEnergy",
        "PhononEnergy",
        "LostEnergy",
        "Tdam",
        "Tdam_LSS",
        "Vnrt",
        "Vnrt_LSS",
        "X"
    };

    return (i < std_tallies) ? names[i] : names[std_tallies];



}


