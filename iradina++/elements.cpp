#include "elements.h"

#include <cstring>

/*************** Adapted from constants.c ***********************/

/* Extracted from NIST data in May 2007 (physics.nist.gov/PhysRefData/Compositions/index.html)
Developers and Contributors:
J. S. Coursey, D. J. Schwab, and R. A. Dragoset
NIST, Physics Laboratory, Office of Electronic Commerce in Scientific and Engineering Data
(There are 100 data but only 92 are used with SRIM 2006) */

float elements::mostAbundantIsotope(int atomic_num)
{

    float MostAbundantIsotope[]=
        { 92,
         1.0078250321f,    4.0026032497f,       7.0160040f,       9.0121821f,      11.0093055f,
         12.0000000f,   14.0030740052f,   15.9949146221f,     18.99840320f,   19.9924401759f,
         22.98976967f,     23.98504190f,     26.98153844f,   27.9769265327f,     30.97376151f,
         31.97207069f,     34.96885271f,    39.962383123f,      38.9637069f,      39.9625912f,
         44.9559102f,      47.9479471f,      50.9439637f,      51.9405119f,      54.9380496f,
         55.9349421f,      58.9332002f,      57.9353479f,      62.9296011f,      63.9291466f,
         68.925581f,      73.9211782f,      74.9215964f,      79.9165218f,      78.9183376f,
         83.911507f,      84.9117893f,      87.9056143f,      88.9058479f,      89.9047037f,
         92.9063775f,      97.9054078f,       97.907216f,     101.9043495f,      102.905504f,
         105.903483f,      106.905093f,     113.9033581f,      114.903878f,     119.9021966f,
         120.9038180f,     129.9062228f,      126.904468f,     131.9041545f,      132.905447f,
         137.905241f,      138.906348f,      139.905434f,      140.907648f,      141.907719f,
         144.912744f,      151.919728f,      152.921226f,      157.924101f,      158.925343f,
         163.929171f,      164.930319f,      165.930290f,      168.934211f,     173.9388581f,
         174.9407679f,     179.9465488f,      180.947996f,     183.9509326f,     186.9557508f,
         191.961479f,      192.962924f,      194.964774f,      196.966552f,      201.970626f,
         204.974412f,      207.976636f,      208.980383f,      208.982416f,      209.987131f,
         222.0175705f,     223.0197307f,     226.0254026f,     227.0277470f,     232.0380504f,
         231.0358789f,     238.0507826f,     237.0481673f,      244.064198f,     243.0613727f,
         247.070347f,      247.070299f,      251.079580f,      252.082970f,     257.0950990f  };

    return (atomic_num<=max_atomic_num && atomic_num>0) ? MostAbundantIsotope[atomic_num] : 0.f;
}

float elements::mass(int atomic_num)
{
    float AtomicMass[]=
        { 92,
         1.00794f,        4.002602f,           6.941f,        9.012182f,          10.811f,
         12.0107f,         14.0067f,         15.9994f,      18.9984032f,         20.1797f,
         22.989770f,         24.3050f,       26.981538f,         28.0855f,       30.973761f,
         32.065f,          35.453f,          39.948f,         39.0983f,          40.078f,
         44.955910f,          47.867f,         50.9415f,         51.9961f,       54.938049f,
         55.845f,       58.933200f,         58.6934f,          63.546f,          65.409f,
         69.723f,           72.64f,        74.92160f,           78.96f,          79.904f,
         83.798f,         85.4678f,           87.62f,        88.90585f,          91.224f,
         92.90638f,           95.94f,             98.f,          101.07f,       102.90550f,
         106.42f,        107.8682f,         112.411f,         114.818f,         118.710f,
         121.760f,          127.60f,       126.90447f,         131.293f,       132.90545f,
         137.327f,        138.9055f,         140.116f,       140.90765f,          144.24f,
         145.f,          150.36f,         151.964f,          157.25f,       158.92534f,
         162.500f,       164.93032f,         167.259f,       168.93421f,          173.04f,
         174.967f,          178.49f,        180.9479f,          183.84f,         186.207f,
         190.23f,         192.217f,         195.078f,       196.96655f,          200.59f,
         204.3833f,           207.2f,       208.98038f,            209.f,            210.f,
         222.f,            223.f,            226.f,            227.f,        232.0381f,
         231.03588f,       238.02891f,            237.f,            244.f,            243.f,
         247.f,            247.f,            251.f,            252.f,            257.f  };

    return (atomic_num<=max_atomic_num && atomic_num>0) ? AtomicMass[atomic_num] : 0.f;
}

static const char* AtomicNames[] =
    {"X",
     "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
     "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
     "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
     "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
     "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
     "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf"};

const char* elements::name(int atomic_num)
{
    return (atomic_num<=max_atomic_num && atomic_num>0) ? AtomicNames[atomic_num] : AtomicNames[0];
}

int elements::atomicNum(const char* name) {
    int Z = 1;
    while (Z <= max_atomic_num && strcmp(name,AtomicNames[Z])) Z++;
    return (Z <= max_atomic_num) ? Z : -1;
}
