#include "simulation.h"
#include "elements.h"
#include "xs.h"

#include <cmath>

#include <cfloat>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <limits>

#include <inicpp.h>

int Fe_2MeV_on_Fe();

using std::cout;
using std::endl;

typedef float RealType;

int main(int argc, char* argv[])
{
    return Fe_2MeV_on_Fe();
}



int Fe_2MeV_on_Fe()
{
    //SimZBLMagic_MSRAND S("Test");
    SimCorteo4bit_MSRAND S("Test");
    //S.setStragglingModel(simulation_base::NoStraggling);
    S.setRandomVarType(simulation_base::Tabulated);

    int Zfe = elements::atomicNum("Fe");
    float Mfe = elements::mass(Zfe);

    material* Fe = S.addMaterial("Fe", 7.8658);
    Fe->addAtom(Zfe, Mfe, 1.f, 40.f,
                  1.f, 1.f, 40.f);
    S.setProjectile(Zfe, Mfe, 2E6);

    grid3D& g = S.grid();
    float L = 1200;
    g.setX(0,L,100);
    g.setY(0,L,1,true);
    g.setZ(0,L,1,true);

    box3D box = g.box();
    S.fill(box,Fe);

    S.init();
    S.run(10000,"testrun.h5");

    cout << "ms/ion = " << S.ms_per_ion() << endl;
    cout << "Ions = " << S.ion_histories() << endl;
    cout << "PKA/Ion = " << 1.f*S.pkas()/S.ion_histories() << endl;
    cout << "Recoils/PKA = " << 1.f*S.recoils()/S.pkas() << endl;

    return 0;
}

int testini(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cerr << "usage: load_ini_file [FILE_PATh]" << std::endl;
        return 1;
    }

    std::string path = argv[1];

    // load the file
    ini::IniFile inif;
    inif.load(path);

    // show the parsed contents of the ini file
    std::cout << "Parsed ini contents" << std::endl;
    std::cout << "Has " << inif.size() << " sections" << std::endl;
    for(const auto &sectionPair : inif)
    {
        const std::string &sectionName = sectionPair.first;
        const ini::IniSection &section = sectionPair.second;
        std::cout << "Section '" << sectionName << "' has " << section.size() << " fields" << std::endl;

        for(const auto &fieldPair : sectionPair.second)
        {
            const std::string &fieldName = fieldPair.first;
            const ini::IniField &field = fieldPair.second;
            std::cout << "  Field '" << fieldName << "' Value '" << field.as<std::string>() << "'" << std::endl;
        }
    }

    return 0;
}



void write_corteo_tst4()
{
    XS_zbl_magic mg;
    XS_corteo4bit cxs;
    XSquad< screeningZBL > xs;

    {
        std::ofstream ostrm1("corteo_tst4_tbl.dat");
        std::ofstream ostrm2("corteo_tst4_quad.dat");
        std::ofstream ostrm3("corteo_tst4_magic.dat");
        float epsilon = 1.e-5f;
        for(int i=0; i<7; i++) {
            float s = 1.e-9f;
            ostrm1 << cxs.sin2Thetaby2(epsilon,s);
            ostrm2 << xs.sin2Thetaby2(epsilon,s);
            ostrm3 << mg.sin2Thetaby2(epsilon,s);
            for(int j=1; j<=1000; j++) {
                s = 50.f*j/1000;
                ostrm1 << '\t' << cxs.sin2Thetaby2(epsilon,s);
                ostrm2 << '\t' << xs.sin2Thetaby2(epsilon,s);
                ostrm3 << '\t' << mg.sin2Thetaby2(epsilon,s);
            }
            ostrm1 << std::endl;
            ostrm2 << std::endl;
            ostrm3 << std::endl;
            epsilon *= 10;
        }
    }
}
