#include "simulation.h"
#include "elements.h"
#include "xs.h"
#include "dedx.h"

#include <cmath>

#include <cfloat>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <limits>

#include <inicpp.h>

void write_corteo();
void write_corteo_comp(int n);
void write_corteo_tst4();
int test_run();

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{

    return test_run();
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

int test_run()
{
    //reducedXS_quad xs;
    //reducedXS_corteo4bit* cXS = new reducedXS_corteo4bit(xs);
    reducedXS cXS(reducedXS::ZBL, reducedXS::Corteo4bitTable);
    cXS.init(&std::cout);

    simulation S("Test",cXS);

    material* Fe = S.addMaterial("FeCr", 7.8658);
    int Zfe = elements::atomicNum("Fe");
    Fe->addAtom(Zfe, elements::mass(Zfe), 1.f, 40.f,
                  1.f, 1.f, 40.f);

//    material* FeCr = S.addMaterial("FeCr", 7.8658);
//    int Zfe = elements::atomicNum("Fe");
//    FeCr->addAtom(Zfe, elements::mass(Zfe), 1.f, 40.f,
//               1.f, 1.f, 40.f);
//    int Zcr = elements::atomicNum("Cr");
//    FeCr->addAtom(Zcr, elements::mass(Zcr), 1.f, 40.f,
//                1.f, 1.f, 40.f);
//    material* Si = S.addMaterial("Si", 7.8);
//    int Zsi = elements::atomicNum("Si");
//    Si->addAtom(Zsi, elements::mass(Zsi), 1.f, 25.f,
//                1.f, 1.f, 25.f);

    S.setProjectile(1, elements::mostAbundantIsotope(1), 5E6);

    grid3D& g = S.grid();
    g.setX(0,100000,100);
    g.setY(0,100000,1,true);
    g.setZ(0,100000,1,true);

    box3D box = g.box();
    S.fill(box,Fe);
    // auto& X0 = box.min();
    // X0.x() = 500;
    // S.fill(box,Si);

    S.init();
    S.run(100,"testrun.h5");

    return 0;
}

void test_dedx()
{
    const float *dedx, *de_stragg;

    simulation S("Test");

    S.setProjectile(1, elements::mostAbundantIsotope(1), 5E6);

    material* Fe = S.addMaterial("FeCr", 7.8658);
    int Zfe = elements::atomicNum("Fe");
    Fe->addAtom(Zfe, elements::mass(Zfe), 1.f, 40.f,
                1.f, 1.f, 40.f);


    cout << "H in Fe" << endl;
    S.getInventory().getDEtables(S.getSource().projectile(), Fe, dedx, de_stragg);
    for(dedx_index ie; ie!=ie.end(); ie++)
        cout << *ie << '\t' << dedx[ie] << '\t' << de_stragg[ie] << endl;

    cout << endl << "Fe in Fe" << endl;
    S.getInventory().getDEtables(Fe->atoms().front(), Fe, dedx, de_stragg);
    for(dedx_index ie; ie!=ie.end(); ie++)
        cout << *ie << '\t' << dedx[ie] << '\t' << de_stragg[ie] << endl;

}

void write_corteo_tst4()
{
    reducedXS xs(reducedXS::ZBL, reducedXS::GaussMehlerQuad);
    reducedXS cxs(reducedXS::ZBL, reducedXS::Corteo4bitTable);
    reducedXS mg(reducedXS::ZBL, reducedXS::ZBL_Magick);
    cxs.init(&std::cout);

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
