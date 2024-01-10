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

typedef float RealType;

int main(int argc, char* argv[])
{

    std::random_device rd;
    std::mt19937 gen(rd());

    const std::size_t Bits = 10;
    int b = std::min(Bits, std::size_t {std::numeric_limits<RealType>::digits});
    cout << b << endl << endl;

    for (int n = 0; n < 10; ++n)
        cout << std::generate_canonical<double, Bits>(gen) << endl;


//    unsigned long w1(2147483647UL), w2(4294967295UL);
//    unsigned long i;
//    cout << "i\t31bitS\t31bitD\t32bitS\t32bitD" << endl;
//    for(i=0; i<1024; i++) {

//        cout << i << '\t'
//             << 1.f*(w1-i)/w1 - 1.f << '\t'
//             << 1.*(w1-i)/w1 -1. << '\t'
//             << 1.f*(w2-i)/w2 - 1.f << '\t'
//             << 1.*(w2-i)/w2 -1.<< endl;
//    }
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

int test_run()
{
    //SimZBLMagic S("Test");
    SimCorteo6bit S("Test");

    material* Fe = S.addMaterial("Fe", 7.8658);
    int Zfe = elements::atomicNum("Fe");
    float Mfe = elements::mass(Zfe);
    Fe->addAtom(Zfe, Mfe, 1.f, 40.f,
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

    S.setProjectile(Zfe, Mfe, 2E6);
    S.setStragglingModel(simulation_base::NoStraggling);

    grid3D& g = S.grid();
    float L = 1200;
    g.setX(0,L,100);
    g.setY(0,L,1,true);
    g.setZ(0,L,1,true);

    box3D box = g.box();
    S.fill(box,Fe);
    // auto& X0 = box.min();
    // X0.x() = 500;
    // S.fill(box,Si);

    S.init();
    S.run(100,"testrun.h5");

    cout << "ms/ion =" << S.ms_per_ion() << endl;

    return 0;
}

void test_dedx()
{
    const float *dedx, *de_stragg;

    SimCorteo4bit S("Test");

    S.setProjectile(1, elements::mostAbundantIsotope(1), 5E6);

    material* Fe = S.addMaterial("FeCr", 7.8658);
    int Zfe = elements::atomicNum("Fe");
    Fe->addAtom(Zfe, elements::mass(Zfe), 1.f, 40.f,
                1.f, 1.f, 40.f);


    cout << "H in Fe" << endl;
    S.getDEtables(S.getSource().projectile(), Fe, dedx, de_stragg);
    for(dedx_index ie; ie!=ie.end(); ie++)
        cout << *ie << '\t' << dedx[ie] << '\t' << de_stragg[ie] << endl;

    cout << endl << "Fe in Fe" << endl;
    S.getDEtables(Fe->atoms().front(), Fe, dedx, de_stragg);
    for(dedx_index ie; ie!=ie.end(); ie++)
        cout << *ie << '\t' << dedx[ie] << '\t' << de_stragg[ie] << endl;

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
