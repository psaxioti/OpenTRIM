#include "simulation.h"
#include "elements.h"
#include "xs.h"
#include "settings.h"

#include <cmath>

#include <cfloat>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <limits>

#include <thread>

int Fe_2MeV_on_Fe();
void print_ini_template();

using std::cout;
using std::cerr;
using std::endl;

typedef float RealType;

#include <bitset>

typedef corteo_index<4,0,2> tsti;

int parse(std::string& fname, simulation_base* &S);

int test0(int argc, char* argv[]);
int test1(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    return test1(argc,argv);
}

int test0(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.ini]" << endl;
        return -1;
    }

    settings s;
    if (s.parse(argv[1],true)!=0) return -1;

    cout << endl << endl << "Starting simulation ..." << endl << endl;

    controller C;
    C.run_simulation(s);

    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << C.ips << " (total), ";
    cout << C.ips/s.psim_.threads << " (per thread)" << endl;
    cout << "Ions = " << s.psim_.max_no_ions << endl;


    return 0;
}

int test1(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.ini]" << endl;
        return -1;
    }

    settings s;
    if (s.parse(argv[1],true)!=0) return -1;

    cout << endl << endl << "Starting simulation ..." << endl << endl;

    simulation_base* S = s.createSimulation();
    S->init();
    S->run();
    S->saveTallys();

    tally t = S->getTally();
    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << S->ips() << " (total), ";
    cout << S->ips()/S->getParameters().threads << " (per thread)" << endl;
    cout << "Ions = " << t.Nions() << endl;


    return 0;
}



int Fe_2MeV_on_Fe()
{
    simulation_base::parameters p;
    p.max_no_ions = 100;
    p.random_generator_type = simulation_base::MinStd;
    p.random_var_type = simulation_base::Tabulated;

    int Zfe = elements::atomicNum("Fe");
    float Mfe = elements::mass(Zfe);

    ion_beam::parameters bp;
    bp.ionZ_ = Zfe;
    bp.ionM_ = Mfe;
    bp.ionE0_ = 2e6f;

    simulation_base* S = simulation_base::fromParameters(p);

    S->setIonBeam(bp);

    material* Fe = S->addMaterial("Fe");
    Fe->setMassDensity(7.8658f);
    Fe->addAtom({.Z=Zfe, .M=Mfe, .Ed=40.f, .El=1.f, .Es=1.f, .Er=40.f}, 1.f);

    grid3D& g = S->grid();
    float L = 1200;
    g.setX(0,L,100);
    g.setY(0,L,1,true);
    g.setZ(0,L,1,true);

    box3D box = g.box();
    S->fill(box,Fe);

    S->init();
    S->run();

    tally t = S->getTally();

    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << S->ips() << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Nrecoils()/t.Npkas() << endl;

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
