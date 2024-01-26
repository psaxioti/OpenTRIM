#include "simulation.h"
#include "elements.h"
#include "settings.h"

#include <iostream>
#include <fstream>

int Fe_2MeV_on_Fe();

using std::cout;
using std::cerr;
using std::endl;

int test0(int argc, char* argv[]);
int Fe_2MeV_on_Fe();

void my_progress_callback(const std::vector<uint>& v)
{
    static int i = 0;

    if (i==0) {
        for(int i=0; i<v.size()-1; i++)
            cout << "#" << i+1 << '\t';
        cout << "Total" << endl;
        for(int i=0; i<v.size()-1; i++)
            cout << "-------" << '\t';
        cout << "-------" << endl;
        i++;
    }
    for(int i=0; i<v.size(); i++)
        cout << v[i] << '\t';
    cout << endl;
}

int main(int argc, char* argv[])
{
    return test0(argc,argv);
    //return Fe_2MeV_on_Fe();
    //return settings::test();
}

int test0(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.ini]" << endl;
        return -1;
    }

    cout << "Parsing config file ... ";
    std::ifstream is(argv[1]);
    json j = json::parse(is);
    settings s = j.template get<settings>();

//    if (s.parse(is)!=0) {
//        cout << "error!" << endl;
//        return -1;
//    }

    cout << "done." << endl << endl;
    cout << "Starting simulation ..." << endl << endl;

    simulation_base* S = s.createSimulation();
    S->init();
    S->exec(my_progress_callback);
    S->saveTallys();

    tally t = S->getTally();
    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << S->ips() << " (total), ";
    cout << S->ips()/S->nThreads() << " (per thread)" << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Nrecoils()/t.Npkas() << endl;

    return 0;
}

int Fe_2MeV_on_Fe()
{
    simulation<XS_corteo4bit, std::mt19937> S("2MeV Fe on Fe");

    S.setRandomVarType(simulation_base::Sampled);
    S.setMaxIons(10);
    S.setNThreads(1);

    int Zfe = elements::atomicNum("Fe");
    float Mfe = elements::mass(Zfe);
    S.setIonBeam({.ionZ = Zfe,
                  .ionM = Mfe,
                  .ionE0 = 2e6});


    material* Fe = S.addMaterial("Fe");
    Fe->setMassDensity(7.8658f);
    Fe->addAtom({.Z=Zfe, .M=Mfe, .Ed=40.f, .El=0.f, .Es=0.f, .Er=40.f}, 1.f);

    grid3D& g = S.grid();
    float L = 1200;
    g.setX(0,L,100);
    g.setY(0,L,1,true);
    g.setZ(0,L,1,true);

    box3D box = g.box();
    S.fill(box,Fe);

    S.init();
    S.exec(my_progress_callback);
    S.saveTallys();

    tally t = S.getTally();
    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << S.ips() << " (total), ";
    cout << S.ips()/S.nThreads() << " (per thread)" << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Nrecoils()/t.Npkas() << endl;

    return 0;
}


