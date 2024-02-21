#include "simulation.h"
#include "elements.h"
#include "options.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;

int Fe_2MeV_on_Fe()
{
    simulation<XS_corteo4bit, std::mt19937> S("2MeV Fe on Fe");

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
    g.setX(0,L,100,false);
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

int genTestJSON()
{
    simulation<XS_corteo4bit, std::mt19937> S("2MeV Fe on Fe");
    S.setOutputFileBaseName("FeFC");
    S.store_pka();

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
    g.setX(0,L,100,false);
    g.setY(0,L,1,true);
    g.setZ(0,L,1,true);

    box3D box = g.box();
    S.fill(box,Fe);

    options opt;
    S.getOptions(opt);
    opt.printJSON(cout);

    return 0;
}

void xs_tst()
{
    XSquad< screeningZBL > xs1;
    XSquad2< screeningZBL > xs2;
    XS_zbl_magic xs3;

    {
        std::ofstream ostrm1("xs1.dat");
        std::ofstream ostrm2("xs2.dat");
        std::ofstream ostrm3("xs3.dat");
        float e = 1.e-5f;
        for(int i=0; i<7; i++) {
            cout << "E = " << e << endl;
            float s = 1.e-9f;
            //cout << "  S = " << s << endl;
            ostrm1 << xs1.sin2Thetaby2(e,s);
            ostrm2 << xs2.sin2Thetaby2(e,s);
            ostrm3 << xs3.sin2Thetaby2(e,s);
            for(int j=1; j<=1000; j++) {
                s = 50.f*j/1000;
                //cout << "  S = " << s << endl;
                ostrm1 << '\t' << xs1.sin2Thetaby2(e,s);
                ostrm2 << '\t' << xs2.sin2Thetaby2(e,s);
                ostrm3 << '\t' << xs3.sin2Thetaby2(e,s);
            }
            ostrm1 << std::endl;
            ostrm2 << std::endl;
            ostrm3 << std::endl;
            e *= 10;
        }
    }
}

int gencorteo4bit2()
{
    XSquad2< screeningZBL > xs_quad_zbl;


    cout << "Computing 4-bit corteo scattering table";
    cout.flush();
    // compute matrix for each reduced energy, reduced impact parameter pair

    {
        std::ofstream ofs("corteo4bit.dat");
        for(corteo4bit::e_index ie; ie!=ie.end(); ie++) {
            corteo4bit::s_index is;
            double sin2ThetaBy2 = xs_quad_zbl.sin2Thetaby2(*ie, *is++);
            ofs << std::setprecision(9);
            ofs << sin2ThetaBy2;
            for(; is!=is.end(); is++) {
                // calculations (summation) made using double to decrease numerical noise
                sin2ThetaBy2 = xs_quad_zbl.sin2Thetaby2(*ie, *is);
                ofs << '\t' << sin2ThetaBy2;
            }
            ofs << endl;
        }
    }
    {
        std::ofstream ofs("e.dat");
         for(corteo4bit::e_index i; i!=i.end(); i++) ofs << *i << endl;
    }
    {
        std::ofstream ofs("s.dat");
        for(corteo4bit::s_index i; i!=i.end(); i++) ofs << *i << endl;
    }

    cout << " done.\n";
    return 0;
}