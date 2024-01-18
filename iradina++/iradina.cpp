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

int test1(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    return test1(argc,argv);
}

struct runner {
    simulation_base* s;
    int run() { return s->run(); }
};

int test1(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.ini]" << endl;
        return 1;
    }

    std::string path = argv[1];
    simulation_base* S;

    int ret = parse(path, S);
    if (ret!=0) return ret;

    int nthreads = S->getParameters().threads;
    unsigned int N = S->getParameters().max_no_ions;
    unsigned int Nth = N/nthreads;

    cout << endl << endl << "Starting simulation ..." << endl << endl;

    // TIMING
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    S->setMaxIons(Nth);
    ret = S->init();
    if (ret!=0) return ret;

    std::vector< runner > sims(nthreads);
    std::seed_seq sseq{1, 2, 3, 4, 5};
    std::vector<std::uint32_t> seeds(nthreads);
    sseq.generate(seeds.begin(), seeds.end());
    sims[0].s = S;
    S->seed(seeds[0]);
    for(int i=1; i<nthreads; i++) {
        parse(path, sims[i].s);
        sims[i].s->setMaxIons(Nth);
        sims[i].s->seed(seeds[i]);
        sims[i].s->init();
    }

    std::vector< std::thread* > threads;
    for(int i=0; i<nthreads; i++) {
        threads.push_back(new std::thread(&runner::run, &sims[i]));
    }

    // waiting for threads to finish...
    for(int i=0; i<nthreads; i++) threads[i]->join();

    // consolidate results
    for(int i=1; i<nthreads; i++)
        S->addTally(sims[i].s->getTally());

    tally t = S->getTally();

    // CALC TIME/ion CLOCK_PROCESS_CPUTIME_ID
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double ms_per_ion_ = (end.tv_sec - start.tv_sec) * 1.e3 / t.Nions() / nthreads;
    ms_per_ion_ += 1.e-6*(end.tv_nsec - start.tv_nsec) / t.Nions() / nthreads;



    S->saveTallys();



    cout << endl << endl << "Completed." << endl;
    cout << "ms/ion = " << ms_per_ion_ << endl;
    cout << "ion/s = " << floor(1000/ms_per_ion_) << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Nrecoils()/t.Npkas() << endl;

    // delete threads
    for(int i=0; i<nthreads; i++) {
        delete threads[i];
        delete sims[i].s;
    }

    return ret ? 0 : -1;
}


int test0(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.ini]" << endl;
        return 1;
    }

    std::string path = argv[1];
    simulation_base* S;

    int ret = parse(path, S);
    if (ret!=0) return ret;

    cout << endl << endl << "Starting simulation ..." << endl << endl;

    ret = (S->init()==0) && (S->run()==0);

    S->saveTallys();

    tally t = S->getTally();

    cout << endl << endl << "Completed." << endl;
    cout << "ms/ion = " << S->ms_per_ion() << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Nrecoils()/t.Npkas() << endl;

    return ret ? 0 : -1;
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

    material* Fe = S->addMaterial("Fe", 7.8658);
    Fe->addAtom(Zfe, Mfe, 1.f, 40.f,
                  1.f, 1.f, 40.f);

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
    cout << "ms/ion = " << S->ms_per_ion() << endl;
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
