#include "simulation.h"
#include "elements.h"
#include "options.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;

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
    if(argc != 2)
    {
        cerr << "usage: iradina++ [ConfigFile.json]" << endl;
        return -1;
    }

    cout << "Parsing config file ... " << endl;
    std::ifstream is(argv[1]);
    options s;
    if (s.parseJSON(is)!=0) return -1;

    cout << endl << endl;
    cout << "Starting simulation ..." << endl << endl;

    simulation* S = s.createSimulation();
    S->init();
    S->exec(my_progress_callback);
    S->saveTallys();

    tally t = S->getTally();
    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << S->ips() << " (total), ";
    cout << S->ips()/S->nThreads() << " (per thread)" << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Ndisp()/t.Npkas() << endl;

    return 0;
}





