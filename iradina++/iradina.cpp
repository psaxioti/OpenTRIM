#include "mcdriver.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char* argv[])
{
    if(argc == 2)
    {
        static const char* usage = "iradina++ v0.1.0\n"
                                   "Usage: iradina++ [options]\n"
                                   "Reads json config from stdin.\n"
                                   "Use iradina++ < [file.json] to read from a file.\n"
                                   "Options:\n"
                                   "  -h: display usage\n";

        cerr << usage << endl;
        return -1;
    }

    cout << "Parsing json config ... " << endl;
    options s;
    if (s.parseJSON(cin)!=0) return -1;

    cout << endl << endl;
    cout << "Starting simulation ..." << endl << endl;

    mcdriver D;
    D.setOptions(s);
    D.exec();
    D.saveTallys();

    tally t = D.getSim()->getTally();
    cout << endl << endl << "Completed." << endl;
    cout << "ion/s = " << D.ips() << " (total), ";
    cout << D.ips()/D.nThreads() << " (per thread)" << endl;
    cout << "Ions = " << t.Nions() << endl;
    cout << "PKA/Ion = " << 1.f*t.Npkas()/t.Nions() << endl;
    cout << "Recoils/PKA = " << 1.f*t.Ndisp()/t.Npkas() << endl;

    return 0;
}





