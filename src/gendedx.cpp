/*
 * Adapted from Corteo by F. Schiettekatte
 * http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
 *
 * gendedx
 *
 * Utility to generate dE/dx files for all elements in all elements for Z=1(H) .. 92(U)
 * using the SRIM utility SRModule
 *
 * Usage:
 *   gendedx srmodule_folder srmodule_cmd output_folder
 *
 * where
 *   srmodule_folder : the path to the folder of SRModule.exe
 *   srmodule_cmd : the command used to run SRModule. In Linux this could be "wine SRModule.exe"
 *   output_folder : a folder to store the output files
 *
 * Example:
 *   ./gendedx "/home/george/.local/SRIM-2013/srmodule" \
 *             'WINEPREFIX=/home/george/.local/share/wineprefixes/Wine32 WINEARCH="win32" wine SRModule.exe' \
 *             dedx
 *
 * The contents of the output folder should be copied to iradinapp/dedx/ in order to build the iondedx library
 *
 *
 */

#include "dedx.h"
#include "elements.h"

#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

#ifdef WIN32
const char endline[] = "\n";
#else
// in Linux or mac, where SRModule is called using Wine, we must feed an input file containing return+linefeed
const char endline[] = "\r\n";
#endif

namespace fs = std::filesystem;
using namespace std;


string srmodule_folder;
string srmodule_command;
string out_folder;

int parse(int argc, char* argv[]);
int doit();

int main(int argc, char* argv[])
{
    if (parse(argc,argv)<0) return -1;
    return doit();
}

int parse(int argc, char* argv[])
{
    const char* usage = "Usage: gendedx srmodule_folder srmodule_cmd out_folder";
    if (argc!=4) {
        cout << usage << endl;
        return -1;
    }
    srmodule_folder = string(argv[1]);
    srmodule_command = string(argv[2]);
    out_folder = string(argv[3]);
    return 0;
}

int run_srim(int Z1, int Z2, vector<float>& data)
{
    // store current folder and change to srim srmodule folder
    fs::path work_dir(fs::current_path());
    fs::current_path(srmodule_folder);

    // create SR.IN file
    {
        // write SRModule input file
        ofstream os("SR.IN");
        os << endline << endline
           << "SR_OUTPUT.txt"
           << endline << endline
           << Z1 << ' ' << elements::mostAbundantIsotope(Z1)
           << endline << endline
           << " 0 1 0"
           << endline << endline
           << "1"
           << endline << endline
           << Z2 << " \"X\" 1 1"
           << endline << endline
           << " 7"
           << endline << endline
           << " 0 0 " << endline;

        for(const float erg : dedx_index()) os << erg/1000 << endline;

        os << " 0 ";
    }

    // call srmodule
    int ret = system(srmodule_command.c_str());

    // if ok read output
    if (ret==0)
    {
        ifstream is("SR_OUTPUT.txt");
        char buff[4096];
        is.getline(buff,4096);
        is.getline(buff,4096);
        is.getline(buff,4096);
        is.getline(buff,4096);

        int i=0;
        while(is.good() && i<data.size())
        {
            double d, dedxval;
            is >> d >> data[i];
            is.getline(buff,4096);
            i++;
        }

        if (!is.good()) {
            cerr << "Failed to read SR_OUTPUT.txt" << endl;
            ret = -1;
        }

    }

    // return to workdir
    fs::current_path(work_dir);
    return ret;
}

// print out a 32bit float so that it is read back the same
std::ostream& printfloat(std::ostream& os, float x)
{
    char buff[32];
    // get number of digits for a float -> text -> float round-trip
    static constexpr auto d = std::numeric_limits<float>::max_digits10;
    std::sprintf(buff,"%.*g", d, x);
    os << buff;
    return os;
}

int print_cpp_matrix(ofstream& Z1cpp, int Z1, int Z2, const vector<float>& data)
{
    Z1cpp << "static const float "
          << elements::name(Z1)
          << "_on_"
          << elements::name(Z2)
          << "[] = {";
    const int val_per_line = 10;
    for(int i=0; i<data.size()-1; ++i) {
        if (i % val_per_line == 0) Z1cpp << endl << "    ";
        printfloat(Z1cpp, data[i]);
        Z1cpp << ", ";
    }
    Z1cpp << data.back() << "};" << endl << endl;
    return 0;
}

int print_cpp_element_matrix(ofstream& Z1cpp, int Z1, int Z2max)
{
    Z1cpp << "const float* dedx"
          << elements::name(Z1)
          << "[] = { nullptr,";
    const int val_per_line = 5;
    for(int Z2=1; Z2<=Z2max; Z2++) {
        if (((Z2-1) % val_per_line) == 0) Z1cpp << endl << "    ";
        Z1cpp << elements::name(Z1)
              << "_on_"
              << elements::name(Z2);
        if (Z2<Z2max) Z1cpp << ", ";
        else Z1cpp << "};" << endl << endl;
    }
    return 0;
}

int doit()
{
    const char* preample = "// file generated by gendedx - do not edit\n\n";

    vector<float> data(dedx_index::size);

    int Zmax = elements::max_atomic_num;
    cout << "Generating dedx data:" << endl;
    for(int Z1=1; Z1<=Zmax; Z1++) {
        cout << elements::name(Z1) << ' ';
        cout.flush();
        string fname(out_folder);
        fname += "/";
        fname += elements::name(Z1);
        fname += ".cpp";
        ofstream Z1cpp(fname);
        // Z1cpp << fixed << setprecision(4);
        Z1cpp << preample;
        for(int Z2=1; Z2<=Zmax; Z2++) {
            if (Z2 % 10 == 0) {
                cout << '.'; cout.flush();
            }
            if (run_srim(Z1,Z2,data)<0) return -1;
            print_cpp_matrix(Z1cpp,Z1,Z2,data);
        }
        print_cpp_element_matrix(Z1cpp, Z1, Zmax);
        cout << endl;
    }

    { // generate master file dedx_data.cpp
        string fname(out_folder);
        fname += "/dedx_data.cpp";
        ofstream cpp(fname);
        cpp << preample;

        for(int Z1=1; Z1<=Zmax; Z1++) {
            cpp << "extern const float* dedx"
                << elements::name(Z1) << "[];" << endl;
        }
        cpp << endl;

        cpp << "const float** dedx_data[] = { nullptr, ";
        const int val_per_line = 5;
        for(int Z1=1; Z1<=Zmax; Z1++) {
            if (((Z1-1) % val_per_line) == 0) cpp << endl << "    ";
            cpp << "dedx" << elements::name(Z1);
            if (Z1<Zmax) cpp << ", ";
            else cpp << "};" << endl << endl;
        }
     }

    return 0;
}



