#include "mccore.h"
#include "dedx.h"
#include "mcdriver.h"

#include <iomanip>
#include <iostream>
#include <H5Cpp.h>

using namespace H5;

using std::cerr;
using std::endl;

template <typename T> struct h5traits;

template<>
struct h5traits<double> {
    static const PredType& predType() { return PredType::NATIVE_DOUBLE; }
};
template<>
struct h5traits<float> {
    static const PredType& predType() { return PredType::NATIVE_FLOAT; }
};
template<>
struct h5traits<unsigned int> {
    static const PredType& predType() { return PredType::NATIVE_UINT; }
};
template<>
struct h5traits<int> {
    static const PredType& predType() { return PredType::NATIVE_INT; }
};

template<typename T>
int save_scalar(H5File* f, const char* name, const T& data)
{
    try
    {   
        DataSpace fspace; // default = scalar
        LinkCreatPropList lcpl;
        lcpl.setCreateIntermediateGroup(true);
        DataSet dataset = f->createDataSet(name, h5traits<T>::predType(), fspace,
                                           DSetCreatPropList::DEFAULT, DSetAccPropList::DEFAULT,
                                           lcpl );
        dataset.write( &data, h5traits<T>::predType() );

    }  // end of try block

    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    return 0;
}

int save_string(H5File* f, const char* name, const std::string& data)
{
    int len = data.size();
    if (len==0) return 0;

    try
    {
        //Exception::dontPrint();
        hsize_t dims = 1;
        DataSpace fspace(1,&dims); // default = scalar
        StrType type(PredType::C_S1, len);
        LinkCreatPropList lcpl;
        lcpl.setCreateIntermediateGroup(true);
        DataSet dataset = f->createDataSet(name, type, fspace,
                                           DSetCreatPropList::DEFAULT, DSetAccPropList::DEFAULT,
                                           lcpl);
        dataset.write( data.c_str(), type );
    }
    catch( FileIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    return 0;
}

int save_string(H5File* f, const char* name, const std::vector<std::string>& strvec)
{
    int n = strvec.size();
    int stride = strvec[0].size();
    for (int i=0; i<n; ++i)
        if (strvec[i].size()> stride) stride = strvec[i].size();
    stride++; // +1 for the null term
    std::vector<char> buff(n*stride,'\0');
    for (int i=0; i<n; ++i)
        memcpy(buff.data() + i*stride,
               strvec[i].data(),
               strvec[i].size());

    try
    {
        //Exception::dontPrint();
        hsize_t dims = n;
        DataSpace fspace(1,&dims); // default = scalar
        StrType type(PredType::C_S1, stride);
        //type.setStrpad(H5T_STR_NULLPAD);
        LinkCreatPropList lcpl;
        lcpl.setCreateIntermediateGroup(true);
        DataSet dataset = f->createDataSet(name, type, fspace,
                                           DSetCreatPropList::DEFAULT, DSetAccPropList::DEFAULT,
                                           lcpl);

        dataset.write( buff.data(), type );
    }
    catch( FileIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    return 0;
}

template<typename T>
int save_array(H5File* f, const char* name,
               const T* data,
               const std::vector<hsize_t>& dims)
{
    try
    {
        DataSpace fspace( dims.size(), dims.data() );
        LinkCreatPropList lcpl;
        lcpl.setCreateIntermediateGroup(true);
        DataSet dataset = f->createDataSet(name, h5traits<T>::predType(), fspace,
                                           DSetCreatPropList::DEFAULT, DSetAccPropList::DEFAULT,
                                           lcpl);

        dataset.write( data, h5traits<T>::predType() );

    }  // end of try block

    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
        cerr << error.getDetailMsg() << endl;
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( GroupIException error )
    {
        cerr << error.getDetailMsg() << endl;
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        cerr << error.getDetailMsg() << endl;
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        cerr << error.getDetailMsg() << endl;
        error.printErrorStack();
        return -1;
    }
    return 0;
}

template <typename T, class _A> struct array_traits;


template <>
struct array_traits<float, grid1D > {
    static std::vector<hsize_t> dims(const grid1D& A) {
        std::vector<hsize_t> d(1);
        d[0] = A.size();
        return d;
    }
    typedef h5traits<float> scalar_traits;
};

template <>
struct array_traits<float, ArrayND<float> > {
    static std::vector<hsize_t> dims(const ArrayND<float>& A) {
        std::vector<hsize_t> d(A.ndim());
        for(int i=0; i<d.size(); i++) d[i] = A.dim()[i];
        return d;
    }
    typedef h5traits<float> scalar_traits;
};

template <>
struct array_traits<double, ArrayND<double> > {
    static std::vector<hsize_t> dims(const ArrayND<double>& A) {
        std::vector<hsize_t> d(A.ndim());
        for(int i=0; i<d.size(); i++) d[i] = A.dim()[i];
        return d;
    }
    typedef h5traits<double> scalar_traits;
};

template <typename T>
struct array_traits<T, std::vector<T> > {
    static std::vector<hsize_t> dims(const std::vector<T>& A) {
        std::vector<hsize_t> d(1);
        d[0] = A.size();
        return d;
    }
    typedef h5traits<T> scalar_traits;
};

template<typename T, class _ArrT = std::vector<T> >
int save_array(H5File* f, const char* name, const _ArrT& A) {
    std::vector<hsize_t> dims = array_traits<T, _ArrT>::dims(A);
    return save_array(f, name, A.data(), dims);
}

template<typename T>
int save_array_nd(H5File* f, const char* name, const ArrayND<T>& A) {
    std::vector<hsize_t> dims(A.ndim());
    for(int i=0; i<dims.size(); i++) dims[i] = A.dim()[i];
    return save_array(f, name, A.data(), dims);
}

template<typename T>
int save_array_normalized(H5File* f, const char* name, const ArrayND<T>& A, const uint& N) {
    std::vector<T> a(A.size());
    for(int i=0; i<A.size(); i++) a[i] = A[i]/N;
    std::vector<hsize_t> dims(A.ndim());
    for(int i=0; i<dims.size(); i++) dims[i] = A.dim()[i];
    return save_array(f, name, a.data(), dims);
}

template<typename T>
int save_array_normalized(H5File* f, const char* name, const char* dname,
                          const ArrayND<T>& A, const ArrayND<T>& dA, const uint& N) {
    assert(A.size()==dA.size());
    std::vector<T> a(A.size()), da(A.size());
    for(int i=0; i<A.size(); i++) {
        a[i] = A[i]/N;
        // error in the mean
        da[i] = std::sqrt((dA[i]/N-a[i]*a[i])/(N-1));
    }
    std::vector<hsize_t> dims(A.ndim());
    for(int i=0; i<dims.size(); i++) dims[i] = A.dim()[i];
    return save_array(f, name, a.data(), dims) + save_array(f, dname, da.data(), dims);
}

int mcdriver::save()
{
    std::string fname(out_opts_.OutputFileBaseName);
    fname += ".h5";
    H5::H5File* h5f(nullptr);
    try
    {
        h5f = new H5File(fname, H5F_ACC_TRUNC);
    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        cerr << error.getDetailMsg() << endl;
        error.printErrorStack();
        return -1;
    }

    options opt;
    getOptions(opt);

    // variable list
    std::stringstream var_list;

    // title
    save_string(h5f, "Title", opt.Output.title);
    var_list << "Title" << '\t' << "string" <<  endl;

    // save options
    {
        std::stringstream ss;
        opt.printJSON(ss);
        save_string(h5f, "config_json", ss.str());
    }
    var_list << "config_json" << '\t' << "string" << '\t' << "JSON formatted simulation options" << endl;
    var_list << endl;

    // save iradina++ version info
    var_list << "Version Info" << endl;
    save_string(h5f, "/version_info/version", IRADINAPP_VERSION);
    save_string(h5f, "/version_info/compiler", COMPILER_ID);
    save_string(h5f, "/version_info/compiler_version", COMPILER_VERSION);
    save_string(h5f, "/version_info/build_system", SYSTEM_ID);
    save_string(h5f, "/version_info/build_time", BUILD_TIME);
    var_list << "/version_info/version" << '\t' << "string" << '\t' << "iradina++ version" << endl;
    var_list << "/version_info/compiler" << '\t' << "string" << '\t' << "compiler" << endl;
    var_list << "/version_info/compiler_version" << '\t' << "string" << '\t' << "compiler version" << endl;
    var_list << "/version_info/build_system" << '\t' << "string" << '\t' << "build system" << endl;
    var_list << "/version_info/build_time" << '\t' << "string" << '\t' << "build timestamp" << endl;
    var_list << endl;

    // save run statistics
    var_list << "Run statistics" << endl;
    const tally& t = s_->getTally();
    const tally& dt = s_->getTallyVar();
    save_scalar(h5f, "/run_stat/Nh", t.Nions());
    var_list << "/run_stat/Nh" << '\t' << "Scalar" << '\t' << "# of histories" << endl;
    save_scalar(h5f, "/run_stat/ips", ips_);
    var_list << "/run_stat/ips" << '\t' << "Scalar" << '\t' << "ions/s" << endl;
    save_scalar(h5f, "/run_stat/cpu_time", t.Nions()/ips_);
    var_list << "/run_stat/cpu_time" << '\t' << "Scalar" << '\t' << "total cpu time [s]" << endl;
    {
        std::stringstream ss;
        ss << std::put_time(std::localtime(&start_time_), "%c %Z");
        save_string(h5f, "/run_stat/start_time", ss.str());
    }
    var_list << "/run_stat/start_time" << '\t' << "string" << '\t' << "start time/date" << endl;
    {
        std::stringstream ss;
        ss << std::put_time(std::localtime(&end_time_), "%c %Z");
        save_string(h5f, "/run_stat/end_time", ss.str());
    }
    var_list << "/run_stat/end_time" << '\t' << "string" << '\t' << "start time/date" << endl;
    var_list << endl;

    // save grid
    var_list << "Spatial Grid" << endl;
    auto grid = s_->getTarget().grid();
    save_array<float, grid1D>(h5f, "/grid/X", grid.x());
    var_list << "/grid/X" << '\t' << "1x" << grid.x().size() << '\t' << "x-axis grid" << endl;
    save_array<float, grid1D>(h5f, "/grid/Y", grid.y());
    var_list << "/grid/Y" << '\t' << "1x" << grid.y().size() << '\t' << "y-axis grid" << endl;
    save_array<float, grid1D>(h5f, "/grid/Z", grid.z());
    var_list << "/grid/Z" << '\t' << "1x" << grid.z().size() << '\t' << "z-axis grid" << endl;
    { // save xyz of each cell center
        int rows = grid.x().size()-1;
        int cols = grid.y().size()-1;
        int layers = grid.z().size()-1;
        ArrayND<double> buff(3,grid.ncells());

        for(int i=0; i<rows; i++)
            for(int j=0; j<cols; j++)
                for(int k=0; k<layers; k++)
                {
                    int l = (i*cols+j)*layers+k;
                    buff(0,l) = 0.5f*(grid.x()[i] + grid.x()[i+1]);
                    buff(1,l) = 0.5f*(grid.y()[j] + grid.y()[j+1]);
                    buff(2,l) = 0.5f*(grid.z()[k] + grid.z()[k+1]);
                }

        save_array_nd(h5f, "/grid/cell_xyz", buff);
        var_list << "/grid/cell_xyz" << '\t' << "[3x" << rows << "]\t" << "cell center coordinates" << endl;
    }
    var_list << endl;

    // save atoms & materials
    var_list << "Atom data" << endl;
    auto atoms = s_->getTarget().atoms();
    std::vector<std::string> atom_labels(atoms.size());
    for(int i=0; i<atoms.size(); i++) {
        std::string& s = atom_labels[i];
        const material* m = atoms[i]->mat();
        s = atoms[i]->name();
        if (m) {
            s += " in ";
            s += m->name();
        } else s += " ion";
    }
    save_string(h5f, "/atom/label", atom_labels);
    var_list << "/atom/label" << '\t' << "string array" << "\t" << "label = [Atom] in [Material]" << endl;
    for(int i=0; i<atoms.size(); i++) atom_labels[i] = atoms[i]->name();
    save_string(h5f, "/atom/name", atom_labels);
    var_list << "/atom/name" << '\t' << "string array" << "\t" << "chemical names" << endl;
    {
        std::vector<float> A(atoms.size());
        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->Z();
        save_array<float>(h5f, "/atom/Z", A);
        var_list << "/atom/Z" << '\t' << "array" << "\t" << "atomic numbers" << endl;

        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->M();
        save_array<float>(h5f, "/atom/M", A);
        var_list << "/atom/M" << '\t' << "array" << "\t" << "atomic masses" << endl;

        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->Ed();
        save_array<float>(h5f, "/atom/Ed", A);
        var_list << "/atom/Ed" << '\t' << "array" << "\t" << "displacement energies" << endl;

        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->El();
        save_array<float>(h5f, "/atom/El", A);
        var_list << "/atom/El" << '\t' << "array" << "\t" << "lattice binding energies" << endl;

        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->Es();
        save_array<float>(h5f, "/atom/Es", A);
        var_list << "/atom/Es" << '\t' << "array" << "\t" << "surface binding energies" << endl;

        for(int i=0; i<atoms.size(); i++) A[i] = atoms[i]->Er();
        save_array<float>(h5f, "/atom/Er", A);
        var_list << "/atom/Er" << '\t' << "array" << "\t" << "replacement energies" << endl;
    }
    var_list << endl;



    // save tallys
    var_list << "Tallies" << endl;
    var_list << "  Results are mean values over the ion histories. " << endl
             << "  [VarName]_std is the standard error of the mean." << endl;
    // PKAs
    uint N = t.Nions();
    double vm = t.Npkas()/N, dvm = std::sqrt((dt.Npkas()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Npkas", vm);
    var_list << "/tally/totals/Npkas" << '\t' << "Scalar" << '\t' << "# of PKAs" << endl;
    save_scalar(h5f, "/tally/totals/Npkas_std", dvm);
    var_list << "/tally/totals/Npkas_std" << '\t' << "Scalar" << '\t' << "(std dev) # of PKAs" << endl;

    // Disp
    vm = t.Ndisp()/N; dvm = std::sqrt((dt.Ndisp()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Ndisp", vm);
    var_list << "/tally/totals/Ndisp" << '\t' << "Scalar" << '\t' << "# of displacements" << endl;
    save_scalar(h5f, "/tally/totals/Ndisp_sig", dvm);
    var_list << "/tally/totals/Ndisp_sig" << '\t' << "Scalar" << '\t' << "(std dev) # of displacements" << endl;

    // Repl
    vm = t.Nrepl()/N; dvm = std::sqrt((dt.Nrepl()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Nrepl", vm);
    var_list << "/tally/totals/Nrepl" << '\t' << "Scalar" << '\t' << "# of replacements" << endl;
    save_scalar(h5f, "/tally/totals/Nrepl_sig", dvm);
    var_list << "/tally/totals/Nrepl_sig" << '\t' << "Scalar" << '\t' << "(std dev) # of replacements" << endl;

    // Nimpl
    vm = t.Nimpl()/N; dvm = std::sqrt((dt.Nimpl()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Nimpl", vm);
    var_list << "/tally/totals/Nimpl" << '\t' << "Scalar" << '\t' << "# of implanted/interstitial ions" << endl;
    save_scalar(h5f, "/tally/totals/Nimpl_sig", dvm);
    var_list << "/tally/totals/Nimpl_sig" << '\t' << "Scalar" << '\t' << "(std dev) # of implanted/interstitial ions" << endl;

    // Nvac
    vm = t.Nvac()/N; dvm = std::sqrt((dt.Nvac()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Nvac", vm);
    var_list << "/tally/totals/Nvac" << '\t' << "Scalar" << '\t' << "# of vacancies" << endl;
    save_scalar(h5f, "/tally/totals/Nvac_sig", dvm);
    var_list << "/tally/totals/Nvac_sig" << '\t' << "Scalar" << '\t' << "(std dev) # of vacancies" << endl;

    // Nlost
    vm = t.Nlost()/N; dvm = std::sqrt((dt.Nlost()/N-vm*vm)/(N-1));
    save_scalar(h5f, "/tally/totals/Nlost", vm);
    var_list << "/tally/totals/Nlost" << '\t' << "Scalar" << '\t' << "# of lost ions" << endl;
    save_scalar(h5f, "/tally/totals/Nlost_sig", dvm);
    var_list << "/tally/totals/Nlost_sig" << '\t' << "Scalar" << '\t' << "(std dev) # of lost ions" << endl;

    if (opt.Simulation.simulation_type == mccore::FullCascade) {

        bool ret = true;
        int k = 1;

        while(ret && k<tally::std_tallies) {
            std::string name("/tally/");
            name += tally::arrayGroup(k);
            name += "/";
            name += tally::arrayName(k);
            std::string dname(name);
            dname += "_std";
            ret = ret && save_array_normalized(h5f, name.c_str(), dname.c_str(),
                                               t.at(k), dt.at(k), t.Nions())==0;
            if (ret) {
                auto dim = t.at(k).dim();
                var_list << name << '\t' << '[' << dim[0];
                for(int i=1; i<dim.size(); ++i)
                    var_list << 'x' << dim[i];
                var_list << ']' << '\t' << tally::arrayDescription(k) << endl;

                var_list << dname << '\t' << '[' << dim[0];
                for(int i=1; i<dim.size(); ++i)
                    var_list << 'x' << dim[i];
                var_list << ']' << '\t' << "(std. dev) " << tally::arrayDescription(k) << endl;

            }
            k++;
        }

        if (!ret) return -1;

    } else {

        bool ret = true;

        if (!ret) return -1;
    }

    if (opt.Output.store_dedx) {
        ArrayNDf A(dedx_index::size);
        for(dedx_index i; i<i.end(); i++) A(i) = *i;
        save_array_nd(h5f,"/eels/dEdx_erg",A);

        ArrayND<mccore::dedx_interp_t*> D = s_->dedx();
        A = ArrayNDf(D.dim()[0],D.dim()[1],dedx_index::size);
        for(int i=0; i<A.dim()[0]; i++)
            for(int j=0; j<A.dim()[1]; j++)
                memcpy(&A(i,j,0),D(i,j)->data(),dedx_index::size);

        save_array_nd(h5f,"/eels/dEdx",A);
        save_array_nd(h5f,"/eels/dEstrag",s_->de_strag());
        save_array_nd(h5f,"/eels/mfp",s_->mfp());
        save_array_nd(h5f,"/eels/ipmax",s_->ipmax());
    }


    save_string(h5f, "variable_list", var_list.str());

    h5f->close();

    return 0;
}
