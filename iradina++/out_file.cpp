#include "out_file.h"
#include "mccore.h"
#include "dedx.h"
#include "mcdriver.h"

#include <H5Cpp.h>

using namespace H5;

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

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {

        /*
     * Turn off the auto-printing when failure occurs so that we can
     * handle the errors appropriately
     */
        Exception::dontPrint();

        /*
     * Create dataspace for the dataset in the file.
     */
        DataSpace fspace; // default = scalar

        /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
        DataSet dataset = f->createDataSet( name, h5traits<T>::predType(), fspace );

        /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
        dataset.write( &data, h5traits<T>::predType() );

    }  // end of try block

    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        //error.printError();
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
        Exception::dontPrint();
        hsize_t dims = 1;
        DataSpace fspace(1,&dims); // default = scalar
        StrType type(PredType::C_S1, len);
        DataSet dataset = f->createDataSet( name, type, fspace );
        dataset.write( data.c_str(), type );
    }
    catch( FileIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        //error.printError();
        return -1;
    }
    return 0;
}

template<typename T>
int save_array(H5File* f, const char* name,
               const T* data,
               const std::vector<hsize_t>& dims)
{

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {

    /*
     * Turn off the auto-printing when failure occurs so that we can
     * handle the errors appropriately
     */
    Exception::dontPrint();

    /*
     * Create dataspace for the dataset in the file.
     */
    DataSpace fspace( dims.size(), dims.data() );

    /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
    DataSet dataset = f->createDataSet( name, h5traits<T>::predType(), fspace );

    /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
    dataset.write( data, h5traits<T>::predType() );

    }  // end of try block

    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        //error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        //error.printError();
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

template<typename T, class _ArrT>
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

out_file::out_file(const mccore *s) :
    sim_(s), h5f(nullptr)
{

}

out_file::~out_file()
{
    close();
    if (h5f) delete h5f;
}

int out_file::open(const char* fname)
{


    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
        /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
        Exception::dontPrint();
        /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
       h5f = new H5File(fname, H5F_ACC_TRUNC);
    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        return -1;
    }

    return 0;  // successfully terminated

}

int out_file::save(const options &opt)
{
    // save options
    std::stringstream ss;
    opt.printJSON(ss);
    save_string(h5f, "json_options", ss.str());


    const tally& t = sim_->getTally();
    save_scalar(h5f, "Nions", t.Nions());
    save_scalar(h5f, "Npkas", t.Npkas());
    save_scalar(h5f, "Ndisp", t.Ndisp());
    save_scalar(h5f, "Nrepl", t.Nrepl());
    save_scalar(h5f, "Nimpl", t.Nimpl());
    save_scalar(h5f, "Nvac", t.Nvac());
    save_scalar(h5f, "Nlost", t.Nlost());

    // save grid
    auto grid = sim_->getTarget().grid();
    save_array<float, grid1D>(h5f, "X", grid.x());
    save_array<float, grid1D>(h5f, "Y", grid.y());
    save_array<float, grid1D>(h5f, "Z", grid.z());
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

        save_array_nd(h5f, "cell_xyz", buff);
    }


    // save tallys
    if (opt.Simulation.simulation_type == mccore::FullCascade) {

        bool ret = true;
        int k = 0;
        while(ret && k<tally::std_tallies) {
            ret = ret && save_array_nd(h5f, tally::arrayName(k), t.at(k))==0;
            k++;
        }

        if (!ret) return -1;

    } else {

        bool ret = true;

        if (!ret) return -1;
    }

    if (opt.Output.store_dedx) {
        save_array_nd(h5f,"dEdx",sim_->dedx());
        save_array_nd(h5f,"dEstrag",sim_->de_strag());
        ArrayND<float> dedx_erg(dedx_index::size);
        for(dedx_index i; i<i.end(); i++) dedx_erg(i) = *i;
        save_array_nd(h5f,"dEdx_erg",dedx_erg);
        save_array_nd(h5f,"maxImpactPar",sim_->max_impact_par());
        save_array_nd(h5f,"max_fp",sim_->max_fp());
    }

    return 0;
}

void out_file::close()
{
    if (h5f) h5f->close();
}
