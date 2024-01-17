#include "out_file.h"
#include "simulation.h"

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

template <typename T>
struct array_traits<T, Array1D<T> > {
    static std::vector<hsize_t> dims(const Array1D<T>& A) {
        std::vector<hsize_t> d(1);
        d[0] = A.rows();
        return d;
    }
    typedef h5traits<T> scalar_traits;
};
template <typename T>
struct array_traits<T, Array2D<T> > {
    static std::vector<hsize_t> dims(const Array2D<T>& A) {
        std::vector<hsize_t> d(2);
        d[0] = A.rows();
        d[1] = A.cols();
        return d;
    }
    typedef h5traits<T> scalar_traits;
};
template <typename T>
struct array_traits<T, Array3D<T> > {
    static std::vector<hsize_t> dims(const Array3D<T>& A) {
        std::vector<hsize_t> d(3);
        d[0] = A.rows();
        d[1] = A.cols();
        d[2] = A.layers();
        return d;
    }
    typedef h5traits<T> scalar_traits;
};
template <>
struct array_traits<float, grid1D > {
    static std::vector<hsize_t> dims(const grid1D& A) {
        std::vector<hsize_t> d(1);
        d[0] = A.size();
        return d;
    }
    typedef h5traits<float> scalar_traits;
};

template<typename T, class _ArrT>
int save_array(H5File* f, const char* name, const _ArrT& A) {
    std::vector<hsize_t> dims = array_traits<T, _ArrT>::dims(A);
    return save_array(f, name, A.data(), dims);
}

out_file::out_file(const simulation_base *s) :
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

int out_file::save()
{
    const tally& t = sim_->getTally();
    save_scalar(h5f, "Nions", t.Nions());
    save_scalar(h5f, "Npkas", t.Npkas());
    save_scalar(h5f, "Nrecoils", t.Nrecoils());
    save_scalar(h5f, "Nrepl", t.Nrepl());
    save_scalar(h5f, "Nimpl", t.Nimpl());
    save_scalar(h5f, "Nvac", t.Nvac());

    // save grid
    save_array<float, grid1D>(h5f, "X", sim_->grid().x());
    save_array<float, grid1D>(h5f, "Y", sim_->grid().y());
    save_array<float, grid1D>(h5f, "Z", sim_->grid().z());


    // save tallys
    if (sim_->simulationType() == simulation_base::FullCascade) {

        bool ret =
            save_array<unsigned int, Array2Dui>(h5f, "Implantations", t.implantations())==0 &&
            save_array<unsigned int, Array2Dui>(h5f, "Replacements", t.replacements())==0 &&
            save_array<unsigned int, Array2Dui>(h5f, "Vacancies", t.vacancies())==0 &&
            save_array<double, Array2Dd>(h5f, "IonizationEnergy", t.ionization())==0 &&
            save_array<double, Array2Dd>(h5f, "PhononEnergy", t.phonons())==0;

        if (!ret) return -1;

    } else {

        bool ret = true;

        if (!ret) return -1;
    }
    return 0;
}

void out_file::close()
{
    if (h5f) h5f->close();
}
