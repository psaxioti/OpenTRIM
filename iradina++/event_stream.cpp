#include "event_stream.h"

#include <H5Cpp.h>

using namespace H5;

int event_stream::open(const char* fname, int cols)
{
    close();
    fname_ = fname;
    cols_ = cols;
    rows_ = 0;
    fs_.open(fname_, std::ios::binary);
    return fs_.is_open() ? 0 : -1;
}

void event_stream::close()
{
    if (fs_.is_open()) fs_.close();
}
void event_stream::write(const event* ev)
{
    if (fs_.is_open()) {
        fs_.write((char*)ev->data(),ev->size()*sizeof(float));
        rows_++;
    }
}

int event_stream::merge(const std::vector<event_stream *> ev, const char* fname, const char* ds_name)
{
    uint cols = ev[0]->cols();
    uint rows = ev[0]->rows();
    for(int i = 1; i<ev.size(); i++) {
        assert(ev[i]->cols()==cols);
        rows += ev[i]->rows();
    }
    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
        //Exception::dontPrint();
        // Open the file and dataset.
        H5File f(fname, H5F_ACC_TRUNC);

        hsize_t dims[2];
        dims[0] = rows;
        dims[1] = cols;
        DataSpace filespace(2,dims);
        // Create the dataset.
        DataSet dataset = f.createDataSet(ds_name, PredType::NATIVE_FLOAT, filespace);

        // mem buffer ~1MB
        dims[0] = (1 << 10);
        std::vector<float> buff(dims[0]*cols);

        hsize_t offset[2] = {0, 0};
        hsize_t dims1[2];
        dims1[1] = cols;

        for(const event_stream* e : ev) {

            std::ifstream ifs(e->fileName(), std::ios::binary);

            // copy data in chunks
            hsize_t nrows = e->rows(); // total rows to copy
            while (nrows) {
                dims1[0] = std::min(nrows,dims[0]); // # of rows to copy in this iter
                nrows -= dims1[0];

                // read from raw buffer
                ifs.read((char*)buff.data(), dims1[0]*cols*sizeof(float));

                // write to file
                DataSpace memspace(2, dims1);
                // Select a hyperslab in dataset.
                filespace.selectHyperslab(H5S_SELECT_SET, dims1, offset);
                dataset.write(buff.data(), PredType::NATIVE_FLOAT, memspace, filespace);
                offset[0] += dims1[0];

            }

            ifs.close();

            std::remove(e->fileName().c_str());
        }
    }  // end of try block

    // catch failure caused by the H5File operations
    catch (FileIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException error) {
        error.printErrorStack();
        return -1;
    }

    return 0;  // successfully terminated
}
