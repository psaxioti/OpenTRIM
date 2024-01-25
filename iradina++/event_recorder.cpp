#include "event_recorder.h"

#include <H5Cpp.h>

using namespace H5;

int event_recorder::init(event* ev, uint buff_rows)
{
    close();
    fname_ = "";
    if (ev_) delete ev_;
    ev_ = ev;
    cols_ = ev->size();
    rows_ = buff_rows;
    buff_.resize(rows_*cols_,0.f);
    file_row_cnt_ = 0;
    buff_idx_ = 0;
    return 0;
}

int event_recorder::open(const char* fname)
{
    fname_ = fname;
    fname_ += ".bin";
    file_row_cnt_ = 0;
    ofs.open(fname_,std::ios::binary);
    return ofs.is_open() ? 0 : -1;
}

void event_recorder::flush()
{
    if (!ofs.is_open()) return;
    file_row_cnt_ += buff_idx_;
    ofs.write((const char*)buff_.data(), buff_idx_*cols_*sizeof(float));
}

void event_recorder::close()
{
    flush();
    if (ofs.is_open()) ofs.close();
}

int event_recorder::merge(const std::vector<event_recorder *> ev, const char* fname, const char* ds_name)
{
    uint cols = ev[0]->cols();
    uint rows = ev[0]->file_rows();
    for(int i = 1; i<ev.size(); i++) {
        assert(ev[i]->cols()==cols);
        rows += ev[i]->file_rows();
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

        for(const event_recorder* e : ev) {

            std::ifstream ifs(e->name(), std::ios::binary);

            // copy data in chunks
            hsize_t nrows = e->file_rows(); // total rows to copy
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

            std::remove(e->name().c_str());
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
