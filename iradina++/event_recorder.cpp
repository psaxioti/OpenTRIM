#include "event_recorder.h"

#include <H5Cpp.h>

using namespace H5;

int event_recorder::open(const char* fname)
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

        /*
         * Create dataspace for the dataset in the file.
         */
        file_idx_ = 0;
        hsize_t dims[2];
        dims[0] = 0;
        dims[1] = cols_;
        hsize_t maxdims[2];
        maxdims[0] = H5S_UNLIMITED; // unlimited rows
        maxdims[1] = cols_;

        //h5sp = new  DataSpace(2, dims, maxdims);
        DataSpace dspace(2, dims, maxdims);

        // Modify dataset creation property to enable chunking
        DSetCreatPropList plist;
        dims[0] = rows_;
        plist.setChunk(2,dims);

        // Create the chunked dataset.  Note the use of pointer.
        h5ds = new DataSet(
            h5f->createDataSet(ev_->name(), PredType::NATIVE_FLOAT, dspace, plist)
            );

    }  // end of try block

    // catch failure caused by the H5File operations
    catch (FileIException error) {
        error.printErrorStack();
        close();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException error) {
        error.printErrorStack();
        close();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException error) {
        error.printErrorStack();
        close();
        return -1;
    }

    return 0;  // successfully terminated
}

void event_recorder::flush()
{
    if (h5f) {

        int nwrite = buff_idx_;

        // this points to the end of the data
        hsize_t offset[2];
        offset[0] = file_idx_;
        offset[1] = 0;

        // extend dataset
        file_idx_ += nwrite;
        hsize_t dims[2];
        dims[0] = file_idx_;
        dims[1] = cols_;
        h5ds->extend(dims);

        // Select a hyperslab in extended portion of the dataset.
        DataSpace filespace(h5ds->getSpace());
        dims[0] = nwrite;
        dims[1] = cols_;
        filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);

        // Define memory space.
        DataSpace memspace(2, dims);

        // Write data to the extended portion of the dataset.
        h5ds->write(buff_.data(), PredType::NATIVE_FLOAT,
                    memspace, filespace);
    }
}

void event_recorder::close()
{
    if (h5ds) {
        delete h5ds;
        h5ds = nullptr;
    }
    if (h5f) {
        delete h5f;
        h5f = nullptr;
    }
}

int event_recorder::merge(const char* fname1, const char* fname2, const char* ds_name)
{
    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
        Exception::dontPrint();

        hsize_t dims1[2], dims2[2], offset1[2], offset2[2], chunk_dims2[2];

        // Open the file and dataset.
        H5File f1(fname1, H5F_ACC_RDWR);
        H5File f2(fname2, H5F_ACC_RDONLY);

        DataSet ds1(f1.openDataSet(ds_name));
        DataSet ds2(f2.openDataSet(ds_name));

        // Get the dataset's dataspace and creation property list.
        {
            DataSpace fspace1(ds1.getSpace());
            DataSpace fspace2(ds2.getSpace());
            DSetCreatPropList prop = ds2.getCreatePlist();

            // Get information to obtain memory dataspace.
            fspace1.getSimpleExtentDims(dims1);
            fspace2.getSimpleExtentDims(dims2);


            if (dims1[1]!=dims2[1]) return -1; // wrong # of cols
            if (dims2[0]==0) return 0; // empty 2nd file

            if (H5D_CHUNKED == prop.getLayout()) {
                int rank_chunk = prop.getChunk(2, chunk_dims2);
                assert(rank_chunk == 2);
            } else {
                chunk_dims2[0] = 32;
                chunk_dims2[1] = dims1[1];
            }
        }

        // offset1 = end of data, offset2 = start of data
        offset1[0] = dims1[0]; offset2[0] = 0;
        offset1[1] = offset2[1] = 0;

        // extend dataset in file 1
        dims1[0] += dims2[0];
        ds1.extend(dims1);

        // copy data in chunks
        hsize_t nrows = dims2[0]; // total rows to copy
        hsize_t crows = std::max(hsize_t(32), chunk_dims2[0]); // rows in chunk
        std::vector<float> buff_(dims2[1]*crows); // memory buffer
        DataSpace fspace1(ds1.getSpace());
        DataSpace fspace2(ds2.getSpace());
        while (nrows) {
            size_t n = std::min(nrows,crows); // # of rows to copy in this iter
            nrows -=n;

            dims1[0] = n;
            DataSpace memspace(2, dims1, NULL);

            // read from file2
            // Select a hyperslab in extended portion of dataset1.
            fspace2.selectHyperslab(H5S_SELECT_SET, dims1, offset2);
            ds2.read(buff_.data(), PredType::NATIVE_FLOAT, memspace, fspace2);
            offset2[0] += n;

            // write to file1
            // Select a hyperslab in extended portion of dataset1.
            fspace1.selectHyperslab(H5S_SELECT_SET, dims1, offset1);
            ds1.write(buff_.data(), PredType::NATIVE_FLOAT, memspace, fspace1);
            offset1[0] += n;

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
