#ifndef OUT_FILE_H
#define OUT_FILE_H

namespace H5 {
class H5File;
};

class simulation_base;

class out_file
{

    const simulation_base* sim_;
    H5::H5File * h5f;

public:
    out_file(const simulation_base* s);
    ~out_file();

    int open(const char* fname);

    int save();

    void close();
};

#endif // OUT_FILE_H
