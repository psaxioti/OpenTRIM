#ifndef OUT_FILE_H
#define OUT_FILE_H

namespace H5 {
class H5File;
};

class simulation;
class options;

class out_file
{

    const simulation* sim_;
    H5::H5File * h5f;

public:
    out_file(const simulation* s);
    ~out_file();

    int open(const char* fname);

    int save(const options& opt);

    void close();
};

#endif // OUT_FILE_H
