#ifndef OUT_FILE_H
#define OUT_FILE_H

namespace H5 {
class H5File;
};

class mccore;
class options;

class out_file
{

    const mccore* sim_;
    H5::H5File * h5f;

public:
    out_file(const mccore* s);
    ~out_file();

    int open(const char* fname);

    int save(const options& opt);

    void close();
};

#endif // OUT_FILE_H
