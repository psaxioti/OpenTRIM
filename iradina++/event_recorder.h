#ifndef EVENT_RECORDER_H
#define EVENT_RECORDER_H

#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

class event_recorder;

class event
{
protected:
    std::string name_;
    float* buff_;
    int size_;

    void setBuff(float* q) { buff_ = q; }

    friend class event_recorder;

public:
    explicit event(const char* name, int sz) :
        name_(name),
        buff_(nullptr),
        size_(sz)
    {}
    virtual ~event()
    {}
    int size() const { return size_; }
    const std::string& name() { return name_; }
};

class event_recorder
{
public:
    typedef unsigned long long ullong;

protected:
    std::vector<float> buff_;
    uint rows_, cols_;
    uint file_row_cnt_;
    uint buff_idx_;
    std::ofstream ofs;
    std::string fname_;
    event* ev_;

public:
    event_recorder() :
        rows_(0), cols_(0),
        file_row_cnt_(0), buff_idx_(0),
        ev_(nullptr)
    {}
    virtual ~event_recorder()
    {
        close();
        if (ev_) delete ev_;
    }
    int init(event* ev, uint buff_rows);
    int open(const char* fname);
    uint buff_rows() const { return rows_; }
    uint file_rows() const { return file_row_cnt_; }
    uint cols() const { return cols_; }
    event& get()
    {
        if (buff_idx_ == rows_) {
            flush();
            // zero-out buffer
            std::memset(buff_.data(), 0, buff_.size()*sizeof(float));
            ev_->setBuff(buff_.data());
            buff_idx_ = 1;
        } else {
            ev_->setBuff(buff_.data() + buff_idx_*cols_);
            buff_idx_++;
        }
        return *ev_;
    }
    void flush();
    const std::string& name() const { return fname_; }
    void close();
    static int merge(const std::vector<event_recorder*> ev, const char* fname, const char* ds_name);
};

class pka_event_recorder;

class pka_event : public event
{
    int natoms_;

    typedef enum {
        ofIonId = 0,
        ofAtomId = 1,
        ofCellId = 2,
        ofErg = 3,
        ofTdam = 4,
        ofVac = 5
    } offset_t;

    friend class pka_event_recorder;

public:

    explicit pka_event(int natoms) :
        event("pka",event_size(natoms)),
        natoms_(natoms)
    {}

    static int event_size(int natoms)
    { return 5 + natoms*3; }
    void init(uint ion_id, uint atom_id, uint cell_id, float erg)
    {
        buff_[ofIonId] = ion_id;
        buff_[ofAtomId] = atom_id;
        buff_[ofCellId] = cell_id;
        buff_[ofErg] = erg;
        addVac(atom_id-1);
    }
    float& Tdam() { return buff_[ofTdam]; }
    const float& Tdam() const { return buff_[ofTdam]; }
    void addVac(int atom_id)
    {
        assert(atom_id >= 0);
        buff_[ofVac + atom_id] += 1.f;
    }
    const float& Vac(int atom_id) const { return buff_[ofVac + atom_id]; }
    void addRepl(int atom_id) {
        assert(atom_id >= 0);
        buff_[ofVac + natoms_ + atom_id] += 1.f;
    }
    const float& Repl(int atom_id) const { return buff_[ofVac + natoms_ + atom_id]; }
    void addImpl(int atom_id) {
        assert(atom_id >= 0);
        buff_[ofVac + 2*natoms_ + atom_id] += 1.f;
    }
    const float& Impl(int atom_id) const { return buff_[ofVac + 2*natoms_ + atom_id]; }
};

class pka_event_recorder : public event_recorder
{
public:
    pka_event_recorder() :
        event_recorder()
    {}
    pka_event& get() {
       event_recorder::get();
       return *(reinterpret_cast<pka_event*>(ev_));
    }
};

#endif // EVENT_RECORDER_H
