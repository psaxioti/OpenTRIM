#ifndef EVENT_STREAM_H
#define EVENT_STREAM_H

#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

class event_stream;
class ion;

class event
{
protected:
    std::vector<float> buff_;
public:
    event()
    {}
    virtual ~event()
    {}
    void resize(int sz) { buff_.resize(sz); }
    void reset() {
        if (!buff_.empty())
            std::memset(buff_.data(),0,buff_.size()*sizeof(float));
    }
    int size() const { return buff_.size(); }
    const float* data() const { return buff_.data(); }
};

class event_stream
{
protected:
    uint rows_, cols_;
    std::ofstream fs_;
    std::string fname_;

public:
    event_stream() :
        rows_(0), cols_(0)
    {}
    virtual ~event_stream()
    {
        close();
    }
    int open(const char* fname, int cols);
    uint rows() const { return rows_; }
    uint cols() const { return cols_; }
    const std::string& fileName() const { return fname_; }
     void close();
    void write(const event* ev);
    static int merge(const std::vector<event_stream*>& ev,
                     const char* fname, const char* ds_name);
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

public:

    pka_event() :
        event(),
        natoms_(0)
    {}

    static int event_size(int natoms)
    { return 5 + natoms*3; }
    void setNatoms(int n) {
        natoms_ = n;
        resize(event_size(n));
    }
    int ionid() const { return buff_[ofIonId]; }
    int atomid() const { return buff_[ofAtomId]; }
    int cellid() const { return buff_[ofCellId]; }
    float recoilE() const { return buff_[ofErg]; }
    void init(const ion* i);
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

class exit_event : public event
{

    typedef enum {
        ofIonId = 0,
        ofAtomId = 1,
        ofCellId = 2,
        ofErg = 3,
        ofPos = 4,
        ofDir = 7,
        ofS = 10,
        ofEnd = 11
    } offset_t;

    // IonId AtomId CellId erg pos(3) dir(3) s

public:

    exit_event() : event()
    {
        resize(ofEnd);
    }

    void set(const ion* i, int cellid, float s);
};

#endif // EVENT_STREAM_H
