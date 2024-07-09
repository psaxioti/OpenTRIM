#ifndef EVENT_STREAM_H
#define EVENT_STREAM_H

#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

typedef unsigned int uint;

class event_stream;
class ion;

class event
{
protected:
    std::vector<float> buff_;
    std::vector<std::string> columnNames_;
    std::vector<std::string> columnDescriptions_;
public:
    event()
    {}
    event(size_t n,
          const std::vector<std::string> &names,
          const std::vector<std::string> &descs) : 
          buff_(n),
          columnNames_(names),
          columnDescriptions_(descs)
    {}
    event(const event& ev) :
    buff_(ev.buff_), 
    columnNames_(ev.columnNames_),
    columnDescriptions_(ev.columnDescriptions_)
    {}
    virtual ~event()
    {}
    void resize(size_t sz) { buff_.resize(sz); }
    void reset() {
        if (!buff_.empty())
            std::memset(buff_.data(),0,buff_.size()*sizeof(float));
    }
    size_t size() const { return buff_.size(); }
    const float* data() const { return buff_.data(); }
    const std::vector<std::string>& columnNames() const
    { return columnNames_;}
    const std::vector<std::string>& columnDescriptions() const
    { return columnDescriptions_;}
};

class event_stream
{
protected:
    size_t rows_, cols_;
    std::ofstream fs_;
    std::string fname_;
    event event_proto_;

public:
    event_stream() :
        rows_(0), cols_(0)
    {}
    virtual ~event_stream()
    {
        close();
        remove();
    }
    int open(const std::string& fname, const event& ev);
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    const std::string& fileName() const { return fname_; }
    void close();
    void remove();
    void write(const event* ev);
    int merge(const event_stream& ev);
    bool is_open() const { return fs_.is_open(); }
    const event& event_prototype() const { return event_proto_; }
};

class pka_event_recorder;

class pka_event : public event
{
    int natoms_;

    enum offset_t {
        ofIonId = 0,
        ofAtomId = 1,
        ofCellId = 2,
        ofErg = 3,
        ofTdam = 4,
        ofVac = 5
    };

public:

    pka_event() :
        event(),
        natoms_(0)
    {}

    static int event_size(int natoms)
    { return 5 + natoms*3; }
    void setNatoms(int n, const std::vector<std::string>& labels);
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

    enum offset_t {
        ofIonId = 0,
        ofAtomId = 1,
        ofCellId = 2,
        ofErg = 3,
        ofPos = 4,
        ofDir = 7,
        ofEnd = 10
    };

    // IonId AtomId CellId erg pos(3) dir(3) s

public:

    exit_event();
    void set(const ion* i, int cellid);
};

#endif // EVENT_STREAM_H
