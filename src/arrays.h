#ifndef ARRAYS_H
#define ARRAYS_H

#include <memory>
#include <cstring>
#include <cassert>
#include <vector>

template<typename Scalar>
class ArrayND
{
    typedef ArrayND<Scalar> _Self;

    struct Private
    {
        std::vector<size_t> dim;
        std::vector<Scalar> buffer;

        static size_t calc_size(const std::vector<size_t>& d)
        {
            size_t S(1);
            for(auto s : d) S *= s;
            return S;
        }

        Private() : dim(), buffer()
        {}

        explicit Private(size_t i) :
            dim(1),
            buffer(i,Scalar(0))
        {
            dim[0] = i;
        }

        Private(size_t i, size_t j) :
            dim(2),
            buffer(i*j,Scalar(0))
        {
            dim[0] = i; dim[1] = j;
        }

        Private(size_t i, size_t j, size_t k) :
            dim(3),
            buffer(i*j*k,Scalar(0))
        {
            dim[0] = i; dim[1] = j; dim[2] = k;
        }

        Private(const std::vector<size_t>& d) :
            dim(d), buffer(calc_size(d),Scalar(0))
        {}

        Private(const Private& p) : dim(p.dim), buffer(p.buffer)
        {}

        void resize(const std::vector<size_t>& d) {
            dim = d;
            buffer.resize(calc_size(d));
        }

        size_t idx(int i, int j) {
            return (dim.size() == 2) ? i*dim[1] + j : 0;
        }
        size_t idx(int i, int j, int k) {
            return (dim.size() == 3) ? (i*dim[1] + j)*dim[2] + k : 0;
        }
        size_t idx(const std::vector<size_t>& d) {
            assert(dim.size() == d.size());
            size_t k = d[0];
            for(int i = 1; i<d.size(); i++)  k = k*dim[i] + d[i];
            return k;
        }
    };

    std::shared_ptr<Private> P_;

    ArrayND(const Private& p) :
        P_(new Private(p))
    {}

public:

    ArrayND() = default;
    explicit ArrayND(const std::vector<size_t>& dim) :
        P_(new Private(dim))
    {}
    explicit ArrayND(int i) :
        P_(new Private(i))
    {}
    explicit ArrayND(int i, int j) :
        P_(new Private(i,j))
    {}
    explicit ArrayND(int i, int j, int k) :
        P_(new Private(i,j,k))
    {}
    ArrayND(const ArrayND& other) : P_(other.P_)
    {}
    ArrayND& operator=(const ArrayND& other)
    { P_ = other.P_; return *this; }

    ArrayND copy() const
    {
        return P_.get() ? ArrayND(*P_) : ArrayND();
    }

    bool isNull() const { return P_ == nullptr; }

    const std::vector<size_t>& dim() const { return P_->dim; }
    size_t size() const { return P_->buffer.size(); }
    int ndim() const { return P_->dim.size(); }

    Scalar& operator[](int i) { return P_->buffer[i]; }
    const Scalar& operator[](int i) const { return P_->buffer[i]; }
    Scalar& operator()(int i) { return P_->buffer[i]; }
    const Scalar& operator()(int i) const { return P_->buffer[i]; }
    Scalar& operator()(int i, int j) { return P_->buffer[P_->idx(i,j)]; }
    const Scalar& operator()(int i, int j) const { return P_->buffer[P_->idx(i,j)]; }
    Scalar& operator()(int i, int j, int k) { return P_->buffer[P_->idx(i,j,k)]; }
    const Scalar& operator()(int i, int j, int k) const { return P_->buffer[P_->idx(i,j,k)]; }

    Scalar* data() { return P_->buffer.data(); }
    const Scalar* data() const { return P_->buffer.data(); }

    ArrayND& operator+=(const ArrayND& a) {
        if (size()==a.size()) {
            Scalar* p = data();
            Scalar* pend = p + size();
            const Scalar* q = a.data();
            while(p<pend) *p++ += *q++;
        }
        return *this;
    }

    void addSquared(const ArrayND& a) {
        if (size()==a.size()) {
            Scalar* p = data();
            Scalar* pend = p + size();
            const Scalar* q = a.data();
            while(p<pend) *p++ += (*q)*(*q++);
        }
    }

    void clear() {
        Scalar* p = data();
        std::memset(p,0,size()*sizeof(Scalar));
    }

};

typedef ArrayND<double> ArrayNDd;
typedef ArrayND<float>  ArrayNDf;

#endif // ARRAYS_H
