#ifndef ARRAYS_H
#define ARRAYS_H

#include <memory>
#include <vector>

template<typename Scalar>
class Array1D
{
    typedef Array1D<Scalar> _Self;

    struct Private
    {
        std::vector<Scalar> buffer;

        Private() : buffer()
        {}

        Private(int r) : buffer(r,Scalar(0))
        {}

        Private(const Private& p) : buffer(p.buffer)
        {}

        void resize(int r) {
            buffer.resize(r);
        }
    };

    std::shared_ptr<Private> P_;

    Array1D(const Private& p) :
        P_(new Private(p))
    {}

public:

    Array1D() :
        P_()
    {}
    explicit Array1D(int rows) :
        P_(new Private(rows))
    {}
    Array1D(const Array1D& other) : P_(other.P_)
    {}
    Array1D& operator=(const Array1D& other)
    { P_ = other.P_; return *this; }

    Array1D copy() const
    {
        return P_.get() ? Array1D(*P_) : Array1D();
    }

    bool isNull() const { return P_ == nullptr; }

    int rows() const { return P_->buffer.size(); }
    int size() const { return P_->buffer.size(); }

    Scalar& operator()(int i) { return P_->p[i]; }
    const Scalar& operator()(int i) const { return P_->p[i]; }

    Scalar& operator[](int i) { return P_->p[i]; }
    const Scalar& operator[](int i) const { return P_->p[i]; }
    Scalar* data() { return P_->buffer.data(); }
    const Scalar* data() const { return P_->buffer.data(); }

    Array1D& operator+=(const Array1D& a) {
        if (size()==a.size()) {
            Scalar* p = data();
            Scalar* pend = p + size();
            const Scalar* q = a.data();
            while(p<pend) *p++ += *q++;
        }
        return *this;
    }

};

template<typename Scalar>
class Array2D
{
    typedef Array2D<Scalar> _Self;

    struct Private
    {
        std::vector<Scalar> buffer;
        std::vector<Scalar*> p;
        int rows, cols;

        Private() :
            rows(0), cols(0), buffer(), p()
        {}
        Private(int r, int c) :
            rows(r), cols(c),
            buffer(r*c,Scalar(0)), p(r)
        {
            init();
        }
        Private(const Private& p) :
            rows(p.rows), cols(p.cols),
            buffer(p.buffer), p(p.p)
        {}

        void resize(int r, int c) {
            rows = r; cols = c;
            buffer.resize(r*c);
            p.resize(r);
            init();
        }
        void init() {
            for(int i=0; i<rows; i++) p[i] = buffer.data() + i*cols;
        }
    };

    std::shared_ptr<Private> P_;

    Array2D(const Private& p) :
        P_(new Private(p))
    {}

public:
    Array2D() :
        P_()
    {}
    Array2D(int rows, int cols) :
        P_(new Private(rows,cols))
    {}
    Array2D(const Array2D& other) : P_(other.P_)
    {}
    Array2D& operator=(const Array2D& other)
    { P_ = other.P_; return *this; }

    Array2D copy() const
    {
        return P_.get() ? Array2D(*P_) : Array2D();
    }

    bool isNull() const { return P_ == nullptr; }

    int rows() const { return P_->rows; }
    int cols() const { return P_->cols; }
    int size() const { return P_->buffer.size(); }

    Scalar& operator()(int i, int j) { return P_->p[i][j]; }
    const Scalar& operator()(int i, int j) const { return P_->p[i][j]; }

    Scalar* operator[](int i) { return P_->p[i]; }
    const Scalar* operator[](int i) const { return P_->p[i]; }
    Scalar* data() { return P_->buffer.data(); }
    const Scalar* data() const { return P_->buffer.data(); }

    Array2D& operator+=(const Array2D& a) {
        if (size()==a.size()) {
            Scalar* p = data();
            Scalar* pend = p + size();
            const Scalar* q = a.data();
            while(p<pend) *p++ += *q++;
        }
        return *this;
    }

};

template<typename Scalar>
class Array3D
{
    typedef Array3D<Scalar> _Self;

    struct Private
    {
        std::vector<Scalar> buffer;
        std::vector<Scalar*> p1;
        std::vector<Scalar**> p;
        int rows, cols, layers;
        Private(int r, int c, int l) :
            rows(r), cols(c), layers(l),
            buffer(r*c*l,Scalar(0)), p(r), p1(r*c)
        {
            init();
        }
        Private() :
            rows(0), cols(0), layers(0),
            buffer(), p(), p1()
        {}
        Private(const Private& p) :
            rows(p.rows), cols(p.cols), layers(p.layers),
            buffer(p.buffer), p(p.p), p1(p.p1)
        {}
        void resize(int r, int c, int l) {
            rows = r; cols = c; layers = l;
            buffer.resize(r*c*l);
            p1.resize(r*c);
            p.resize(r);
            init();
        }
        void init() {
            for(int i=0; i<rows; i++) {
                p[i] = p1.data() + i*cols;
                for(int j=0; j<cols; j++) {
                    int k = i*cols+j;
                    p1[k] = buffer.data() + k*layers;
                }
            }
        }
    };

    std::shared_ptr<Private> P_;

    Array3D(const Private& p) :
        P_(new Private(p))
    {}

public:
    Array3D() :
        P_()
    {}
    Array3D(int rows, int cols, int layers) :
        P_(new Private(rows,cols,layers))
    {}
    Array3D(const Array3D& other) : P_(other.P_)
    {}
    Array3D& operator=(const Array3D& other)
    { P_ = other.P_; return *this; }
    Array3D copy() const
    {
        return P_ ? Array3D(*P_) : Array3D();
    }

    bool isNull() const { return !P_; }

    int rows() const { return P_->rows; }
    int cols() const { return P_->cols; }
    int layers() const { return P_->layers; }
    int size() const { return P_->buffer.size(); }

    Scalar& operator()(int i, int j, int k) { return P_->p[i][j][k]; }
    const Scalar& operator()(int i, int j, int k) const { return P_->p[i][j][k]; }

    Scalar** operator[](int i) { return P_->p[i]; }
    const Scalar** operator[](int i) const { return (const Scalar**)(P_->p[i]); }
    Scalar* data() { return P_->buffer.data(); }
    const Scalar* data() const { return P_->buffer.data(); }

    Array3D& operator+=(const Array3D& a) {
        if (size()==a.size()) {
            Scalar* p = data();
            Scalar* pend = p + size();
            const Scalar* q = a.data();
            while(p<pend) *p++ += *q++;
        }
        return *this;
    }

};

typedef Array1D<float> Array1Df;
typedef Array2D<float> Array2Df;
typedef Array3D<float> Array3Df;

typedef Array1D<double> Array1Dd;
typedef Array2D<double> Array2Dd;
typedef Array3D<double> Array3Dd;

typedef Array1D<unsigned int> Array1Dui;
typedef Array2D<unsigned int> Array2Dui;
typedef Array3D<unsigned int> Array3Dui;

#endif // ARRAYS_H
