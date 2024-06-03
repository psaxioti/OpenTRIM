#ifndef CORTEO_INTERP_H
#define CORTEO_INTERP_H

#include "corteo.h"

#include <cmath>
#include <array>

template<typename T, class idx_t>
T corteo_lin_interp(const T* y, const T& x) {

    static_assert(sizeof(T) == sizeof(typename idx_t::RealType),
                  "real types must have the same size");

    if (x <= idx_t::minVal) return y[0];
    if (x >= idx_t::maxVal) return y[idx_t::size-1];
    idx_t i(x), j(i+1);
    T t = *i;
    t = (x-t)/(*j-t);
    return y[i]*(1-t) + y[j]*t;
}

template<typename T, class idx_t>
T corteo_log_interp(const T* y, const T& x) {

    static_assert(sizeof(T) == sizeof(typename idx_t::RealType),
                  "real types must have the same size");

    if (x <= idx_t::minVal) return y[0];
    if (x >= idx_t::maxVal) return y[idx_t::size-1];
    idx_t i(x), j(i+1);
    T t = std::log2(*i);
    t = (std::log2(x)-t)/(std::log2(*j)-t);
    return std::exp2(std::log2(y[i])*(1-t) + std::log2(y[j])*t);
}

template<class idx_t>
class corteo_lin_interp1D {

public:
    typedef typename idx_t::RealType RealType;

    corteo_lin_interp1D()
    {}

    template<class Cont>
    explicit corteo_lin_interp1D(const Cont& y)
    { set(y); }

    template<class Cont>
    void set(const Cont& y) {
        for(idx_t i,j(1); i<i.end()-1; i++,j++) {
            y_[i] = y[i];
            dydx_[i] = (y[j] - y[i]) / (*j - *i);
        }
        y_[idx_t::size-1] = y[idx_t::size-1];
    }

    RealType operator()(const RealType& x) const
    {
        if (x <= idx_t::minVal) return y_.front();
        if (x >= idx_t::maxVal) return y_.back();
        idx_t i(x);
        return y_[i] + dydx_[i]*(x-*i);
    }

    const RealType* data() const { return &y_[0]; }

private:
    std::array<RealType, idx_t::size> y_, dydx_;

};

template<class idx_t>
class corteo_log_interp1D {

public:
    typedef typename idx_t::RealType RealType;

    corteo_log_interp1D()
    {
        for(idx_t i; i<i.end(); i++) logx_[i] = std::log2(*i);
    }

    template<class Cont>
    explicit corteo_log_interp1D(const Cont& y)
    {
        for(idx_t i; i<i.end(); i++) logx_[i] = std::log2(*i);
        set(y);
    }

    template<class Cont>
    void set(const Cont& y) {
        for(idx_t i,j(1); i<i.end()-1; i++,j++) {
            y_[i] = y[i];
            logy_[i] = std::log2(y[i]);
            dydx_[i] = (std::log2(y[j]) - std::log2(y[i])) /
                       (logx_[j] - logx_[i]);
        }
        logy_[idx_t::size-1] = std::log2(y[idx_t::size-1]);
        y_[idx_t::size-1] = y[idx_t::size-1];
    }

    RealType operator()(const RealType& x) const
    {
        if (x <= idx_t::minVal) return y_.front();
        if (x >= idx_t::maxVal) return y_.back();
        idx_t i(x);
        return std::exp2(logy_[i] + dydx_[i]*(std::log2(x)-logx_[i]));
    }

    const RealType* data() const { return &y_[0]; }

private:
    std::array<RealType, idx_t::size> y_, logx_, logy_, dydx_;
};

#endif // CORTEO_INTERP_H
