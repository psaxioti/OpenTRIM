#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <algorithm>

#include <unsupported/Eigen/AlignedVector3>

typedef Eigen::AlignedVector3<float> vector3;

typedef Eigen::AlignedVector3<int> ivector3;

typedef Eigen::AlignedBox1f box1D;
typedef Eigen::AlignedBox2f box2D;
typedef Eigen::AlignedBox3f box3D;
typedef Eigen::AlignedBox1i ibox1D;

struct grid1D : public std::vector<float> {

    typedef std::vector<float> vec_t;

    float w;
    bool equispaced;
    bool periodic;

    void set(const vec_t& x) {
        vec_t::assign(x.begin(),x.end());
        w = back() - front();
        equispaced = false;
    }

    void set(const float& x0, const float& x1, int n) {
        resize(n+1);
        w = x1-x0;
        float step = w/n;
        for(int i=0; i<=n; i++) (*this)[i] = x0 + i*step;
        equispaced = true;
    }

    bool contains(const float& x) const {
        return (x>=front()) && (x<back());
    }

    bool periodic_contains(const float& x) const {
        return periodic ? 1 : (x>=front()) && (x<back());
    }

    bool check_pos(float& x) const {
        if (x<front()) {
            if (periodic) {
                do x += w; while (x<front());
                return true;
            } else return false;
        } else if (x>=back()) {
            if (periodic) {
                do x -= w; while (x>=back());
                return true;
            } else return false;
        } else return true;
    }

    bool cell_contains(int i, const float& x) const {
        return (x>=at(i)) && (x<at(i+1));
    }

    // should be called only if contains==true !!
    int find_index(const float& x) const {
        if (size()==2) return 0;
        if (equispaced)
            return std::floor((x - front())/w*(size()-1));
        else
            return std::upper_bound(begin(), end(), x) - begin() - 1;
    }

    void applyBC(float& x) const
    {
        if (periodic) {
            if (x < front())
                do x += w; while( x < front() );
            else if (x >= back())
                do x -= w; while( x >= back() );
        }
    }

    ibox1D range(const box1D& b) const
    {
        int i1(0), i2(size()-2);
        float x1 = 0.5*(at(i1)+at(i1+1));
        float x2 = 0.5*(at(i2)+at(i2+1));
        while (i1<i2 && !b.contains(box1D::VectorType(x1))) {
            i1++; x1 = 0.5*(at(i1)+at(i1+1));
        }
        while (i2>i1 && !b.contains(box1D::VectorType(x2))) {
            i2--; x2 = 0.5*(at(i2)+at(i2+1));
        }
        return ibox1D(i1,i2);
    }

};



class grid3D
{
    grid1D x_, y_, z_;
    box3D box_;

    void calcBox()
    {
        box_.min() = box3D::VectorType(x_.front(),y_.front(),z_.front());
        box_.max() = box3D::VectorType(x_.back(),y_.back(),z_.back());
    }

public:


    typedef enum {
        X=0,
        Y=1,
        Z=2
    } Axis;


    grid3D()
    {
        x_.set(0.f,1.f,1);
        y_.set(0.f,1.f,1);
        z_.set(0.f,1.f,1);
        calcBox();
    }

    void setX(const float& x0, const float& x1, int n, bool periodic = false)
    { x_.set(x0,x1,n); calcBox(); x_.periodic = periodic; }
    void setY(const float& x0, const float& x1, int n, bool periodic = false)
    { y_.set(x0,x1,n); calcBox(); y_.periodic = periodic; }
    void setZ(const float& x0, const float& x1, int n, bool periodic = false)
    { z_.set(x0,x1,n); calcBox(); z_.periodic = periodic; }

    const grid1D& x() const { return x_; }
    const grid1D& y() const { return y_; }
    const grid1D& z() const { return z_; }
    const box3D& box() const { return box_; }

    bool contains(const vector3& v) const {
        return x_.contains(v.x()) &&
               y_.contains(v.y()) &&
               z_.contains(v.z());
    }

    bool periodic_contains(const vector3& v) const {
        return x_.periodic_contains(v.x()) &&
               y_.periodic_contains(v.y()) &&
               z_.periodic_contains(v.z());
    }

    bool check_pos(vector3& v) const {
        return x_.check_pos(v.x()) &&
               y_.check_pos(v.y()) &&
               z_.check_pos(v.z());
    }

    bool contains(const ivector3& i, const vector3& v) const {
        return x_.cell_contains(i.x(),v.x()) &&
            y_.cell_contains(i.y(),v.y()) &&
            z_.cell_contains(i.z(),v.z());
    }

    // should be called only if contains()==true
    ivector3 pos2cell(const vector3& v) const
    {
        return ivector3(x_.find_index(v.x()),
                        y_.find_index(v.y()),
                        z_.find_index(v.z()));
    }

    int cellid(const ivector3& i) const {
        return (i.x()*(y_.size()-1) + i.y())*(z_.size()-1) + i.z();
    }

    static bool isNull(const ivector3& i)  {
        return (i.x()<0) || (i.y()<0) || (i.z()<0);
    }

    static ivector3 nullcell() {
        return ivector3(-1,-1,-1);
    }

    int ncells() const {
        return (x_.size()-1)*
               (y_.size()-1)*
               (z_.size()-1);
    }

    float volume() const {
        return x_.w * y_.w * z_.w;
    }

    void applyBC(vector3& v) const {
        x_.applyBC(v.x());
        y_.applyBC(v.y());
        z_.applyBC(v.z());
    }
};


#endif // GEOMETRY_H
