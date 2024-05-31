#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <algorithm>

#include <unsupported/Eigen/AlignedVector3>

/**
 * \defgroup Geometry Geometry
 *
 * \brief 3D vectors and spatial grid.
 *
 * @{
 *
 * The geometric features of the simulation are defined on a 3D rectangular grid
 * represented by the class grid3D.
 *
 * Space is divided into rectangular cells.
 *
 * The grid3D class offers functions to check if a particular
 * 3D point is within the grid, to which cell it belongs, etc
 *
 * A grid3D comprises of three objects of class grid1D representing
 * the grid in each of the axes.
 *
 * @ingroup MC
 *
 * @}
 */


/**
 * @brief The class vector3 represents a 3D vector.
 *
 * It is based on Eigen::AlignedVector3 with float as the basic type,
 * so that the code can
 * benefit from SSE cpu instructions where possible. For example,
 * vector addition can be performed in one cpu intstruction.
 *
 * @ingroup Geometry
 */
typedef Eigen::AlignedVector3<float> vector3;

/**
 * @brief The class ivector3 represents a 3D index vector.
 *
 * It is based on Eigen::AlignedVector3 with integer as the basic type.
 * Thus, the code can
 * benefit from SSE cpu instructions where possible. For example,
 * vector addition can be performed in one cpu intstruction.
 *
 * Index vectors are used to specify a cell in the 3D rectangular grid.
 *
 * @ingroup Geometry
 */
typedef Eigen::AlignedVector3<int> ivector3;

/**
 * @brief A 1D range based on Eigen::AlignedBox1f
 * @ingroup Geometry
 */
typedef Eigen::AlignedBox1f box1D;

/**
 * @brief A 1D index range based on Eigen::AlignedBox1i
 * @ingroup Geometry
 */
typedef Eigen::AlignedBox1i ibox1D;

/**
 * @brief A rectangular 3D range based on Eigen::AlignedBox3f
 * @ingroup Geometry
 */
typedef Eigen::AlignedBox3f box3D;

/**
 * @brief The grid1D class represents a 1D spatial partition
 *
 * It is essentially an array of \f$ N \f$ monotonically increasing
 * points, \f$ x_0,x_1,...,x_{N-1} \f$, which divides the region
 * between the 1st and the Nth point in \f$ N-1 \f$ cells.
 *
 * The cells do not have to be of equal width, however, currently only
 * equidistant grids have been employed.
 *
 * @ingroup Geometry
 *
 */
class grid1D : public std::vector<float> {

    typedef std::vector<float> vec_t;

    float w_;
    bool equispaced_;
    bool periodic_;

public:
    /// Default constructor creates an empty grid
    grid1D() : vec_t(), w_(0),
        equispaced_(false), periodic_(false)
    {}
    /// Copy constructor
    grid1D(const grid1D& g) : vec_t(g),
        w_(g.w_), equispaced_(g.equispaced_),
        periodic_(g.periodic_)
    {}

    /// Returns the total width \f$ x_{N-1}-x_0 \f$
    float w() const { return w_; }
    /// Returns true if all cells have the same width
    bool equispaced() const { return equispaced_; }
    /// Returns true if the grid has periodic boundary conditions
    bool periodic() const { return periodic_; }
    /// Set periodic boundary conditions on or off depebing on the value of b
    void setPeriodic(bool b) { periodic_ = b; }

    /// Set the grid to the values of the array x
    void set(const vec_t& x) {
        vec_t::assign(x.begin(),x.end());
        w_ = back() - front();
        equispaced_ = false;
    }
    /// Create an equidistant grid from x0 to x1 divided in n cells
    void set(const float& x0, const float& x1, int n) {
        resize(n+1);
        w_ = x1-x0;
        w_ = w_/n;
        for(int i=0; i<=n; i++) (*this)[i] = x0 + i*w_;
        equispaced_ = true;
    }

    /// Returns true if x is within the grid region
    bool contains(const float& x) const {
        return (x>=front()) && (x<back());
    }

    /**
     * @brief Returns true if x is within the grid region, anticipating periodic boundary conditions
     *
     * If periodic boundary conditions are on and x is outside the grid, then x will be
     * brought within the range by adding or subtracting an integer multiple of the
     * grid period (width). The function will return
     * true.
     *
     * @param x is the posistion to check
     * @return true if x is within the grid
     */
    bool contains_with_bc(float& x) const {
        if (x<front()) {
            if (periodic_) {
                do x += w_; while (x<front());
                assert(x<back());
                return true;
            } else return false;
        } else if (x>=back()) {
            if (periodic_) {
                do x -= w_; while (x>=back());
                assert(x>=front());
                return true;
            } else return false;
        } else return true;
    }

    /// Returns true if x is inside the i-th cell, \f$ x_i \leq x < x_{i+1} \f$
    bool contains(int i, const float& x) const {
        return (x>=at(i)) && (x<at(i+1));
    }

    /// Returns the cell index i for which \f$ x_i \leq x < x_{i+1} \f$. Should be called only if grid1D::contains(x) returns true
    int pos2cell(const float& x) const {
        assert(contains(x));
        if (size()==2) return 0;
        if (equispaced_) {
            return std::floor((x - front())/w_);
        }
        else
            return std::upper_bound(begin(), end(), x) - begin() - 1;
    }

    /**
     * @brief Return the range of cell indexes that are within the spatial range b
     *
     * The function returns the range of cells that have their centers
     * within the interval defined by b.
     *
     *
     * @param b the range as a box1D
     * @return the range of cells as an ibox1D
     */
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


/**
 * @brief The grid3D class represents a 3D rectangular grid.
 *
 * It comprises of three grid1D objects, one for each dimension,
 * which can be accessed by the funtions grid3D::x(), grid3D::y() and
 * grid3D::z().
 *
 * @ingroup Geometry
 */
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

    /// Default constructor creates empty grid
    grid3D()
    {
        x_.set(0.f,1.f,1);
        y_.set(0.f,1.f,1);
        z_.set(0.f,1.f,1);
        calcBox();
    }
    /// Copy constructor
    grid3D(const grid3D& g) :
        x_(g.x_), y_(g.y_), z_(g.z_), box_(g.box_)
    {}

    /// Set the x-axis grid to an equidistant partition of n cells
    void setX(const float& x0, const float& x1, int n, bool periodic)
    { x_.set(x0,x1,n); calcBox(); x_.setPeriodic(periodic); }
    /// Set the y-axis grid to an equidistant partition of n cells
    void setY(const float& x0, const float& x1, int n, bool periodic)
    { y_.set(x0,x1,n); calcBox(); y_.setPeriodic(periodic); }
    /// Set the z-axis grid to an equidistant partition of n cells
    void setZ(const float& x0, const float& x1, int n, bool periodic)
    { z_.set(x0,x1,n); calcBox(); z_.setPeriodic(periodic); }

    /// Return the x-axis grid
    const grid1D& x() const { return x_; }
    /// Return the y-axis grid
    const grid1D& y() const { return y_; }
    /// Return the z-axis grid
    const grid1D& z() const { return z_; }
    /// Return the rectangular box containing the whole 3D grid
    const box3D& box() const { return box_; }

    /// Return the i-th cell rectangular box
    box3D box(const ivector3& i) const {
        vector3 X0(x_[i.x()],y_[i.y()],z_[i.z()]);
        vector3 X1(x_[i.x()+1],y_[i.y()+1],z_[i.z()+1]);
        return box3D(X0, X1);
    }

    /// Return true if v is within the grid
    bool contains(const vector3& v) const {
        return x_.contains(v.x()) &&
               y_.contains(v.y()) &&
               z_.contains(v.z());
    }
    /// Return true if v is within the grid, adjusting v for periodic boundaries
    bool contains_with_bc(vector3& v) const {
        return x_.contains_with_bc(v.x()) &&
               y_.contains_with_bc(v.y()) &&
               z_.contains_with_bc(v.z());
    }
    /// Return true if v is within the cell with indexes i
    bool contains(const ivector3& i, const vector3& v) const {
        return x_.contains(i.x(),v.x()) &&
            y_.contains(i.y(),v.y()) &&
            z_.contains(i.z(),v.z());
    }
    /// Return the cell index for v; call only if contains()==true
    ivector3 pos2cell(const vector3& v) const
    {
        assert(contains(v));
        return ivector3(x_.pos2cell(v.x()),
                        y_.pos2cell(v.y()),
                        z_.pos2cell(v.z()));
    }

    /// Return the id of the i-th cell
    int cellid(const ivector3& i) const {
        return (i.x()*(y_.size()-1) + i.y())*(z_.size()-1) + i.z();
    }

    static bool isNull(const ivector3& i)  {
        return (i.x()<0) || (i.y()<0) || (i.z()<0);
    }

    static ivector3 nullcell() {
        return ivector3(-1,-1,-1);
    }
    /// Returns the total number of cells
    int ncells() const {
        return (x_.size()-1)*
               (y_.size()-1)*
               (z_.size()-1);
    }
    /// Returns the total volume
    float volume() const {
        return x_.w() * y_.w() * z_.w();
    }
};


#endif // GEOMETRY_H
