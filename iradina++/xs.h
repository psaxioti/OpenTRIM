#ifndef _XS_H_
#define _XS_H_

#include <cmath>
#include <vector>

#include "arrays.h"
#include "corteo.h"

#define BOHR_RADIUS 0.05291772108 /* Bohr radius [nm] */
#define SCREENCONST 0.88534 /* Lindhard screening length constant */
#define E2 1.43996445      /* e^2 / 4 pi eps0 = e^2 c^2 in [eV][nm] */
#define AMUbyE 1.036426883E-8 /* amu/e = 1.660538782E-27/1.602176487E-19 */
#define ELEMENTARY_CHARGE 1.602176487E-19

class reducedXS
{
public:

    typedef enum {
        GaussMehlerQuad,
        ZBL_Magick,
        Corteo4bitTable
    } IntegrationMethod;

    typedef enum {
        Unscreened = 0,
        Moliere,
        KrC,
        LenzJensen,
        ZBL
    } ScreeningPotentialType;

    struct reducedXS_impl {
        virtual IntegrationMethod integrationMethod() const = 0;
        virtual ScreeningPotentialType potentialType() const = 0;
        virtual double screeningLength(int Z1, int Z2) const = 0;
        virtual double screeningFunction(const double& x) const = 0;
        virtual double sin2Thetaby2(const double& epsilon, const double& s) const = 0;
        virtual double sin2Thetaby2(int ie, int is) const { return 0.; }
        virtual int init(std::ostream* os) { return -1; }
    };

private:

    std::shared_ptr<reducedXS_impl> xs_;

public:
    reducedXS();
    reducedXS(ScreeningPotentialType S, IntegrationMethod I, unsigned int n = 100);
    reducedXS(const reducedXS& other) : xs_(other.xs_)
    {}

    IntegrationMethod integrationMethod() const
    { return xs_->integrationMethod(); }
    ScreeningPotentialType potentialType() const
    { return xs_->potentialType(); }
    double screeningLength(int Z1, int Z2) const
    { return xs_->screeningLength(Z1,Z2); }
    double screeningFunction(const double& x) const
    { return xs_->screeningFunction(x); }
    double sin2Thetaby2(const double& epsilon, const double& s) const
    { return xs_->sin2Thetaby2(epsilon, s); }
    double sin2Thetaby2(int ie, int is) const
    { return xs_->sin2Thetaby2(ie, is); }
    int init(std::ostream* os = NULL) { return xs_->init(os); }
};

template<reducedXS::ScreeningPotentialType _T>
struct screening_base {
    static reducedXS::ScreeningPotentialType type() { return _T; }
    static const char* name() {
        const char* names[] = {
            "Unscreened Coulomb",
            "Moliere",
            "Kr-C",
            "Lenz-Jensen",
            "Ziegler-Biersack-Littmark (ZBL)"
        };
        return names[(int)_T];
    }
};

template<reducedXS::ScreeningPotentialType _T>
struct screening : public screening_base<_T> {
};

template<>
struct screening<reducedXS::Unscreened> : public screening_base<reducedXS::Unscreened>  {
    static double screeningLength(int Z1, int Z2) { return BOHR_RADIUS; } // Bohr radius in nm
    static double screeningFunction(const double& x) { return 1.; }
};
template<>
struct screening<reducedXS::LenzJensen> :
     public screening_base<reducedXS::LenzJensen>
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/std::sqrt(std::pow(Z1,2./3)+std::pow(Z2, 2./3));
    }
    static double screeningFunction(const double& x) {
        double y = 3.108*std::sqrt(x);
        return exp(-y)*(1.+y*(1.+y*(0.3344+y*(0.0485+2.647e-3*y))));
    }
};
template<>
struct screening<reducedXS::KrC> : public screening_base<reducedXS::KrC>
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Screening function coefficients */
    constexpr static const double C[] =
        {0.190945, 0.473674, 0.335381};
    constexpr static const double A[] =
        {0.278544, 0.637174, 1.919249};

    static double screeningFunction(const double& x) {
        return C[0]*exp(-A[0]*x)+C[1]*exp(-A[1]*x)+C[2]*exp(-A[2]*x);
    }
};
template<>
struct screening<reducedXS::Moliere> : public screening_base<reducedXS::Moliere>
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Screening function coefficients */
    constexpr static const double C[] =
        {0.35, 0.55, 0.10};
    constexpr static const double A[] =
        {0.30, 1.20, 6.00};

    static double screeningFunction(const double& x) {
        return C[0]*exp(-A[0]*x)+C[1]*exp(-A[1]*x)+C[2]*exp(-A[2]*x);
    }
};
template<>
struct screening<reducedXS::ZBL> : public screening_base<reducedXS::ZBL>
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Universal screening function coefficients TRIM85 */
    constexpr static const double C[] =
        {0.18175, 0.50986, 0.28022, 0.028171};
    constexpr static const double A[] =
        {3.19980, 0.94229, 0.40290, 0.201620};

    static double screeningFunction(const double& x) {
        return C[0]*exp(-A[0]*x)+C[1]*exp(-A[1]*x)+C[2]*exp(-A[2]*x)+C[3]*exp(-A[3]*x);
    }
};

/*
 * ion scattering cross section using Gauss-Mehler quadrature
 * see Yuan et al. NIMB1993
 */
template<class PHI = screening<reducedXS::ZBL> >
struct reducedXS_quad : public reducedXS::reducedXS_impl
{
    reducedXS_quad(unsigned int nsum = 100) : NSUM(nsum)
    {}

    typedef PHI screening_function;

    virtual reducedXS::IntegrationMethod integrationMethod() const override
    { return reducedXS::GaussMehlerQuad; }
    virtual reducedXS::ScreeningPotentialType potentialType() const override
    { return PHI::type(); }
    virtual double screeningLength(int Z1, int Z2) const override
    {
        return PHI::screeningLength(Z1, Z2);
    }
    virtual double screeningFunction(const double& x) const override
    {
        return PHI::screeningFunction(x);
    }
    virtual double sin2Thetaby2(const double& epsilon, const double& s) const override
    {
        double v = std::sin(0.5*theta(epsilon,s));
        return v*v;
    }

    // scattering angle calculated using Gauss-Mehler quadrature
    // (return NaN on error)
    double theta(double epsilon, double s) const;

    /* returns the cross section in the center of mass frame considering a screened potential */
    double crossSection(double E, unsigned int Z1, unsigned int Z2, double massRatio, double thetaCM) const;

protected:
    unsigned int NSUM;

    // use bisection (Newton) method to find the minimal approach distance x0
    // (return THETAERR on error)
    static double findX0(double epsilon, double s);
    // Following functions used to compute scattering angle following Gauss-Mehler quadrature
    static double H(double u, double x0, double epsilon, double s);
    // the root of this function (vs x0) is the minimal
    // approach distance at energy epsilon and impact parameter s
    static double funcX0(double x0, double epsilon, double s);
    // use bisection method to find the reduced
    // impact parameter s given epsilon and thetaCM
    // (return NaN on error)
    double finds(double epsilon, double thetaCM) const;
};

struct reducedXS_zbl_magic : public reducedXS::reducedXS_impl
{
    typedef screening<reducedXS::ZBL> screening_function;

    virtual reducedXS::IntegrationMethod integrationMethod() const override
    { return reducedXS::ZBL_Magick; }
    virtual reducedXS::ScreeningPotentialType potentialType() const override
    { return screening_function::type(); }

    static double theta(double epsilon, double s)
    {
        return 2*std::acos(MAGIC(epsilon,s));
    }
    virtual double screeningLength(int Z1, int Z2) const override
    {
        return screening_function::screeningLength(Z1, Z2);
    }
    virtual double screeningFunction(const double& x) const override
    {
        return screening_function::screeningFunction(x);
    }
    virtual double sin2Thetaby2(const double& epsilon, const double& s) const override
    {
        double m = MAGIC(epsilon,s);
        return 1. - m*m;
    }

    /**
     * @brief Implementation of MAGIC formula
     * 
     * @param B reduced impact par
     * @param epsilon reduced center of mass energy
     * @return double cos(theta/2) of the scattering event
     */
    static double MAGIC(const double& epsilon, const double& s);

    /**
     * @brief ZBL potential
     * 
     * Evaluate the ZBL potential and optionally dV/dR 
     * Implementation taken from ZBL85 
     * 
     * @param R reduced radius
     * @param Vprime if a non-NULL pointer is passed, it receives the value of dV/dR
     * @return double the value of the potential
     */
    static double ZBL_and_deri(double R, double* Vprime);
};

template<class _XS_impl>
class reducedXS_corteo4bit : public reducedXS::reducedXS_impl
{
    std::shared_ptr<_XS_impl> xs_impl_;

    /* scattering matrix: sin^2(theta_CM/2) */
    Array2Df  matrix_;

public:

    const static int rows = corteo4bit::e_index::dim;
    const static int cols = corteo4bit::s_index::dim;

    reducedXS_corteo4bit(unsigned int nsum = 100) :
        xs_impl_(new _XS_impl(nsum)),
        matrix_(rows,cols)
    {}
    reducedXS_corteo4bit(const reducedXS_corteo4bit& other) :
        xs_impl_(other.xs_impl_), matrix_(other.matrix_)
    {}

    virtual reducedXS::IntegrationMethod integrationMethod() const override
    { return reducedXS::Corteo4bitTable; }
    virtual reducedXS::ScreeningPotentialType potentialType() const override
    { return xs_impl_->potentialType(); }
    virtual double screeningLength(int Z1, int Z2) const override
    {
        return xs_impl_->screeningLength(Z1, Z2);
    }
    virtual double screeningFunction(const double& x) const override
    {
        return xs_impl_->screeningFunction(x);
    }
    virtual double sin2Thetaby2(const double& epsilon, const double& s) const override
    {
        corteo4bit::e_index ie(epsilon);
        corteo4bit::s_index is(s);
        return matrix_[ie][is];
    }
    virtual double sin2Thetaby2(int ie, int is) const override
    {
        return matrix_[ie][is];
    }

    const float* data() const { return matrix_.data(); }
    const float* operator[](int i) const { return matrix_[i]; }

    int size_epsilon() const { return corteo4bit::rows; }
    int size_s() const { return corteo4bit::cols; }
    int size() const
    { return corteo4bit::rows * corteo4bit::cols; }


    /* compute all the elements of the matrix
    user sets showProgress!=0 to display the progress of this (long) calculation to the console
    return 1 if successful, 0 if not able to write file */
    /*
     * compute all the elements of the xs matrix
     *
     * if a ostream pointer is passed, msgs are written showing the progress
     *
     * return 1 if successful, 0 if not able to write file
     *
     */
    virtual int init(std::ostream* os = NULL) override;

};

class scatteringXS
{
protected:

    float screening_length;    /* screening length from ZBL85,p45, eq.2-60, but in [nm] */
    float inv_screening_length;/* 1/screening length in 1/nm */
    float mass_ratio;          /* ... */
    float sqrt_mass_ratio;     /* we will need this occasionally */
    float kfactor_m;           /* mass part of the kinematic factor. This is called EC in the TRIM code */
    float red_E_conv;          /* reduced energy conversion factor. This is called FI in TRIM */

    reducedXS xs;

    Array2Df cosTable, sinTable;

    void fillCosSinTable();

public:
    explicit scatteringXS(const reducedXS& axs = reducedXS());

    float sqrtMassRatio() const { return sqrt_mass_ratio; }

    virtual void init(float Z1, float M1, float Z2, float M2);

    virtual void scatter(float erg, float P,
                         float &recoil_erg, float &sintheta, float &costheta) const;

};

#endif
