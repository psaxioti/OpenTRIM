#ifndef _XS_H_
#define _XS_H_

#include <cmath>
#include <vector>

#define BOHR_RADIUS 0.05291772108 /* Bohr radius in nm */
#define SCREENCONST 0.8853 /* screening length constant [A] */
#define E2 1.43996445      /* e^2 / 4 pi eps0 = e^2 c^2 in [eV][nm] */
#define AMUbyE 1.036426883E-8 /* amu/e = 1.660538782E-27/1.602176487E-19 */
#define ELEMENTARY_CHARGE 1.602176487E-19


struct screening_none {
    static double screeningLength(int Z1, int Z2) { return BOHR_RADIUS; } // Bohr radius in nm
    static double screeningFunction(double x) { return 1.; }
};
struct screening_lj {
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/std::sqrt(std::pow(Z1,2./3)+std::pow(Z2, 2./3));
    }
    static double screeningFunction(double x) {
        double y = 3.108*std::sqrt(x);
        return exp(-y)*(1.+y*(1.+y*(0.3344+y*(0.0485+2.647e-3*y))));
    }
};
struct screening_univ {
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Universal screening function coefficients TRIM85 */
    constexpr static const double C[] =
        {0.18175, 0.50986, 0.28022, 0.028171};
    constexpr static const double A[] =
        {-3.1998, -0.94229, -0.4029, -0.20162};

    static double screeningFunction(double x) {
        return C[0]*exp(A[0]*x)+C[1]*exp(A[1]*x)+C[2]*exp(A[2]*x)+C[3]*exp(A[3]*x);
    }
};

struct xs_base
{
    virtual double screeningLength(int Z1, int Z2) const = 0;
    virtual double screeningFunction(double x) const = 0;
    virtual double sin2Thetaby2(double epsilon, double s) const = 0;
};

/*
 * ion scattering cross section using Gauss-Mehler quadrature
 * see Yuan et al. NIMB1993
 */
template<class PHI = screening_univ>
struct xs_quad : public xs_base
{
    xs_quad(unsigned int nsum = 100) : NSUM(nsum)
    {}

    typedef PHI screening_function;

    virtual double screeningLength(int Z1, int Z2) const override
    {
        return PHI::screeningLength(Z1, Z2);
    }
    virtual double screeningFunction(double x) const override
    {
        return PHI::screeningFunction(x);
    }
    virtual double sin2Thetaby2(double epsilon, double s) const override
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

struct xs_zbl : public xs_base
{

    static double theta(double epsilon, double s)
    {
        return 2*std::acos(MAGIC(epsilon,s));
    }
    virtual double screeningLength(int Z1, int Z2) const override
    {
        return screening_univ::screeningLength(Z1, Z2);
    }
    virtual double screeningFunction(double x) const override
    {
        return screening_univ::screeningFunction(x);
    }
    virtual double sin2Thetaby2(double epsilon, double s) const override
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
    static double MAGIC(double epsilon, double B);

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

class scattering
{
protected:
    float screening_length;    /* screening length from ZBL85,p45, eq.2-60, but in [nm] */
    float inv_screening_length;/* 1/screening length in 1/nm */
    float mass_ratio;          /* ... */
    float sqrt_mass_ratio;     /* we will need this occasionally */
    float kfactor_m;           /* mass part of the kinematic factor. This is called EC in the TRIM code */
    float red_E_conv;          /* reduced energy conversion factor. This is called FI in TRIM */

    const xs_base& xs;

public:
    explicit scattering(const xs_base& axs) : xs(axs)
    {}

    virtual void init(float Z1, float M1, float Z2, float M2);

    virtual void scatter(float erg, float P,
                         float &recoil_erg, float &sintheta, float &costheta);

};

class corteo;

class scattering_corteo : public scattering
{

    const corteo& corteo_;
    std::vector<float> cosTable, sinTable;
    void fillCosSinTable();

public:

    explicit scattering_corteo(const xs_base& axs, const corteo& c);

    virtual void init(float Z1, float M1, float Z2, float M2) override;

    virtual void scatter(float erg, float P,
                         float &recoil_erg, float &sintheta, float &costheta) override;
};

#endif
