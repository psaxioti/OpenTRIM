#ifndef _XS_H_
#define _XS_H_

#include <cmath>
#include <vector>
#include <array>

#include "corteo.h"
#include "corteo4bit.h"
#include "corteo6bit.h"

#define BOHR_RADIUS 0.05291772108 /* Bohr radius [nm] */
#define SCREENCONST 0.88534 /* Lindhard screening length constant */
#define E2 1.43996445      /* e^2 / 4 pi eps0 = e^2 c^2 in [eV][nm] */
#define AMUbyE 1.036426883E-8 /* amu/e = 1.660538782E-27/1.602176487E-19 */
#define ELEMENTARY_CHARGE 1.602176487E-19

//const char* names[] = {
//    "Unscreened Coulomb",
//    "Moliere",
//    "Kr-C",
//    "Lenz-Jensen",
//    "Ziegler-Biersack-Littmark (ZBL)"};

struct screeningNone {
    static double screeningLength(int Z1, int Z2) { return BOHR_RADIUS; } // Bohr radius in nm
    static double screeningFunction(const double& x) { return 1.; }
};

struct screeningLenzJensen
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/std::sqrt(std::pow(Z1,2./3)+std::pow(Z2, 2./3));
    }
    static double screeningFunction(const double& x) {
        double y = 3.108*std::sqrt(x);
        return exp(-y)*(1.+y*(1.+y*(0.3344+y*(0.0485+2.647e-3*y))));
    }
};
struct screeningKrC
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
struct screeningMoliere
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
struct screeningZBL
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
        return C[0]*exp(-A[0]*x)+
               C[1]*exp(-A[1]*x)+
               C[2]*exp(-A[2]*x)+
               C[3]*exp(-A[3]*x);
    }
};

/*
 * ion scattering cross section using Gauss-Mehler quadrature
 * see Yuan et al. NIMB1993
 */
template<class PHI = screeningZBL >
struct XSquad
{
    XSquad(unsigned int nsum = 100) : NSUM(nsum)
    {}

    typedef PHI screening_function;

    static double screeningLength(int Z1, int Z2)
    {
        return PHI::screeningLength(Z1, Z2);
    }
    static double screeningFunction(const double& x)
    {
        return PHI::screeningFunction(x);
    }
    double sin2Thetaby2(const double& e, const double& s)
    {
        double v = std::sin(0.5*theta(e,s));
        return v*v;
    }

    // scattering angle calculated using Gauss-Mehler quadrature
    // (return NaN on error)
    double theta(double e, double s) const;

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

struct XS_zbl_magic
{
    typedef screeningZBL screening_function;

    static double theta(double epsilon, double s)
    {
        return 2*std::acos(MAGIC(epsilon,s));
    }
    static double screeningLength(int Z1, int Z2)
    {
        return screening_function::screeningLength(Z1, Z2);
    }
    static double screeningFunction(const double& x)
    {
        return screening_function::screeningFunction(x);
    }
    static double sin2Thetaby2(const double& epsilon, const double& s)
    {
        double m = MAGIC(epsilon,s);
        return 1. - m*m;
    }

    /**
     * @brief Implementation of MAGIC formula
     *
     * @param e is the reduced energy in center of mass system
     * @param s is the reduced impact parameter
     * @return double cos(theta/2) of the scattering event
     */
    static double MAGIC(const double& e, const double& s);

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

struct XS_corteo4bit
{
    typedef screeningZBL screening_function;
    typedef corteo4bit corteo_idx_t;

    const static int rows = corteo4bit::e_index::dim;
    const static int cols = corteo4bit::s_index::dim;

    static double screeningLength(int Z1, int Z2)
    {
        return screening_function::screeningLength(Z1, Z2);
    }
    static double screeningFunction(const double& x)
    {
        return screening_function::screeningFunction(x);
    }
    static double sin2Thetaby2(const double& e, const double& s)
    {
        const float* p = corteo4bitdata();
        int i = corteo4bit::table_index(e,s);
        return p[i];
    }
    static float sin2Thetaby2(int ie, int is)
    {
        const float* p = corteo4bitdata();
        return p[ie*cols + is];
    }
    static const float* data() { return corteo4bitdata(); }
};

struct XS_corteo6bit
{
    typedef screeningZBL screening_function;
    typedef corteo6bit corteo_idx_t;

    const static int rows = corteo6bit::e_index::dim;
    const static int cols = corteo6bit::s_index::dim;

    static double screeningLength(int Z1, int Z2)
    {
        return screening_function::screeningLength(Z1, Z2);
    }
    static double screeningFunction(const double& x)
    {
        return screening_function::screeningFunction(x);
    }
    static double sin2Thetaby2(const double& e, const double& s)
    {
        const float* p = corteo6bitdata();
        int i = corteo6bit::table_index(e,s);
        return p[i];
    }
    static float sin2Thetaby2(int ie, int is)
    {
        const float* p = corteo6bitdata();
        return p[ie*cols + is];
    }
    static const float* data() { return corteo6bitdata(); }
};


struct cm_pars {
    float screening_length;    /* screening length from ZBL85,p45, eq.2-60, but in [nm] */
    float inv_screening_length;/* 1/screening length in 1/nm */
    float mass_ratio;          /* ... */
    float sqrt_mass_ratio;     /* we will need this occasionally */
    float kfactor_m;           /* mass part of the kinematic factor. This is called EC in the TRIM code */
    float red_E_conv;          /* reduced energy conversion factor. This is called FI in TRIM */

    template<class _XScm>
    void init(float Z1, float M1, float Z2, float M2)
    {
        /* Adapted from the corteo code */
        screening_length     = _XScm::screeningLength(Z1,Z2); /* This is in nm! */
        inv_screening_length = 1.0f/screening_length;
        mass_ratio           = M1/M2;
        sqrt_mass_ratio      = std::sqrt(mass_ratio);
        kfactor_m            = 4*mass_ratio / ((mass_ratio+1) * (mass_ratio+1)); /* k factor without the angle part */
        red_E_conv           = screening_length / ((mass_ratio+1) * Z1 * Z2 * E2);
    }
};

template<class _XScm> class XSlab;

template<>
class XSlab<XS_zbl_magic>
{
    cm_pars P_;
public:
    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
    void init(float Z1, float M1, float Z2, float M2) {
        P_.init<XS_zbl_magic>(Z1,M1,Z2,M2);
    }
    void scatter(float e, float s,
                 float &recoil_erg, float &sintheta, float &costheta) const
    {
        float sin2thetaby2 = XS_zbl_magic::sin2Thetaby2(e*P_.red_E_conv, s*P_.inv_screening_length);
        costheta = 1.f - 2*sin2thetaby2;
        sintheta = std::sqrt(1.f-costheta*costheta);
        /* now conversion to lab frame of reference: */
        float theta=atan( sintheta/(costheta+P_.mass_ratio) );
        /*if(isnan(theta)){theta=0;}*/
        sintheta=sin(theta);
        costheta=cos(theta);
        recoil_erg = e*P_.kfactor_m*sin2thetaby2;
    }
};

template<class _XS>
class xs_corteo_impl_ {
    typedef typename _XS::corteo_idx_t _My_t;
    cm_pars P_;
    std::array<float, _My_t::rows * _My_t::cols > sinTable, cosTable;
public:
    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
    void init(float Z1, float M1, float Z2, float M2) {
        P_.init<_XS>(Z1,M1,Z2,M2);
        /* compute scattering angle components */
        double sin2thetaby2, costheta, costhetaLab, sinthetaLab;
        double mr = P_.mass_ratio;
        for (typename _My_t::e_index ie; ie!=ie.end(); ie++)
            for (typename _My_t::s_index is; is!=is.end(); is++)
            {
                sin2thetaby2 = _XS::sin2Thetaby2(ie,is);
                costheta = 1.-2.*sin2thetaby2;

                if(costheta==-1.0 && mr==1.0) {
                    costhetaLab = 0.0;  /* peculiar case of head-on collision of identical masses */
                    sinthetaLab = 1.0;
                } else {
                    costhetaLab = (costheta+mr)/sqrt(1.+2.*mr*costheta+mr*mr);
                    sinthetaLab = sqrt(1.-costhetaLab*costhetaLab);
                }

                int k = ie*_My_t::cols + is;
                cosTable[k] = costhetaLab;
                sinTable[k] = sinthetaLab;

                /* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
                /* In some rare cases, when cos=1, then sin becomes "Not a Number". To prevent this, I will set the sine to 0 in those cases. */
                if( std::isnan(sinTable[k]) ) {
                    cosTable[k]=0.f;
                    sinTable[k]=1.f;
                }
            }
    }
    void scatter(float e, float s,
                 float &recoil_erg, float &sintheta, float &costheta) const
    {
        int k = _My_t::table_index(e*P_.red_E_conv, s*P_.inv_screening_length);
        const float* p = _XS::data();
        float sin2thetaby2 = p[k];
        recoil_erg = e*P_.kfactor_m*sin2thetaby2;
        sintheta=sinTable[k];
        costheta=cosTable[k];
    }
};

template<>
class XSlab<XS_corteo4bit> : public xs_corteo_impl_< XS_corteo4bit >
{};

template<>
class XSlab<XS_corteo6bit> : public xs_corteo_impl_< XS_corteo6bit >
{};

//template<>
//class XSlab<XS_corteo4bit>
//{
//    cm_pars P_;
//    std::array<float, corteo4bit::rows*corteo4bit::cols > sinTable, cosTable;
//public:
//    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
//    void init(float Z1, float M1, float Z2, float M2) {
//        P_.init<XS_corteo4bit>(Z1,M1,Z2,M2);
//        /* compute scattering angle components */
//        double sin2thetaby2, costheta, costhetaLab, sinthetaLab;
//        double mr = P_.mass_ratio;
//        for (corteo4bit::e_index ie; ie!=ie.end(); ie++)
//            for (corteo4bit::s_index is; is!=is.end(); is++)
//            {
//                sin2thetaby2 = XS_corteo4bit::sin2Thetaby2(ie,is);
//                costheta = 1.-2.*sin2thetaby2;

//                if(costheta==-1.0 && mr==1.0) {
//                    costhetaLab = 0.0;  /* peculiar case of head-on collision of identical masses */
//                    sinthetaLab = 1.0;
//                } else {
//                    costhetaLab = (costheta+mr)/sqrt(1.+2.*mr*costheta+mr*mr);
//                    sinthetaLab = sqrt(1.-costhetaLab*costhetaLab);
//                }

//                int k = ie*corteo4bit::cols + is;
//                cosTable[k] = costhetaLab;
//                sinTable[k] = sinthetaLab;

//                /* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
//                /* In some rare cases, when cos=1, then sin becomes "Not a Number". To prevent this, I will set the sine to 0 in those cases. */
//                if( std::isnan(sinTable[k]) ) {
//                    cosTable[k]=0.f;
//                    sinTable[k]=1.f;
//                }
//            }
//    }
//    void scatter(float e, float s,
//                 float &recoil_erg, float &sintheta, float &costheta) const
//    {
//        int k = corteo4bit::table_index(e*P_.red_E_conv, s*P_.inv_screening_length);
//        const float* p = corteo6bitdata();
//        float sin2thetaby2 = p[k];
//        recoil_erg = e*P_.kfactor_m*sin2thetaby2;
//        sintheta=sinTable[k];
//        costheta=cosTable[k];
//    }

//};

//template<>
//class XSlab<XS_corteo6bit>
//{
//    cm_pars P_;
//    std::array<float, corteo6bit::rows*corteo6bit::cols > sinTable, cosTable;
//public:
//    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
//    void init(float Z1, float M1, float Z2, float M2) {
//        P_.init<XS_corteo6bit>(Z1,M1,Z2,M2);
//        /* compute scattering angle components */
//        double sin2thetaby2, costheta, costhetaLab, sinthetaLab;
//        double mr = P_.mass_ratio;
//        for (corteo6bit::e_index ie; ie!=ie.end(); ie++)
//            for (corteo6bit::s_index is; is!=is.end(); is++)
//            {
//                sin2thetaby2 = XS_corteo6bit::sin2Thetaby2(ie,is);
//                costheta = 1.-2.*sin2thetaby2;

//                if(costheta==-1.0 && mr==1.0) {
//                    costhetaLab = 0.0;  /* peculiar case of head-on collision of identical masses */
//                    sinthetaLab = 1.0;
//                } else {
//                    costhetaLab = (costheta+mr)/sqrt(1.+2.*mr*costheta+mr*mr);
//                    sinthetaLab = sqrt(1.-costhetaLab*costhetaLab);
//                }

//                int k = ie*corteo6bit::cols + is;
//                cosTable[k] = costhetaLab;
//                sinTable[k] = sinthetaLab;

//                /* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
//                /* In some rare cases, when cos=1, then sin becomes "Not a Number". To prevent this, I will set the sine to 0 in those cases. */
//                if( std::isnan(sinTable[k]) ) {
//                    cosTable[k]=0.f;
//                    sinTable[k]=1.f;
//                }
//            }
//    }
//    void scatter(float e, float s,
//                 float &recoil_erg, float &sintheta, float &costheta) const
//    {
//        int k = corteo6bit::table_index(e*P_.red_E_conv, s*P_.inv_screening_length);
//        const float* p = corteo6bitdata();
//        float sin2thetaby2 = p[k];
//        recoil_erg = e*P_.kfactor_m*sin2thetaby2;
//        sintheta=sinTable[k];
//        costheta=cosTable[k];
//    }

//};

#endif
