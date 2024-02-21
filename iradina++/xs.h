#ifndef _XS_H_
#define _XS_H_

#include <cmath>
#include <vector>
#include <array>

#include <Eigen/Dense>

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

    static double F(double x, double e, double s) {
        double sx = s/x;
        return 1.-PHI::screeningFunction(x)/(x*e)-sx*sx;
    }

    static double H(double u, double x0, double e, double s) {
        return std::sqrt((1-u*u)/F(x0/u,e,s));
    }

    // use bisection (Newton) method to find the minimal approach distance x0
    // (return THETAERR on error)
    static double minApproach(double e, double s) {
        double x2 = 1.0/(2.0*e);
        x2 = x2+sqrt(x2+s*s); // inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15
        double x1 = x2/10.;
        double f1 = F(x1, e, s);  // should be always negative
        double f2 = F(x2, e, s);  // should be always positive (starting ~1.0)

        if(f1>=0.) {
            // initial guess for x1 too optimistic, start with a safe value (much longer)
            x1 = 1.e-8;
            f1 = F(x1, e, s);  // should be always negative
        }

        assert(f1<0.0 && f2>0.0); // values should be on each side of 0

        double xm = 0.5*(x1+x2);
        double fm = F(xm, e, s);
        do {
            if (fm<0.) x1 = xm;
            else x2 = xm;
            xm = 0.5*(x1+x2);
            fm = F(xm, e, s);
        } while (fabs(fm) > (10*std::numeric_limits<double>::epsilon()) );
        // 1.-XXX wont be more precise than 2^(-52)=2e-16
        // but using 1e-10 saves time and is enough precise for float conversion later on

        return xm;
    }

    // scattering angle calculated using Gauss-Mehler quadrature
    // (return THETAERR on error)
    static double theta(double e, double s, int nsum = 100)
    {
        double x0 = minApproach(e, s);

        double sum(0.);
        int m = nsum >> 1;
        double a = M_PI/2/m;
        double b = a/2;
        for(unsigned int j=0; j<m; j++) {
            double uj = cos(a*j+b);
            sum += H(uj, x0, e, s);
        }

        return M_PI-2.0*s/x0*a*sum;
    }

    // use bisection method to find the reduced
    // impact parameter s given epsilon and thetaCM
    // (return THETAERR on error)
    static double findS(double e, double thetaCM) {

        // inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
        double gamma = (M_PI-thetaCM)/M_PI;
        double x0 = minApproach((1.0-gamma*gamma)*e, 1.e-8);
        double x1 = 0.7*gamma*x0;
        double x2 = 1.0/(2.0*e*tan(thetaCM/2.0));

        double f1 = thetaCM-theta(e, x1);  // should be always negative
        double f2 = thetaCM-theta(e, x2);  // should be always positive (starting ~1.0)

        assert(f1<0.0 && f2>0.0);    // values should be on each side of 0

        double xm = 0.5*(x1+x2);
        double fm = thetaCM-theta(e,xm);
        int k = 0;
        do {
            if (fm <0. ) x1 = xm;
            else x2 = xm;
            xm = 0.5*(x1+x2);
            fm = thetaCM-theta(e,xm);
            k++;
        } while (fabs(fm)>1e-6 && k<100);

        return xm;
    }

    /* returns the cross section in the center of mass frame considering a screened potential */
    static double crossSection(double e, double thetaCM)
    {
        // find corresponfding reduced impact parameter
        double s = findS(e, thetaCM);

        // ds/dTheta using five-point stencil
        double ds = s*0.001;
        double dsdTheta = (12.0*ds)/(-theta(e,s+2.0*ds)
                                         +8.0*theta(e,s+ds)
                                         -8.0*theta(e,s-ds)
                                         +theta(e,s-2.0*ds));

        // return "non-reduced" center-of-mass cross section
        return s/sin(thetaCM)*fabs(dsdTheta);
    }

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
    static double MAGIC(const double& e, const double& s)
    {
        double cost2;  /* cos(theta/2)*/
        double RoC,Delta,R,RR,A,G,alpha,beta,gamma,V,V1,FR,FR1,Q;
        double SQE;

        /* TRIM 85:  */
        static const double C[] = {0., 0.99229, 0.011615, 0.0071222, 14.813, 9.3066};

        /* Initial guess for R: */
        R=s;
        RR=-2.7*log(e*s);
        if(RR>=s){
            /*   if(RR<B) calc potential; */
            RR=-2.7*log(e*RR);
            if(RR>=s){
                R=RR;
            }
        }
        /* TRIM85: 330 */
        do{
            /* Calculate potential and its derivative */
            V = ZBL_and_deriv(R,&V1);
            FR  = s * s / R + V*R/e - R;
            FR1 = - s * s / (R * R) + (V+V1*R)/e - 1.0;
            Q   = FR/FR1;
            R   = R-Q;
        } while(fabs(Q/R)>0.001);

        RoC = -2.0 * (e-V)/V1;
        SQE = sqrt(e);

        alpha = 1+ C[1]/SQE;
        beta  = (C[2]+SQE) / (C[3]+SQE);           /* TRIM85: CC */
        gamma = (C[4]+e)/(C[5]+e);
        A     = 2*alpha*e*pow(s,beta);
        G     = gamma / ( sqrt((1.0+A*A))-A  );    /* TRIM85: 1/FF */
        Delta = A * (R-s)/(1+G);

        cost2=(s+RoC+Delta)/(R+RoC);
        return cost2;
    }

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
    static double ZBL_and_deriv(double R, double* Vprime)
    {
        auto &C =  screening_function::C;
        auto &A =  screening_function::A;

        double EX1 = C[0]*exp(-A[0]*R);
        double EX2 = C[1]*exp(-A[1]*R);
        double EX3 = C[2]*exp(-A[2]*R);
        double EX4 = C[3]*exp(-A[3]*R);

        double V=(EX1+EX2+EX3+EX4)/R;
        if (Vprime)
            *Vprime = -(V + A[0]*EX1 + A[1]*EX2 + A[2]*EX3 + A[3]*EX4)/R;

        return V;
    }

    // use bisection method to find the reduced
    // impact parameter s given epsilon and thetaCM
    // (return THETAERR on error)
    static double findS(double e, double thetaCM) {

        // inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
        double gamma = (M_PI-thetaCM)/M_PI;
        double x0 = XSquad< screeningZBL >::minApproach((1.0-gamma*gamma)*e, 1.e-8);
        double x1 = 0.7*gamma*x0;
        double x2 = 1.0/(2.0*e*tan(thetaCM/2.0));

        double f1 = thetaCM-theta(e, x1);  // should be always negative
        double f2 = thetaCM-theta(e, x2);  // should be always positive (starting ~1.0)

        assert(f1<0.0 && f2>0.0);    // values should be on each side of 0

        double xm = 0.5*(x1+x2);
        double fm = thetaCM-theta(e,xm);
        do {
            if (fm <0. ) x1 = xm;
            else x2 = xm;
            xm = 0.5*(x1+x2);
            fm = thetaCM-theta(e,xm);
        } while (fabs(fm)>1e-6);

        return xm;
    }
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
    XSlab() : P_()
    {}
    XSlab(const XSlab& x) : P_(x.P_)
    {}
    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
    void init(float Z1, float M1, float Z2, float M2) {
        P_.init<XS_zbl_magic>(Z1,M1,Z2,M2);
    }
    void scatter(float E, float S,
                 float &recoil_erg, float &sintheta, float &costheta) const
    {
        float sin2thetaby2 = XS_zbl_magic::sin2Thetaby2(E*P_.red_E_conv, S*P_.inv_screening_length);
        costheta = 1.f - 2*sin2thetaby2;
        sintheta = std::sqrt(1.f-costheta*costheta);
        /* now conversion to lab frame of reference: */
        float theta=atan( sintheta/(costheta+P_.mass_ratio) );
        /*if(isnan(theta)){theta=0;}*/
        sintheta=sin(theta);
        costheta=cos(theta);
        recoil_erg = E*P_.kfactor_m*sin2thetaby2;
    }
    float impactPar(float E, float T)
    {
        double thetaCM = T/E/P_.kfactor_m;
        thetaCM = 2.*std::asin(std::sqrt(thetaCM));
        return XS_zbl_magic::findS(E*P_.red_E_conv,thetaCM)*P_.screening_length;

    }
};

template<class _XS>
class xs_corteo_impl_ {
    typedef typename _XS::corteo_idx_t _My_t;
    cm_pars P_;
    std::array<float, _My_t::rows * _My_t::cols > sinTable, cosTable;
public:
    xs_corteo_impl_()
    {}
    xs_corteo_impl_(const xs_corteo_impl_& x) :
        P_(x.P_),
        sinTable(x.sinTable),
        cosTable(x.cosTable)
    {}
    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
    void init(float Z1, float M1, float Z2, float M2) {
        P_.init<_XS>(Z1,M1,Z2,M2);
        /* compute scattering angle components */
        double costhetaLab, sinthetaLab;
        double mr = P_.mass_ratio;
        for (typename _My_t::e_index ie; ie!=ie.end(); ie++)
            for (typename _My_t::s_index is; is!=is.end(); is++)
            {
                double s2 = _XS::sin2Thetaby2(ie,is);
                double costheta = 1.-2.*s2;

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
//        int k = _My_t::table_index(e*P_.red_E_conv, s*P_.inv_screening_length);
//        const float* p = _XS::data();
//        recoil_erg = e*P_.kfactor_m*p[k];
//        sintheta=sinTable[k];
//        costheta=cosTable[k];

        // bilinear interp
        const float* p = _XS::data();
        recoil_erg = e*P_.kfactor_m;
        Eigen::Vector4f sinth, costh, sin2thetaby2;
        e *= P_.red_E_conv;
        s *= P_.inv_screening_length;
        typename _My_t::e_index i1(e), i2(i1);
        typename _My_t::s_index j1(s), j2(j1);
        i2++; j2++;
        int k = i1 * _My_t::cols + j1;
        sinth[0] = sinTable[k]; costh[0] = cosTable[k]; sin2thetaby2[0] = p[k];
        k = i2 * _My_t::cols + j1;
        sinth[1] = sinTable[k]; costh[1] = cosTable[k]; sin2thetaby2[1] = p[k];
        k = i2 * _My_t::cols + j2;
        sinth[2] = sinTable[k]; costh[2] = cosTable[k]; sin2thetaby2[2] = p[k];
        k = i1 * _My_t::cols + j2;
        sinth[3] = sinTable[k]; costh[3] = cosTable[k]; sin2thetaby2[3] = p[k];
        float t = (e - *i1)/(*i2 - *i1);
        float u = (s - *j1)/(*j2 - *j1);
        Eigen::Vector4f coeff = { (1-t)*(1-u), t*(1-u), t*u, (1-t)*u };
        sintheta = coeff.dot(sinth);
        costheta = coeff.dot(costh);
        recoil_erg *= coeff.dot(sin2thetaby2);
    }
    float impactPar(float E, float T)
    {
        double thetaCM = T/E/P_.kfactor_m;
        thetaCM = 2.*std::asin(std::sqrt(thetaCM));
        return XSquad< screeningZBL >::findS(E*P_.red_E_conv,thetaCM)*P_.screening_length;

    }
};

template<>
class XSlab<XS_corteo4bit> : public xs_corteo_impl_< XS_corteo4bit >
{};

template<>
class XSlab<XS_corteo6bit> : public xs_corteo_impl_< XS_corteo6bit >
{};

#endif
