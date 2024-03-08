#ifndef _XS_H_
#define _XS_H_

#include <cmath>
#include <vector>
#include <array>

#include <Eigen/Dense>

#include "corteo.h"

#define BOHR_RADIUS 0.05291772108 /* Bohr radius [nm] */
#define SCREENCONST 0.88534 /* Lindhard screening length constant */
#define E2 1.43996445      /* e^2 / 4 pi eps0 = e^2 c^2 in [eV][nm] */
#define AMUbyE 1.036426883E-8 /* amu/e = 1.660538782E-27/1.602176487E-19 */
#define ELEMENTARY_CHARGE 1.602176487E-19

/**
 * \defgroup XS Screened Coulomb potentials and cross-sections
 * @{
 *
 * A set of C++ objects for classical scattering calculations of screened
 * Coulomb interatomic interactions.
 *
 * The general form of the potential is
 * \f[
 * V(r) = \frac{Z_1 Z_2 e^2}{r} \Phi(r/a)
 * \f]
 * where \f$ \Phi \f$ is the screening function and \f$ a \f$ the screening length.
 *
 * Different types of screening functions used for ion scattering are coded, see \ref ScreeningF.
 *
 * The scattering angle in the center-of-mass system can be obtained by
 * the scattering integral
 * \f[
 * \theta = \pi - 2 s \int_{x_0}^\infty {x^{-2}F^{-1/2}(x)\,dx}
 * \f]
 * where
 * \f[
 * F(x) = 1 - \frac{\Phi(x)}{x\,\epsilon} - \frac{s^2}{x^2}
 * \f]
 * and \f$ x_0 \f$ is the distance of closest approach which satisfies
 * \f$ F(x_0)=0 \f$.
 *
 * \f$  \epsilon = E_{CM} a/Z_1 Z_2 e^2\f$ and \f$ s = p/a\f$ are the
 * reduced center-of-mass energy and impact parameter.
 *
 * The integral can be evaluated by quadrature, which is performed by the
 * \ref xs_quad class.
 *
 * As quadrature is a costly operation, the scattering integrals are typically tabulated
 * for use in Monte-Carlo codes.
 * Here we adapt the tabulation method of the program Corteo, see \ref CorteoIdx.
 *
 * Alternatively, SRIM's MAGIC approximation for the scattering integrals of the ZBL potential
 * is also implemented
 * for comparison in the class \ref xs_zbl_magic.
 *
 * Ref.: Yuan et al. NIMB1993
 */

/**
 * @brief Screening function types
 */
typedef enum {
    ScreeningNone,          /**< Unscreened Coulomb potential */
    ScreeningLenzJensen,    /**< Lenz-Jensen */
    ScreeningKrC,           /**< Kr-C */
    ScreeningMoliere,       /**< Moliere */
    ScreeningZBL            /**< Ziegler-Biersack-Littmark (ZBL) Universal */
} screening_function_t;

/**
 * @brief Screening function definition
 *
 * A templated structure that includes the definition of the screening
 * function (screening_function::operator()) and the screening length
 * (screening_function::screeningLength).
 *
 * Additionally, it returns the name and type of screening as a string or enum, respectively.
 *
 * @tparam ScreeningType (screening_function_t) the type of screening function
 */
template<screening_function_t ScreeningType>
struct screening_function
{
    /**
     * @brief Screening length
     * @param Z1 is the atomic number of the projectile
     * @param Z2 is the atomic number of the target
     * @return the screening length [nm]
     */
    static double screeningLength(int Z1, int Z2);

    /**
     * @brief operator () implements the actual screening function
     * @param x is the radial distance in units of the screening length
     * @return the value of \f$ \Phi(x) \f$
     */
    double operator()(const double& x) const;

    /// The name of the screening function
    static const char* screeningName();

    /// The type of screening as a screening_function_t enum
    static screening_function_t screeningType() { return ScreeningType; }
};

/**@}*/  // XS

// Explicit specialization for the unscreened Coulomb interaction
template<>
struct screening_function< ScreeningNone >
{
    static double screeningLength(int Z1, int Z2) { return BOHR_RADIUS; } // Bohr radius in nm
    double operator()(const double& x) const { return 1.; }
    static const char* screeningName() { return "Unscreened Coulomb"; }
    static screening_function_t screeningType() { return ScreeningNone; }
};
// Explicit specialization for the Lenz-Jensen potential
template<>
struct screening_function< ScreeningLenzJensen >
{
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/std::sqrt(std::pow(Z1,2./3)+std::pow(Z2, 2./3));
    }
    double operator()(const double& x) const {
        double y = 3.108*std::sqrt(x);
        return exp(-y)*(1.+y*(1.+y*(0.3344+y*(0.0485+2.647e-3*y))));
    }
    static const char* screeningName() { return "Lenz-Jensen"; }
    static screening_function_t screeningType() { return ScreeningLenzJensen; }
};
// Explicit specialization for the Kr-C potential
template<>
struct screening_function< ScreeningKrC >
{
    static screening_function_t screeningType() { return ScreeningKrC; }
    static const char* screeningName() { return "Kr-C"; }
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Screening function coefficients */
    constexpr static const double C[] =
        {0.190945, 0.473674, 0.335381};
    constexpr static const double A[] =
        {0.278544, 0.637174, 1.919249};


    double operator()(const double& x) const {
        return C[0]*exp(-A[0]*x)+C[1]*exp(-A[1]*x)+C[2]*exp(-A[2]*x);
    }
};
// Explicit specialization for the Moliere potential
template<>
struct screening_function< ScreeningMoliere >
{
    static screening_function_t screeningType() { return ScreeningMoliere; }
    static const char* screeningName() { return "Moliere"; }
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Screening function coefficients */
    constexpr static const double C[] =
        {0.35, 0.55, 0.10};
    constexpr static const double A[] =
        {0.30, 1.20, 6.00};

    double operator()(const double& x) const {
        return C[0]*exp(-A[0]*x)+C[1]*exp(-A[1]*x)+C[2]*exp(-A[2]*x);
    }
};
// Explicit specialization for the Ziegler-Biersack-Littmark (ZBL) potential
template<>
struct screening_function< ScreeningZBL >
{
    static screening_function_t screeningType() { return ScreeningZBL; }
    static const char* screeningName() { return "Ziegler-Biersack-Littmark (ZBL)"; }
    static double screeningLength(int Z1, int Z2) {
        return SCREENCONST*BOHR_RADIUS/(std::pow(Z1,0.23)+std::pow(Z2,  0.23));
    }

    /* Universal screening function coefficients TRIM85 */
    constexpr static const double C[] =
        {0.18175, 0.50986, 0.28022, 0.028171};
    constexpr static const double A[] =
        {3.19980, 0.94229, 0.40290, 0.201620};

    double operator()(const double& x) const {
        return C[0]*exp(-A[0]*x)+
               C[1]*exp(-A[1]*x)+
               C[2]*exp(-A[2]*x)+
               C[3]*exp(-A[3]*x);
    }
};

/**
 * @brief The xs_base class forms the basis of all cross-section classes
 *
 * @tparam
 *   ScreeningType the type of screening function
 *
 * @ingroup XS
 */
template<screening_function_t ScreeningType>
struct xs_base : public screening_function< ScreeningType >
{
    typedef screening_function< ScreeningType > Phi;

    /**
     * @brief The function \f$ F(x) \f$ under the scattering integral
     *
     * \f[
     * F(x) = 1 - \frac{\Phi(x)}{x\,e} - \frac{s^2}{x^2}
     * \f]
     *
     * @param x integral variable
     * @param e reduced energy
     * @param s reduced impact parameter
     * @return the value of the function
     */
    static double F(double x, double e, double s) {
        double sx = s/x;
        Phi p;
        return 1. - p(x)/(x*e)-sx*sx;
    }

    /**
     * @brief Finds the minimal approach distance of a scattered projectile
     *
     * For the unscreened Coulomb potential the minimal approach is
     * \f[
     * x_0 = \frac{1}{2\epsilon} \left[ 1 + \sqrt{1 + (2\,\epsilon\,s)^2} \right]
     * \f]
     *
     * For screened potentials \f$ x_0 \f$ is found by numerically solving
     * \f$ F(x_0)=0 \f$
     * using bysection.
     *
     * @param e reduced energy of the scattered particle
     * @param s reduced impact parameter
     * @return the minimal approach distance (in screening length units)
     */
    static double minApproach(double e, double s) {
        double x2 = 1.0/(2.0*e);
        x2 = x2+sqrt(x2*x2+s*s); // inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15

        double x1 = x2/10.;
        double f1 = F(x1, e, s);  // should be always negative
        double f2 = F(x2, e, s);  // should be always positive (starting ~1.0)

        if(f1>=0.) {
            // initial guess for x1 too optimistic, start with a safe value (much longer)
            x1 = 1.e-8;
            f1 = F(x1, e, s);  // should be always negative
        }
        if (f2 <=0. ) {
            do {
                x2 *= 1.001;
                f2 = F(x2, e, s);
            } while (f2 <= 0.);
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

    /**
     * @brief Scattering angle in the Impulse approximation
     *
     * From Lehmann & Leibfried, Z. Physic 172 (1963) 465, the 1st order term
     * in the "momentum approx." for the screened potential \f$ V(r) = (A/r)e^{-r/a} \f$ is
     *
     * \f[
     * \theta_1 = \epsilon^{-1}\,K_1(s)
     * \f]
     *
     * where \f$ K_1(x) \f$ is the modified Hankel function (or modified
     * Bessel function of 2nd kind).
     *
     * The above formula can be generalized for a screening functions
     * expressed as a sum of exponentials, i.e., for the ZBL, Kr-C & Moliere potentials.
     *
     * In all other cases the function returns the unscreened Coulomb
     * exact scattering angle
     *
     * \f[
     * \theta = 2 \arcsin \left[ 1 + (2\epsilon\, s)^2\right]^{-1/2}
     * \f]
     *
     * @param e reduced energy of the scattered particle
     * @param s reduced impact factor
     * @return the scattering angle [rad]
     */
    static double theta_impulse_approx(double e, double s) {
        // exact unscreened Coulomb theta
        double x = 2*e*s;
        return 2*std::asin(std::sqrt(1./(1+x*x)));
    }
};

// Unscreened Coulomb minimal approach distance
template<>
inline double xs_base<ScreeningNone>::minApproach(double e, double s)
{
    double x0 = 1.0/(2*e);
    return x0+std::sqrt(x0*x0+s*s);
}
// ZBL potential impulse aprox
template<>
inline double xs_base<ScreeningZBL>::theta_impulse_approx(double e, double s)
{
    auto &C =  Phi::C;
    auto &A =  Phi::A;
    double th = C[0]*A[0]*std::cyl_bessel_k(1,A[0]*s) +
                C[1]*A[1]*std::cyl_bessel_k(1,A[1]*s) +
                C[2]*A[2]*std::cyl_bessel_k(1,A[2]*s) +
                C[3]*A[3]*std::cyl_bessel_k(1,A[3]*s);
    return th/e;
}
// KrC potential impulse aprox
template<>
inline double xs_base<ScreeningKrC>::theta_impulse_approx(double e, double s)
{
    auto &C =  Phi::C;
    auto &A =  Phi::A;
    double th = C[0]*A[0]*std::cyl_bessel_k(1,A[0]*s) +
                C[1]*A[1]*std::cyl_bessel_k(1,A[1]*s) +
                C[2]*A[2]*std::cyl_bessel_k(1,A[2]*s);
    return th/e;
}
// Moliere potential impulse aprox
template<>
inline double xs_base<ScreeningMoliere>::theta_impulse_approx(double e, double s)
{
    auto &C =  Phi::C;
    auto &A =  Phi::A;
    double th = C[0]*A[0]*std::cyl_bessel_k(1,A[0]*s) +
                C[1]*A[1]*std::cyl_bessel_k(1,A[1]*s) +
                C[2]*A[2]*std::cyl_bessel_k(1,A[2]*s);
    return th/e;
}

/**
 * @brief The xs_quad class implements screened potential scattering integral evaluation by quadrature
 *
 * Using Gauss–Chebyshev quadrature (https://dlmf.nist.gov/3.5#v) the integral is evaluated as
 *
 * \f[
 *   \theta = \pi - \frac{2\,s}{x_0} \frac{\pi}{N}
 *   \sum_{j=0}^{N/2-1}{H\left[ \cos\left( \frac{\pi}{N}\,j + \frac{\pi}{2N} \right) \right]}
 * \f]
 * where
 * \f[
 * H(u) = \sqrt{\frac{1-u^2}{F(x_0/u)}}
 * \f]
 *
 * @tparam
 *   ScreeningType specifies the type of screening function
 *
 * @ingroup XS
 */
template<screening_function_t ScreeningType = ScreeningZBL>
struct xs_quad : public xs_base< ScreeningType >
{
    using xs_base< ScreeningType >::F;
    using xs_base< ScreeningType >::minApproach;
    using xs_base< ScreeningType >::theta_impulse_approx;

    /**
     * @brief Return \f$ \sin^2(\theta/2) \f$
     * @param e the reduced energy
     * @param s the reduced impact parameter
     * @return the value of \f$ \sin^2(\theta/2) \f$
     */
    double sin2Thetaby2(const double& e, const double& s)
    {
        double v = std::sin(0.5*theta(e,s));
        return v*v;
    }

    /**
     * @brief The function \f$ H(u) \f$ used in the quadrature sum
     *
     * \f[
     * H(u) = \sqrt{\frac{1-u^2}{F(x_0/u}}
     * \f]
     *
     * @param u function argument
     * @param x0 closest approach distance
     * @param e reduced energy
     * @param s reduced impact parameter
     * @return the value of H
     */
    static double H(double u, double x0, double e, double s) {
        return std::sqrt((1-u*u) / F(x0/u,e,s));
    }

    /**
     * @brief Scattering angle in center-of-mass (CM) system
     *
     * Calculated using Gauss-Mehler quadrature.
     *
     * \f[
     *   \theta = \pi - \frac{2\,s}{x_0} \frac{\pi}{N}
     *   \sum_{j=0}^{N/2-1}{H\left[ \cos\left( \frac{\pi}{N}\,j + \frac{\pi}{2N} \right) \right]}
     * \f]
     *
     * For \f$ s\cdot \epsilon^{1/6} > 100 \f$ where the scattering angle is very small
     * (\f$ \theta < 10^{-8} \f$) this function
     * returns the impulse approx xs_base::theta_impulse_approx
     *
     *
     * @param e reduced energy of the scattered particle in CM system
     * @param s reduced impact parameter
     * @param nsum number of terms in the Gauss-Mehler sum \f$ N \f$
     * @return
     */
    static double theta(double e, double s, int nsum = 100)
    {    
        // if s > 100/e^(1/6) use impulse approx
        // The limit was found empirically for theta < 1e-9
        double s3 = s*s*s;
        if (e*s3*s3 > 1.e12)
            return theta_impulse_approx(e,s);

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

    /**
     * @brief Return the reduced impact parameter
     *
     * Use bisection to find the reduced impact parameter s
     * given energy and scattering angle
     *
     * @param e reduced energy
     * @param thetaCM scattering angle (rad) in center-of-mass system
     * @return the reduced impact parameter
     */
    static double findS(double e, double thetaCM) {

        // inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
        double gamma = (M_PI-thetaCM)/M_PI;
        double x0 = minApproach((1.0-gamma*gamma)*e, 1.e-8);
        double x1, x2;
        if (e >= 1.) {
            x1 = 0.7*gamma*x0;
            x2 = 1.0/(2.0*e*tan(thetaCM/2.0));
        } else {
            x1 = 0.9*gamma*x0;
            x2 = 1.4*gamma*x0;
        }

        assert(thetaCM-theta(e, x1)<0.0 && thetaCM-theta(e, x2)>0.0);    // values should be on each side of 0

        double xm = 0.5*(x1+x2);
        double fm = thetaCM-theta(e,xm);
        int k = 0;
        do {
            if (fm <0. ) x1 = xm;
            else x2 = xm;
            xm = 0.5*(x1+x2);
            fm = thetaCM-theta(e,xm);
            k++;
        } while (fabs(fm)>1e-9 && k<100); //1e-6

        return xm;
    }

    /**
     * @brief Differential cross-section in center-of-mass system
     *
     * \f[
     *   \frac{d\sigma}{d\Omega} = \frac{s}{\sin(\theta)}\left| \frac{ds}{d\theta}\right|
     * \f]
     *
     * In units of \f$ a^2 \f$
     *
     * @param e is the reduced energy
     * @param thetaCM is the center-of-mass scattering angle (rad)
     * @return the reduced cross-section
     */
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


        return s/sin(thetaCM)*fabs(dsdTheta);
    }
};

/**
 * @brief Screened potential cross-section by the MAGIC formula
 *
 * Calculate the scattering integral for the scattering angle in
 * center-of-mass system and the cross-section for the
 * Ziegler-Biersack-Littmark (ZBL) Universal screening function
 * using the MAGIC interpolation formula of Biersack & Haggmark.
 *
 * Ref.: Biersack & Haggmark NIM1980
 *
 * @ingroup XS
 */
struct xs_zbl_magic : public xs_base< ScreeningZBL >
{
    using xs_base<ScreeningZBL>::minApproach;
    /**
     * @brief Scattering angle in center-of-mass (CM) system
     *
     * @param e reduced energy of the scattered particle in CM system
     * @param s reduced impact parameter
     * @return the scattering angle in rad
     */
    static double theta(double e, double s)
    {
        return 2*std::acos(MAGIC(e,s));
    }

    /**
     * @brief Return \f$ \sin^2(\theta/2) \f$
     * @param e the reduced energy
     * @param s the reduced impact parameter
     * @return the value of \f$ \sin^2(\theta/2) \f$
     */
    static double sin2Thetaby2(const double& e, const double& s)
    {
        double m = MAGIC(e,s);
        return 1. - m*m;
    }

    /**
     * @brief Implementation of MAGIC formula
     *
     * Returns \f$ \cos(\theta/2) \f$ where \f$ \theta \f$ is the center-of-mass
     * scattering angle
     *
     * @param e is the reduced energy in center of mass system
     * @param s is the reduced impact parameter
     * @return \f$ \cos(\theta/2) \f$
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
     * Evaluate the ZBL potential and optionally the 1st derivative dV/dR
     *
     * Implementation taken from ZBL85 
     * 
     * @param R is the reduced radius
     * @param Vprime if a non-NULL pointer is passed, it receives the value of dV/dR
     * @return the value of the potential
     */
    static double ZBL_and_deriv(double R, double* Vprime)
    {
        auto &C =  screening_function<ScreeningZBL>::C;
        auto &A =  screening_function<ScreeningZBL>::A;

        double EX1 = C[0]*exp(-A[0]*R);
        double EX2 = C[1]*exp(-A[1]*R);
        double EX3 = C[2]*exp(-A[2]*R);
        double EX4 = C[3]*exp(-A[3]*R);

        double V=(EX1+EX2+EX3+EX4)/R;
        if (Vprime)
            *Vprime = -(V + A[0]*EX1 + A[1]*EX2 + A[2]*EX3 + A[3]*EX4)/R;

        return V;
    }

    /**
     * @brief Return the reduced impact parameter
     *
     * Use bisection to find the reduced impact parameter s
     * given energy and scattering angle
     *
     * @param e reduced energy
     * @param thetaCM scattering angle (rad) in center-of-mass system
     * @return the reduced impact parameter
     */
    static double findS(double e, double thetaCM) {

        // inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
        double gamma = (M_PI-thetaCM)/M_PI;
        double x0 = minApproach((1.0-gamma*gamma)*e, 1.e-8);
        double x1 = 0.7*gamma*x0;
        double x2 = 1.0/(2.0*e*tan(thetaCM/2.0));

        // values for x1, x2 should be on either side of 0
        assert(thetaCM-theta(e, x1)<0.0 && thetaCM-theta(e, x2)>0.0);

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

    /**
     * @brief Differential cross-section in center-of-mass system
     *
     * \f[
     *   \frac{d\sigma}{d\Omega} = \frac{s}{\sin(\theta)}\left| \frac{ds}{d\theta}\right|
     * \f]
     *
     * In units of \f$ a^2 \f$
     *
     * @param e is the reduced energy
     * @param thetaCM is the center-of-mass scattering angle (rad)
     * @return the reduced cross-section
     */
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


        return s/sin(thetaCM)*fabs(dsdTheta);
    }
};

const float* corteo4bitdata();
// const float* corteo6bitdata();

/**
 * @brief 4-bit corteo-tabulated ZBL screened potential scattering integral
 *
 * A 2-dimensional table of \f$ \sin^2\theta/2 \f$ as a function of reduced energy
 * and reduced impact factor, where \f$ \theta \f$ is
 * the center-of-mass scattering angle for the ZBL screened potential
 * pre-calculated by Gauss–Chebyshev quadrature.
 *
 * The tabulated values are calculated at log-spaced energy
 * and impact factor values according to the corteo method.
 *
 *
 * @ingroup XS
 */
struct xs_corteo4bit : public screening_function< ScreeningZBL >
{
    typedef corteo4bit corteo_idx_t;

    const static int rows = corteo4bit::e_index::dim;
    const static int cols = corteo4bit::s_index::dim;

    /**
     * @brief sin2Thetaby2
     * @param e
     * @param s
     * @return
     */
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

/**
 * @brief The xs_corteo6bit class
 * @ingroup XS
 */
struct xs_corteo6bit : public screening_function< ScreeningZBL >
{
    typedef corteo6bit corteo_idx_t;

    const static int rows = corteo6bit::e_index::dim;
    const static int cols = corteo6bit::s_index::dim;

    static double sin2Thetaby2(const double& e, const double& s)
    {
        const float* p = corteo4bitdata();
        int i = corteo6bit::table_index(e,s);
        return p[i];
    }
    static float sin2Thetaby2(int ie, int is)
    {
        const float* p = corteo4bitdata();
        return p[ie*cols + is];
    }
    static const float* data() { return corteo4bitdata(); }
};

// center-of-mass to lab parameters
struct cm_pars {
    float screening_length;    /* screening length [nm] */
    float mass_ratio;          /* M1/M2 */
    float sqrt_mass_ratio;     /* we will need this occasionally */
    float kfactor_m;           /* 4 M1 M2 / (M1 + M2)^2 */
    float red_E_conv;          /* reduced energy conversion factor */

    template<class _XScm>
    void init(float Z1, float M1, float Z2, float M2)
    {
        /* Adapted from the corteo code */
        screening_length     = _XScm::screeningLength(Z1,Z2);
        mass_ratio           = M1/M2;
        sqrt_mass_ratio      = std::sqrt(mass_ratio);
        kfactor_m            = 4*mass_ratio / ((mass_ratio+1) * (mass_ratio+1));
        red_E_conv           = screening_length / ((mass_ratio+1) * Z1 * Z2 * E2);
    }
};

/**
 * @brief The xs_lab class template describes a cross-section in lab system
 *
 * After initializing the class with XSlab::init for a given combination
 * of projectile and target, the function XSlab::scatter can be used to
 * calculate scattering quantities.
 *
 * @tparam
 *   xs_cm is the class of the center-of-mass cross-section
 *
 * @ingroup XS
 */
template<class xs_cm>
class xs_lab
{
    cm_pars P_;
public:
    xs_lab() : P_()
    {}
    xs_lab(const xs_lab& x) : P_(x.P_)
    {}
    float sqrtMassRatio() const { return P_.sqrt_mass_ratio; }
    /**
     * @brief Initialize the cross-section for a specific projectile-target combination
     * @param Z1 the projectile atomic number
     * @param M1 the projectile mass
     * @param Z2 the target atomic number
     * @param M2 the target mass
     */
    void init(float Z1, float M1, float Z2, float M2) {
        P_.init<xs_cm>(Z1,M1,Z2,M2);
    }
    /**
     * @brief Calculate scattering angle and target recoil energy.
     *
     * Given the initial energy E and impact parameter S of an incoming
     * projectile, scatter calculates the target atom recoil energy and the projectile
     * scattering angle.
     *
     * All quantities refer to the lab system.
     *
     * @param E is the initial projectile energy [eV]
     * @param S is the impact factor [nm]
     * @param recoil_erg is the target recoil energy [eV]
     * @param sintheta the sin of the scattering angle in the lab system
     * @param costheta the cos of the scattering angle in the lab system
     */
    void scatter(float E, float S,
                 float &recoil_erg, float &sintheta, float &costheta) const
    {
        float e = E*P_.red_E_conv;
        float sin2thetaby2 = xs_cm::sin2Thetaby2(e, S/P_.screening_length);
        recoil_erg = e*sin2thetaby2;
        /* convert scattering angle to lab frame of reference: */
        costheta = 1.f - 2*sin2thetaby2;
        sintheta = std::sqrt(1.f-costheta*costheta);
        float theta=atan( sintheta/(costheta+P_.mass_ratio) );
        sintheta=sin(theta);
        costheta=cos(theta);
    }
    /**
     * @brief Calculate the impact parameter with given initial projectile energy and target recoil energy
     * @param E the projectile initial energy [eV]
     * @param T the target recoil energy [eV]
     * @return the corresponding impact factor [nm]
     */
    float impactPar(float E, float T)
    {
        double thetaCM = T/E/P_.kfactor_m;
        thetaCM = 2.*std::asin(std::sqrt(thetaCM));
        return xs_cm::findS(E*P_.red_E_conv,thetaCM)*P_.screening_length;
    }
    /**
     * @brief Differential cross-section \f$ d\sigma(E,T)/dT \f$
     * @param E projectile energy [eV]
     * @param T recoil energy [eV]
     * @return the cross-section [nm^2/eV]
     */
    float crossSection(float E, float T)
    {
        double thetaCM = T/E/P_.kfactor_m;
        thetaCM = 2.*std::asin(std::sqrt(thetaCM));
        return xs_cm::crossSection(E*P_.red_E_conv,thetaCM)*4*M_PI*P_.screening_length*P_.screening_length/E/P_.kfactor_m;
    }
};


// implementation of XSlab for corteo-tabulated center-of-mass XS
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
        // bilinear interp
        const float* p = _XS::data();
        recoil_erg = e*P_.kfactor_m;
        Eigen::Vector4f sinth, costh, sin2thetaby2;
        e *= P_.red_E_conv;
        s /= P_.screening_length;
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
        return xs_quad<ScreeningZBL>::findS(E*P_.red_E_conv,thetaCM)*P_.screening_length;

    }
    float crossSection(float E, float T)
    {
        double thetaCM = T/E/P_.kfactor_m;
        thetaCM = 2.*std::asin(std::sqrt(thetaCM));
        return xs_quad<ScreeningZBL>::crossSection(E*P_.red_E_conv,thetaCM)*4*M_PI*P_.screening_length*P_.screening_length/E/P_.kfactor_m;
    }
};

/**
 * @brief Explicit specialization of xs_lab class for 4-bit corteo-tabulated
 * center-of-mass cross-section
 *
 * @ingroup XS
 */
template<>
class xs_lab<xs_corteo4bit> : public xs_corteo_impl_< xs_corteo4bit >
{};

/**
 * @brief Explicit specialization of xs_lab class for 6-bit corteo-tabulated
 * center-of-mass cross-section
 *
 * @ingroup XS
 */
template<>
class xs_lab<xs_corteo6bit> : public xs_corteo_impl_< xs_corteo6bit >
{};

#endif
