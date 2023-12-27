#include "xs.h"
#include "corteo.h"

// Following functions used to compute scattering angle following Gauss-Mehler quadrature
template<class PHI>
double xs_quad<PHI>::H(double u, double x0, double epsilon, double s) {
	double x0u = x0/u;
	double eps1 = 1./epsilon;
    double ss = s*s;
    return std::sqrt((1-u*u)/(1-PHI::screeningFunction(x0u)/x0u*eps1-ss/x0u/x0u));
}

// the root of this function (vs x0) is the minimal 
// approach distance at energy epsilon and impact parameter s
template<class PHI>
double xs_quad<PHI>::funcX0(double x0, double epsilon, double s) {
	double sx0 = s/x0;
    return 1.-PHI::screeningFunction(x0)/(x0*epsilon)-sx0*sx0;
}

// use bisection (Newton) method to find the minimal approach distance x0
// (return THETAERR on error)
template<class PHI>
double xs_quad<PHI>::findX0(double epsilon, double s) {
	double funcX0val1, funcX0val2, funcX0valm;

	double temp = 1.0/(2.0*epsilon);
	double x2 = temp+sqrt(temp+s*s); // inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15 

	double x1 = x2/10.;
	funcX0val1 = funcX0(x1, epsilon, s);  // should be always negative
	funcX0val2 = funcX0(x2, epsilon, s);  // should be always positive (starting ~1.0)

	if(funcX0val1>=0.) {
		// initial guess for x1 too optimistic, start with a safe value (much longer)
        x1 = 1.e-8;
		funcX0val1 = funcX0(x1, epsilon, s);  // should be always negative
	}

    if(funcX0val1>=0. || funcX0val2<=0.)
        return std::nan("");    // error, values should be on each side of 0

    double xm = 0.5*(x1+x2);
    double fm = funcX0(xm, epsilon, s);
	do {
        if(fm<0.) {
            // funcX0val1 = funcX0valm;
            x1 = xm;
		} else {
            // funcX0val2 = funcX0valm;
            x2 = xm;
		}
        xm = 0.5*(x1+x2);
        fm = funcX0(xm, epsilon, s);
    } while (fabs(fm)>1e-14);
	// 1.-XXX wont be more precise than 2^(-52)=2e-16
	// but using 1e-10 saves time and is enough precise for float conversion later on

    return xm;
}

// use bisection method to find the reduced 
// impact parameter s given epsilon and thetaCM 
// (return THETAERR on error)
template<class PHI>
double xs_quad<PHI>::finds(double epsilon, double thetaCM) const {
	double funcX0val1, funcX0val2, funcX0valm;

	// inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
	double gamma = (M_PI-thetaCM)/M_PI;
    double x0 = findX0((1.0-gamma*gamma)*epsilon, 1.e-8);
	double x1 = 0.7*gamma*x0;
	double x2 = 1.0/(2.0*epsilon*tan(thetaCM/2.0)); 

    funcX0val1 = thetaCM-theta(epsilon, x1);  // should be always negative
    funcX0val2 = thetaCM-theta(epsilon, x2);  // should be always positive (starting ~1.0)
    funcX0valm = thetaCM-theta(epsilon, (x1+x2)/2.);

    if(funcX0val1>=0. || funcX0val2<=0.)
        return std::nan("");    // error, values should be on each side of 0

	do {
		if(funcX0valm<0.) {
			funcX0val1 = funcX0valm;
			x1 = (x1+x2)/2.;
		} else {
			funcX0val2 = funcX0valm;
			x2 = (x1+x2)/2.;
		}
        funcX0valm = thetaCM-theta(epsilon, (x1+x2)/2.);
	} while (fabs(funcX0valm)>1e-5);

	return (x1+x2)/2.;
}

// scattering angle calculated using Gauss-Mehler quadrature
// (return THETAERR on error)
template<class PHI>
double xs_quad<PHI>::theta(double epsilon, double s) const
{
    double x0 = findX0(epsilon, s);
    if(std::isnan(x0)) {
        return x0;
	}

    double sum(0.);
    for(unsigned int j=0; j<NSUM/2; j++)
        sum += H(cos((2.*j-1.)*M_PI/(2.*NSUM)), x0, epsilon, s);

    return M_PI*(1.-2.*s/NSUM/x0*sum);
}

/* returns the cross section in the center of mass frame considering a screened potential */
template<class PHI>
double xs_quad<PHI>::crossSection(double E, unsigned int Z1, unsigned int Z2, double massRatio, double thetaCM) const
{
    double ds, dsdTheta, s, epsilon;
    double screenLength = PHI::screeningLength(Z1,Z2);

	// reduced energy according to screening length
	epsilon = E*screenLength/(Z1*Z2*E2);

	// find corresponfding reduced impact parameter 
    s = finds(epsilon, thetaCM);
    if(std::isnan(s)) return s;

	// ds/dTheta using five-point stencil 
	ds = s*0.01;
    dsdTheta = (12.0*ds)/(-theta(epsilon,s+2.0*ds)
                          +8.0*theta(epsilon,s+ds)
                          -8.0*theta(epsilon,s-ds)
                          +theta(epsilon,s-2.0*ds));

	// return "non-reduced" center-of-mass cross section 
	return screenLength*screenLength*s/sin(thetaCM)*fabs(dsdTheta);

}

// explicit instantiation of all variants
template class xs_quad<screening_none>;
template class xs_quad<screening_lj>;
template class xs_quad<screening_univ>;

double xs_zbl::MAGIC(double epsilon, double B){
    /* B: reduced impact par
     epsilon: reduced center of mass energy
     returns cos(theta/2) of the scattering event */

    double cost2;  /* cos(theta/2)*/
    double RoC,Delta,R,RR,A,G,alpha,beta,gamma,V,V1,FR,FR1,Q;
    double SQE;
//    double C[6];

//    C[1]=0.99229;  /* TRIM 85:  */
//    C[2]=0.011615;
//    C[3]=0.0071222;
//    C[4]=14.813;
//    C[5]=9.3066;

    /* TRIM 85:  */
    static const double C[] = {0.,
                               0.99229,
                               0.011615,
                               0.0071222,
                               14.813,
                               9.3066};

    /* Initial guess for R: */
    R=B;
    RR=-2.7*log(epsilon*B);
    if(RR>=B){
        /*   if(RR<B) calc potential; */
        RR=-2.7*log(epsilon*RR);
        if(RR>=B){
            R=RR;
        }
    }
    /* TRIM85: 330 */
    do{
        /* Calculate potential and its derivative */
        V=xs_zbl::ZBL_and_deri(R,&V1);
        FR  = B * B / R + V*R/epsilon - R;
        FR1 = - B * B / (R * R) + (V+V1*R)/epsilon - 1.0;
        Q   = FR/FR1;
        R   = R-Q;
    } while(fabs(Q/R)>0.001);

    RoC = -2.0 * (epsilon-V)/V1;
    SQE = sqrt(epsilon);

    alpha = 1+ C[1]/SQE;
    beta  = (C[2]+SQE) / (C[3]+SQE);           /* TRIM85: CC */
    gamma = (C[4]+epsilon)/(C[5]+epsilon);
    A     = 2*alpha*epsilon*pow(B,beta);
    G     = gamma / ( sqrt((1.0+A*A))-A  );    /* TRIM85: 1/FF */
    Delta = A * (R-B)/(1+G);

    cost2=(B+RoC+Delta)/(R+RoC);
    return cost2;
}

double xs_zbl::ZBL_and_deri(double R, double* Vprime){
    /* return ZBL potential, and via the pointer Vprime its derivative */
    /* Values are taken from ZBL85 */

    double EX1,EX2,EX3,EX4,V;
    /* corteo:    return 0.1818*exp(-3.*x)+0.5099*exp(-0.9423*x)+0.2802*exp(-0.4028*x)+0.02817*exp(-0.2016*x); */
    /* EX1=0.1818   * exp( -3.0  * R);
     EX2=0.5099   * exp( -0.9423 * R);
     EX3=0.2802   * exp( -0.4028  * R);
     EX4=0.02817  * exp( -0.2016  * R);
     V=(EX1+EX2+EX3+EX4)/R;
     *Vprime = -(V+3.0*EX1+0.9423*EX2 + 0.4028*EX3 + 0.2016*EX4)/R;
     return V;*/

    auto &C =  screening_univ::C;
    auto &A =  screening_univ::A;

    /* TRIM85: */
    // EX1=0.18175  * exp( -3.1998  * R);
    /*  if(R>=7){EX1=0.0;}*/ /* According to TRIM95 */
    // EX2=0.50986  * exp( -0.94229 * R);
    // EX3=0.28022  * exp( -0.4029  * R);
    // EX4=0.028171 * exp( -0.20162 * R);

    EX1 = C[0]*exp(A[0]*R);
    EX2 = C[1]*exp(A[1]*R);
    EX3 = C[2]*exp(A[2]*R);
    EX4 = C[3]*exp(A[3]*R);

    V=(EX1+EX2+EX3+EX4)/R;

    if (Vprime)
        *Vprime = -(V - A[0]*EX1 - A[1]*EX2 - A[2]*EX3 - A[3]*EX4)/R;

    return V;
}

void scattering::init(float Z1, float M1, float Z2, float M2)
{
    /* Creates the matrix with scattering results and stores it to the structure pointed to by ScatMatrix */

    /* Adapted from the corteo code */

    screening_length     = xs.screeningLength(Z1,Z2); /* This is in nm! */
    inv_screening_length = 1.0f/screening_length;
    mass_ratio           = M1/M2;
    sqrt_mass_ratio      = std::sqrt(mass_ratio);
    kfactor_m            = 4*mass_ratio / ((mass_ratio+1) * (mass_ratio+1)); /* k factor without the angle part */
    red_E_conv           = screening_length / ((mass_ratio+1) * Z1 * Z2 * E2);

    /* The maximum reduced impact parameter depends on density and flight length. This
    means that we cannot calculate it as a function of projectile and target element only.
    Thus we need to calculate it somewhere else (as a material property of target). This
    is different from corteo, where the scattering matrices are stored for each target
    layer with known density etc. */

    /* ScatMatrix->max_red_im_par   = 1.0f / ( ScatMatrix->screening_length * sqrtdf( PI * a )  ) */
    /* ion->smax[ilayer][ielem] = 1.0f/(ion->a[ilayer][ielem]*sqrtdf(PI*layerDensity[ilayer])*sqrtMeanFreePath[ilayer]);  */
    /* float max_red_im_par;       Maximum reduced impact parameter */

}

void scattering::scatter(float erg, float P,
                         float& recoil_erg, float& sintheta, float& costheta)
{
    float sin2thetaby2 = xs.sin2Thetaby2(erg*red_E_conv, P*inv_screening_length);
    costheta = 1.f - 2*sin2thetaby2;
    sintheta = std::sqrt(1.f-costheta*costheta);
    /* now conversion to lab frame of reference: */
    float theta=atan( sintheta/(costheta+mass_ratio) );
    /*if(isnan(theta)){theta=0;}*/
    sintheta=sin(theta);
    costheta=cos(theta);
    recoil_erg = erg*kfactor_m*sin2thetaby2;
}

scattering_corteo::scattering_corteo(const xs_base& axs, const corteo& c) :
    scattering(axs), corteo_(c),
    cosTable(c.size()),sinTable(c.size())
{}

/*************** Functions adapted from corteo.c ***********************/

/* fill the table of cos and sin of the scattering angle in the lab frame given a mass ratio */
void scattering_corteo::fillCosSinTable()
{
    double sin2thetaby2, costheta, costhetaLab, sinthetaLab;
    double mr = mass_ratio;

    /* compute scattering angle components */
    const float* corteo_tbl = corteo_.data();
    for(int i=0; i<corteo_.size(); i++) {
        sin2thetaby2 = corteo_tbl[i];
        costheta = 1.-2.*sin2thetaby2;

        if(costheta==-1.0 && mr==1.0) {
            costhetaLab = 0.0;  /* peculiar case of head-on collision of identical masses */
            sinthetaLab = 1.0;
        } else {
            costhetaLab = (costheta+mr)/sqrt(1.+2.*mr*costheta+mr*mr);
            sinthetaLab = sqrt(1.-costhetaLab*costhetaLab);
        }
        cosTable[i] = costhetaLab;
        sinTable[i] = sinthetaLab;

        /* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
        /* In some rare cases, when cos=1, then sin becomes "Not a Number". To prevent this, I will set the sine to 0 in those cases. */
        if( std::isnan(sinTable[i]) ) {
            sinTable[i]=0.f;
            cosTable[i]=1.f;
        }
    }
}

void scattering_corteo::init(float Z1, float M1, float Z2, float M2)
{
    scattering::init(Z1,M1,Z2,M2);
    fillCosSinTable();
}

void scattering_corteo::scatter(float erg, float P,
                         float& recoil_erg, float& sintheta, float& costheta)
{
    int i = corteo_.table_index(erg*red_E_conv, P*inv_screening_length);
    float sin2thetaby2 = corteo_.data()[i];
    recoil_erg = erg*kfactor_m*sin2thetaby2;
    sintheta=sinTable[i];
    costheta=cosTable[i];
}
