# Flight Path Selection {#flightpath}

The Monte-Carlo simulation of particle transport typically proceeds in the following manner:

- Sample the distance to the next collision
- Propagate the particle to the collision site
- Select collision partner
- Sample energy, angle, etc. of the particle after collision
- Repeat 

The distance to the next collision is typically sampled from the Poisson distribution, \f$p(x) = e^{-N\sigma_0 x}=e^{-x/\ell}\f$, where \f$N\f$ is the atomic density of scattering centers, \f$\sigma_0\f$ the total cross-section and \f$\ell = (N\sigma_0)^{-1}\f$ denotes the mean free path (mfp).

There is a difficulty in employing this type of flight path selection for the Coulomb and screened Coulomb interactions because the total cross-section diverges. This can be circumvented by setting a lower cutoff, either w.r.t. the scattering angle or the recoil energy or some other scattering parameter. Scattering events below the cutoff are ignored. The cutoff value is selected so that the result of the calculation is not affected significantly.

Mendenhall & Weller 2005 set such a lower cutoff for the recoil energy, \f$T_c\f$, which corresponds to the lowest recoil energy of interest in the problem under study. M&W suggest a value in the range of 1 - 10 eV for ion penetration in solids. 

\f$T_c\f$ corresponds to a lower bound \f$\theta_c\f$ for the center-of-mass scattering angle, which can be obtained from 
$$
T_c = T_m \sin^2(\theta_c/2)
$$
For a given reduced energy \f$\epsilon\f$ of the particle, we can find the maximum impact parameter, \f$p_{max}\f$, corresponding to \f$\theta_c\f$ and thus define the effective total cross-section
$$
\sigma_0 = \pi\, p_{max}^2(\epsilon) = 
\int_0^{p_{max}}{2\pi\, p\, dp} =
\int_{\theta_c}^{\pi}{d\sigma(\epsilon,\theta)}
$$

## Mendenhall-Weller Algorithm

> - For a given particle energy \f$\epsilon\f$ find the values of \f$p_{max}\f$, \f$\sigma_0\f$ and \f$\ell\f$ corresponding to the recoil energy cutoff \f$T_c\f$
> - Take a random number sample \f$u \in (0,1)\f$
> - If \f$p_{max} < L\f$, where \f$L\f$ is half the interatomic distance, 
>    - \f$p = p_{max} \sqrt{-\log u}\f$
>    - If \f$p>p_{max}\f$ ignore the scattering event 
>    - Propagate particle by \f$\ell\f$
> - if \f$p_{max} \geq L\f$ then 
>    - \f$p = \sqrt{u /(\pi N L)}\f$
>    - Propagate particle by \f$L\f$ 

This is essentially the same as the algorithm used by SRIM.

The 2nd part of the algorithm, for \f$p_{max} \geq L\f$, ensures that impact parameters larger than the interatomic distance do not occur. This is accepted in both SRIM and M&W2005. M&W justify it as follows:

> ... impact parameters larger than half the lattice spacing do not occur, since then one is closer to the adjacent atom.

while Biersack1980 and the SRIM manual (Ch.7) write:

> This
> procedure maintains the atomic density in the target
> without correlating the lateral positions of successive
> target atoms (neglection of lattice structure).

We believe that these propositions are not justified and one can assume a mfp smaller than the interatomic distance without any problems.

## Standard algorithm

We propose to employ a more "standard" procedure for path selection similar to what is done for neutron or photon transport simulations:

> - For a given particle energy \f$\epsilon\f$ we have \f$p_{max}\f$, \f$\sigma_0=\pi p_{max}^2\f$ and \f$\ell = (N\sigma_0)^{-1}\f$
> - Take 2 random samples \f$u_{1,2} \in (0,1)\f$
> - Path to the next collision \f$x = -\ell\,\log u_1\f$ [Poisson distr.]
> - \f$p = p_{max}\sqrt{u_2}\f$
> - Propagate particle by \f$x\f$

Although here we need 2 random numbers instead of 1 in the MW/SRIM algorithm, the penalty is not so high since all outcomes are valid (in M&W algorithm some events are discarted). 

An advantage of this method is that the average mfp of the simulated ions is exactly \f$\ell\f$, which in the MW algorithm is not true. This permits us to preset a value for the mfp. This can be useful, e.g., when simulating thin targets (setting \f$\ell \ll d\f$ ensures scattering in the target)

## Criteria for setting \f$p_{max}\f$, \f$\sigma_0\f$ and \f$\ell\f$

A number of different criteria can be used in parallel:

- **Recoil energy lower cutoff** \n 
  As explained above, a lower cutoff \f$T_c\f$ in the recoil energy results in a maximum impact parameter \f$p_{max} = p(\epsilon,T_c)\f$ and \f$\ell = 1/\sqrt{\pi N p_{max}^2}\f$

- **Upper bound for electronic stopping** \n 
  \f$\Delta E / E = (dE/dx)\,\ell/E< \delta\f$ or \f$\ell < \ell_{max} = \delta \cdot E/(dE/dx)\f$. Typically a \f$\delta\f$ value of 1 - 5% is employed

- **Upper bound for the mfp** \n 
  \f$\ell < \ell_{max}\f$, where \f$\ell_{max}\f$ is set for geometric reasons

A table of these values as a function of incident energy is calculated before starting the main simulation.


