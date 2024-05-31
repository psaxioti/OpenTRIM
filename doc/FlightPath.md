# Flight Path Selection

The Monte-Carlo simulation of particle transport typically proceeds in the following manner:

- Sample the distance to the next collision
- Propagate the particle to the collision site
- Select collision partner
- Sample energy, angle, etc. of the particle after collision
- Repeat 

The distance to the next collision is typically sampled from the Poisson distribution, $p(x) = e^{-N\sigma_0 x}=e^{-x/\ell}$, where $N$ is the atomic density of scattering centers, $\sigma_0$ the total cross-section and $\ell = (N\sigma_0)^{-1}$ denotes the mean free path (mfp).

There is a difficulty in employing this type of flight path selection for the Coulomb and screened Coulomb interactions because the total cross-section diverges. This can be circumvented by setting a cutoff based on scattering angle or recoil energy or some other scattering parameter. Scattering events outside the cutoff are ignored. The cutoff value is selected so that the result of the calculation is not affected significantly.

Mendenhall & Weller 2005 set a such a lower cutoff for the recoil energy, $T_c$, which corresponds to the lowest recoil energy of interest in the problem under study. M&H suggest a value in the range of 1 - 10 eV for ion penetration in solids. 

$T_c$ corresponds to a lower bound $\theta_c$ for the center-of-mass scattering angle, which can be obtained from 
$$
T_c = T_m \sin^2(\theta_c/2)
$$
For a given reduced energy $\epsilon$ of the particle, we can find the maximum impact parameter, $p_{max}$, corresponding to $\theta_c$ and thus define the effective total cross-section
$$
\sigma_0 = \pi\, p_{max}^2(\epsilon) = 
\int_0^{p_{max}}{2\pi\, p\, dp} =
\int_{\theta_c}^{\pi}{d\sigma(\epsilon,\theta)}
$$

## Mendenhall-Weller Algorithm

>- For a given particle energy $\epsilon$ we know the values of $p_{max}$, $\sigma_0$ and $\ell$
>- Take a random number sample $u \in (0,1)$
>- If $p_{max} < L$, where $L$ is half the interatomic distance, 
>    - $p = p_{max} \sqrt{-\log u}$
>    - Propagate particle by $\ell$
>    - If $p>p_{max}$ ignore the scattering event 
>- if $p_{max} \geq L$ then 
>    - $p = \sqrt{u /(\pi N L)}$
>    - Propagate particle by $L$ 

This is essentially the same as the algorithm used by SRIM.

The 2nd part of the algorithm, for $p_{max} \geq L$, ensures impact parameters larger than the interatomic distance do not occur. This is accepted in both SRIM and MH2005. MH justify it as follows:

> ... impact parameters
larger than half the lattice spacing do not occur,
since then one is closer to the adjacent atom.

while Biersack1980 and the SRIM manual (Ch.7) write 

> This
procedure maintains the atomic density in the target
without correlating the lateral positions of successive
target atoms (neglection of lattice structure).

We believe that these propositions are not justified and one can assume a mfp smaller than the interatomic distance without any problems.

## Standard algorithm

We propose to employ a standard procedure for path selection similar to neutron or photon transport simulations:

> - For a given particle energy $\epsilon$ we have $p_{max}$, $\sigma_0=\pi p_{max}^2$ and $\ell = (N\sigma_0)^{-1}$
> - Take 2 random samples $u_{1,2} \in (0,1)$
> - Path to the next collision $x = -\ell\,\log u_1$ [Poisson distr.]
> - $p = p_{max}\sqrt{u_2}$
> - Propagate particle by $x$

Although here we need 2 random numbers instead of 1 in the MH/SRIM algorithm, the penalty is not so high since no outcome is discarded. 

An advantage of this method is that the average mfp of the simulated ions is exactly $\ell$, which in the MH algorithm is not true. This permits to set a particular value for the mfp, which can be useful, e.g., when simulating thin targets (setting $\ell \ll d$ ensures scattering in the target)

## Criteria for setting $p_{max}$, $\sigma_0$ and $\ell$

A number of different criteria can be used in parallel:

- **Recoil energy lower cutoff** \
  As explained above, a lower cutoff $T_c$ in the recoil energy results in a maximum impact parameter $p_{max} = p(\epsilon,T_c)$ and $\ell = 1/\sqrt{\pi N p_{max}^2}$

- **Upper bound for electronic stopping** \
  $\Delta E / E = (dE/dx)\,\ell/E< \delta$ or $\ell < \ell_{max} = \delta \cdot E/(dE/dx)$. Typically a $\delta$ value of 1 - 5% is employed

- **Upper bound for the mfp** \
  $\ell < \ell_{max}$, where $\ell_{max}$ is set for geometric reasons

A table of these values as a function of incident energy is calculated before starting the main simulation.


