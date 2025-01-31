clear


Z1 = 26; M1 = 55.845;
Z2 = 26; M2 = 55.845;
N=84.822; # Fe, at/nm^3

A = M1/M2;
gam = 4*A/(1+A)^2;


E2 = 1.43996445; # e^2/4pi eps0 [eV-nm]
ab = 0.05291772108; # Bohr radius [nm]
a = 0.88534*ab/(Z1^0.23 + Z2^0.23);
c = a/(1+A)/Z1/Z2/E2;

# low E stopping E < 1 eV
# fit S(E)=Sn+Se to a power law S = S1*E^n
# l(E) = int_0^E { dx / S(E) } = E^(1-n)/S1/(1-n)

E=logspace(-3,0,31);

Sn = screened_coulomb_sn(E*c)*pi*a^2*gam/c*N;
[e,d] = dedx(Z1,Z2);
Se = (E/e(1)).^0.5*d(1)*N*0.1;
St = Se+Sn;
p = polyfit(log(E),log(St),1);
S1 = exp(p(2));
n = p(1);

#loglog(E, [St; S1*E.^n])
#return

# 0.1 < E < 100eV

E=logspace(-1,2,31);

Sn = screened_coulomb_sn(E*c)*pi*a^2*gam/c*N;

St = Sn + (E/e(1)).^0.5*d(1)*N*0.1;

figure 1
clf
plot(E, cumtrapz(E,1./St) + 1/S1/(1-n)*E(1).^(1-n))
xlabel('Ion energy E (eV)')
ylabel('Range (nm)')
title('Low energy Fe recoil ions range in Fe')

