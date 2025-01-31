clear

% The phase space of (e,s) used in iradina
% e : reduced energy
% s : reduced impact parameter
e = logspace(-6,6,51);
s = logspace(-8,2,31);
[E,S] = meshgrid(e,s);

% calc scattering angle
Th = screened_coulomb_theta(E,S);

% Plot scattering angle
figure 1
clf
loglog(e,Th,'.-')
xlabel('\epsilon')
ylabel('\theta')
title('Scattering angle \theta(\epsilon,s)')

% get back S from (E,Th)
S1 = screened_coulomb_ip(E,Th);
figure 2
clf
% 100/e^(1/6) is the limit where impulse approx is used
loglog(e,S1,'.-',e,100./e.^(1/6))
xlabel('\epsilon')
ylabel("s'")
title('Get back s'' = s(\epsilon,\theta)')

figure 3
clf
loglog(e,abs((S1-S)./S),'.-')
xlabel('\epsilon')
ylabel("\delta")
title("Error \delta = (s'- s)/s")

% calc cross-section
XS = screened_coulomb_xs(E,Th);
figure 4
clf
loglog(e,XS,'.-')
xlabel('\epsilon')
title('d\sigma/d\Omega (\epsilon,\theta)')

