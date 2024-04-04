clear

% The phase space of (e,s) used in iradina
% e : reduced energy
% s : reduced impact parameter
e = logspace(-6,6,51);
s = logspace(-8,2,31);
[E,S] = meshgrid(e,s);

% calc scattering angle
Th = screened_coulomb_theta(E,S);
figure 1
clf
loglog(e,Th,'.-')

% get back S from (E,Th)
S1 = screened_coulomb_ip(E,Th);
figure 2
clf
loglog(e,S1,'.-',e,100./e.^(1/6))
figure 3
clf
loglog(e,abs((S1-S)./S),'.-')

% calc cross-section 
XS = screened_coulomb_xs(E,Th);
figure 4
clf
loglog(e,XS,'.-')

