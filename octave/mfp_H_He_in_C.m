clear

E0 = logspace(0,6,31); % eV
Z1 = 1;  M1 = 1.00784; % H
Z1 = 2;  M1 = 4.003; % He
Z2 = 6; M2 = 12;
gam = 4*M1*M2/(M1+M2)^2;
Tm = gam*E0;
Tc_abs = 10; % eV
Tc_rel = 0.999;
Tc = Tc_abs*ones(size(Tm)); % eV
i = find( Tc > Tc_rel*Tm );
Tc(i) = Tc_rel*Tm(i);

mfpmax = 10; % 44;


ab = 0.05291772108; % Bohr radius nm
e2 = 1.43996445; % e^2 / 4 pi eps0 = e^2 c^2 in [eV][nm]
a = 0.88534*ab/(Z1^0.23+Z2^0.23);
C = a / (1+M1/M2)/Z1/Z2/e2;
e = E0 * C;

n = 6.02214e2*2.27/M2; % at. density at/nm^3

% electronic stopping
[Ee, Se] = dedx(Z1,Z2);
Se = Se*n*0.1; % eV/nm
i = find(Ee>max(E0));
Ee(i) = [];
Se(i) = [];
if Ee(1) > E0(1),
  Se = [Se(1)*sqrt(E0(1)/Ee(1)) Se];
  Ee = [E0(1) Ee];
end
% nuclear stopping
Sn = screened_coulomb_sn(e)*pi*a^2*gam/C*n;

% atomic distance
Lat = (4*pi*n/3)^(-1/3);
Pat = 1/sqrt(pi*n*Lat);
that=screened_coulomb_theta(e,Pat/a);
Tat = sin(that/2).^2.*Tm;
that = 2*asin(sqrt(Tat./Tm));

% mfp for T>Tc
thc = 2*asin(sqrt(Tc./Tm));
Pc = a*screened_coulomb_ip(e,thc);
Lc = 1/n/pi./Pc.^2;
Snc = screened_coulomb_sn(e,thc)*pi*a^2*gam/C*n;

% path for Se*dx/E < delta
delta = 0.05;
Le = delta*E0./exp(spline(log(Ee),log(Se),log(E0)));
%Le = delta*Ee./Se;
Pe = 1./sqrt(pi*n*Le);
the = screened_coulomb_theta(e,Pe/a);
Te = sin(the/2).^2.*Tm;

% final
L0 = Lc;
i = find(L0 > Le);
L0(i) = Le(i);
i = find(L0 > mfpmax);
L0(i) = mfpmax;
i = find(L0 < Lat);
%L0(i) = Lat;
P0 = 1./sqrt(pi*n*L0);
th0 = screened_coulomb_theta(e,P0/a);
T0 = sin(th0/2).^2.*Tm;
th0 = 2*asin(sqrt(T0./Tm));
Sn0 = screened_coulomb_sn(e,th0)*pi*a^2*gam/C*n;

figure 1
clf
loglog(E0,Lat*ones(size(E0)),'.-',E0,Lc,E0,Le,E0,L0)
ylabel('mfp (nm)')
xlabel('E (eV)')
legend('ML','Tc 1eV','DEe/E < 5%','Total')

##figure 2
##clf
##loglog(E0,Pat*ones(size(E0)),E0,Pc,E0,Pe,E0,P0)
##ylabel('impact parameter (nm)')
##xlabel('E (eV)')
##legend('atomic distance','T cutoff 1eV','DEe/E < 1%','Total')

figure 2
clf
loglog(E0,Tat,E0,Tc.*ones(size(E0)),E0,Te,E0,T0)
ylabel('recoil energy cutoff (eV)')
xlabel('E (eV)')
legend('ML','Tc 1eV','DEe/E < 5%','Total')

figure 3
clf
loglog(E0,that,E0,thc,E0,the,E0,th0)
ylabel('\theta (rad)')
xlabel('E (eV)')
legend('ML','Tc 1eV','DEe/E < 5%','Total')

##figure 4
##clf
##loglog(Ee,Se,E0,Sn,E0,Snc,E0,Sn0)
##ylabel('Se, Sn (eV/nm)')
##xlabel('E (eV)')
##legend('e','n','n T cutoff 1eV','Total')








