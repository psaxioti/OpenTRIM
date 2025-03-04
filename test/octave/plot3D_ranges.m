clear
clf

%% only nb=6,7 of interest
nb = 6;
%monolayer
%pathmnl = ['../srim/mnl_3d/b' num2str(nb) '/'];
%S = srim2mat(pathmnl);

% FC mode
path = ['../srim/b' num2str(nb) '/'];
S = srim2mat(path);

%   Dis:   displacement 3d Displacements/(Angstrom-Ion)
%   EI3:   ionization energy eV/(Angstrom-Ion)
%   EP3:   phonon energy eV /(Angstrom-Ion)
%   Ri3:   rande distribution ions
%   Rr3:   range distribution recoils
%   RC3:   replacement collision Replacements/(Angstrom-Ion)
%   V3r:   target atom vacancies Vacancies/(Angstrom-Ion)
%

# Size of the cell
cellsize_all = [12 6 12 6 3 500 100];
cellsize = cellsize_all(nb) ;

x = S.x/10;
%displacement
figure 1

subplot(1,2,1)
surf(x,x,S.Dis*cellsize*10)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# displacements/ion');
view(40,30)
title('Displacements/ion')
subplot(1,2,2)
contour(x,x,S.Dis*cellsize*10,100)
colorbar


% Ionization
figure 2
subplot(1,2,1)
surf(x,x,S.EI3*cellsize*10)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# ionization energy (eV)/ion');
view(40,30)
subplot(1,2,2)
contour(x,x,S.EI3*cellsize*10,100)
colorbar
title('Phonon energy/ion')

% !!! Ion distribution/implantation
figure 4
subplot(1,2,1)
surf(x,x,S.Ri3.*cellsize*10e-8)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# ions');
colormap('jet')
view(40,30)
subplot(1,2,2)
contour(x,x,S.Ri3.*cellsize*10e-8,100)
colorbar
title('Ion distribution')

%Recoil distribution
figure 5
subplot(1,2,1)
surf(x,x,S.Rr3*cellsize*10e-8)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# recoils');
view(40,30)
subplot(1,2,2)
contour(x,x,S.Rr3*cellsize*10e-8,100)
colorbar
title('Recoil distribution')

% Replacement colisions distribution
figure 6
subplot(1,2,1)
surf(x,x,S.RC3*cellsize*10)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# replacements/ion');
view(40,30)
subplot(1,2,2)
contour(x,x,S.RC3*cellsize*10,100)
colorbar
title('Replacements/ion')

%Vacancy distribution
figure 7
subplot(1,2,1)
surf(x,x,S.V3r*cellsize*10)
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('# ions');
view(40,30)
subplot(1,2,2)
contour(x,x,S.V3r*cellsize*10,100)
colorbar
title('Target vacancies/ion')

