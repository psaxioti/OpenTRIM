clear
clf

## !!!!! EDIT REGION - EDITABLE VALUES  !!!!
# Benchmark Number
nb = 3;

# Number of ions run
### nb = 1,4,5,6,7 Nions=20000
## nb = 2 Nions=20
## nb = 3 Nions=8000

Nions = 8000;

# Size of the cell
### nb = 1,3 cellsize=12
## nb = 2,4 cellsize =6
## nb = 5 cellsize=3
## nb = 6 cellsize =500
## nb = 7 cellsize=100
cellsize = 12 ;

# !!!!! END OF EDIT REGION !!!!!

# Load HDF5 Results from iradina++
ipp = load(['../iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);

#Load data from srim
path = strcat(['../srim/b' num2str(nb) '/']);
%cd (char(path))
%addpath(char(path));
S = srim2mat(path);

#Load data from iradina
path = strcat(['../iradina/b' num2str(nb) '/output/']);
%cd (char(path))
addpath(char(path));
quantities = ['.ions.total'; '.ions.replacements'; '.vac.sum'; '.energy.electronic'; '.energy.phonons' ];

for i=1:size(quantities,1),
fname = strcat('b',num2str(nb),quantities(i,:));
A = load('-ascii',fname );
if i == 1, ions_total = A(:,4)./Nions;,endif;
if i == 2,
  replmnts = A(:,4)./Nions;
  implants = ions_total-replmnts;
endif;
if i == 3,
  if (nb == 3 )|(nb == 4) | (nb == 5),
    B = load('-ascii',strcat('b',num2str(nb),'.vac.z92.m238.029.mat0.elem0'));
    C = load('-ascii',strcat('b',num2str(nb),'.vac.z8.m15.999.mat0.elem1'));
    D = load('-ascii',strcat('b',num2str(nb),'.ions_vac'));
    %vac(:,1) = D(:,5).*5.05e6/Nions; % vac by ions
    vac(:,2) = B(:,4)./Nions; % vac by U
    vac(:,3) = C(:,4)./Nions; % vac by O
    vac(:,4) = A(:,4)/Nions; % vac_sum
    vac(:,1) = vac(:,4)-(vac(:,2)+vac(:,3)) % vac by ions
    %%%
    %for iradina++
   Vac_plus = ipp.tally.defects.Vacancies;
   Vac_all = Vac_plus(:,2) + Vac_plus(:,3); %vacancies by U and O
   else
    vac = A(:,4)/Nions;
    endif;
  endif
if i == 4, ionization = A(:,4)./Nions;,endif;
if i == 5, phonons = A(:,4)./Nions;,endif;
 end
#iradina position
xx = A(:,1)*cellsize;

# Title
titlestr = ipp.Title;
# x axis = cell centers
x = ipp.grid.cell_xyz(:,1);

figure 1
plot(x,ipp.tally.defects.Implantations(:,1),'-o;iradina++;')
hold on
if (nb == 3 )|(nb == 4) | (nb == 5),
  plot(xx,ions_total,'-^;iradina;')
  else
plot(xx,implants,'-^;iradina;')
endif
plot(S.x./10,S.Ri.*cellsize*10e-8,'-d;srim;')
hold off
title([titlestr ' - Impl/I'])


figure 2
if (nb == 3 )|(nb == 4) | (nb == 5),
plot(x,ipp.tally.defects.Vacancies,'-o;iradina++;',x,Vac_all,'-o;iradina++;')
hold on
%plot(S.x./10,S.Vi*cellsize*10,'-d;srim - Xe ions;')
plot(S.x./10,S.Vr(:,1)*cellsize*10,'-d;srim - U ions;')
plot(S.x./10,S.Vr(:,2)*cellsize*10,'-d;srim - O ions;')
plot(S.x./10,(S.Vi+S.Vr(:,1)+S.Vr(:,2))*cellsize*10,'-d;srim - sum all vacancies;')
elseif
plot(x,(ipp.tally.defects.Vacancies(:,1)+ipp.tally.defects.Vacancies(:,2)),'-o;iradina++;')
hold on
plot(S.x./10,(S.Vi+S.Vr)*cellsize*10,'-d;srim - sum all vacancies;')
endif
plot(xx,vac,'-^;iradina;')
hold off
title([titlestr ' - Vacancies'])

figure 3
%if (nb == 3 )|(nb == 4) | (nb == 5),
%  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)+ipp.tally.defects.Replacements(:,3)),'-o;iradina++;')
%elseif
%  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)),'-o;iradina++;')
%end
plot(x,ipp.tally.defects.Replacements(:,1),'-o;iradina++;')
hold on
plot(S.x./10,S.RC*cellsize./Nions,'-d;srim;')
plot(xx,replmnts,'-^;iradina;')
hold off
title([titlestr ' - Replacements/I'])

figure 4
if (nb == 3 )|(nb == 4) | (nb == 5),
  plot(x,(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3)),'-o;Iradina++;')
elseif
  plot(x,(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)),'-o;Iradina++;')
end
hold on
plot(S.x./10,(S.EIi+S.EIr)*cellsize*10,'-d;srim;')
plot(xx,ionization,'-^;iradina;')
hold off
title([titlestr ' - Ionization total'])

figure 5
if (nb == 3 )|(nb == 4) | (nb == 5),
  plot(x,(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)+ipp.tally.energy_deposition.Phonons(:,3)),'-o;Iradina++;')
elseif
  plot(x,(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)),'-o;Iradina++;')
end
hold on
plot(S.x./10,(S.EPi+S.EPr)*cellsize*10,'-d;srim;')
plot(xx,phonons,'-^;iradina;')
hold off
title([titlestr ' - Phonons total'])
