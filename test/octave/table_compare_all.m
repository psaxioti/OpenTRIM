clear
clf

## !!!!! EDIT REGION - EDITABLE VALUES  !!!!
# Benchmark Number
nb = 1;
# !!!!! END OF EDIT REGION !!!!!

# Number of ions run
Nions_all = [20000 20 8000 20000 20000 20000 20000];
Nions = Nions_all(nb);

# Size of the cell
cellsize_all = [12 6 12 6 3 500 100];
cellsize = cellsize_all(nb) ;

Energy_all = [2e6 2e6 3e6 1e6 3e5 3e6 1e6]; %ion energy in eV
E0 = Energy_all(nb);

# Load HDF5 Results from iradina++
ipp = load(['../iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);

#Load data from srim
S = srim2mat(['../srim/b' num2str(nb) '/']);

#Load data from iradina
path = ['../iradina/b' num2str(nb) '/output/'];
quantities = {'.ions.total'; '.repl.sum'; ...
  '.vac.sum'; '.energy.electronic'; '.energy.phonons'};

for i=1:size(quantities,1),
  A = load('-ascii',[path 'b' num2str(nb) quantities{i}]);
  switch i,
    case 1
      ions_total = A(:,4)./Nions;
    case 2
      replmnts = A(:,4)./Nions;
      implants = ions_total-replmnts;
    case 3
      if (nb == 3 ) || (nb == 4) || (nb == 5),
        B = load('-ascii',[path 'b' num2str(nb) '.vac.z92.m238.029.mat0.elem0']);
        C = load('-ascii',[path 'b' num2str(nb) '.vac.z8.m15.999.mat0.elem1']);
        D = load('-ascii',[path 'b' num2str(nb) '.ions_vac']);
        %vac(:,1) = D(:,5).*5.05e6/Nions;
        vac(:,1) = B(:,4)./Nions; % vac by U
        vac(:,2) = C(:,4)./Nions; % vac by O
        vac(:,3) = A(:,4)/Nions; % vac_sum
        %vac(:,1) = vac(:,4)-(vac(:,2)+vac(:,3));
        %%%
        %for iradina++
        Vac_plus = ipp.tally.defects.Vacancies;
        Vac_all = Vac_plus(:,2) + Vac_plus(:,3); %vacancies by U and O
      else
        vac = A(:,4)/Nions;
      end
    case 4
      ionization = A(:,4)./Nions;
    case 5
      phonons = A(:,4)./Nions;
    otherwise
      error('i case not handled');
  endswitch
end

ener = ipp.tally.energy_deposition;

disp(sprintf('Quantity\tSRIM      \tiradina   \tiradina++'))

% Implanted Ions
str = 'Impl. Ions';
X(1) = sum(S.Ri.*cellsize*10e-8);  % srim
X(2) = sum(ions_total); % iradina
X(3) = sum(ipp.tally.defects.Implantations(:,1)+ipp.tally.defects.Replacements(:,1)); %iradina++
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,X))

disp(' ')

%  Phonons
% phonons by ions
str = 'Phon/Ion';
Pi(1) = sum(S.EPi*cellsize*10); % srim
Pi(2) = NaN; %iradina
Pi(3) = sum(ener.Phonons(:,1)); %iradina++
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pi))

% phonons by recoils
str = 'Phon/rec';
Pr(1) = sum(S.EPr*cellsize*10); % srim
Pr(2) = NaN; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Pr(3) = sum(ener.Phonons(:,2)+ener.Phonons(:,3)); % iradina++
else
  Pr(3) = sum(ener.Phonons(:,2)); % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pr))

%total phonons
str = 'Phon Total';
Pt(1) = sum((S.EPr+S.EPi)*cellsize*10); % srim
Pt(2) = sum(phonons); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Pt(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)+ener.Phonons(:,3)); % iradina++
else
  Pt(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)); % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pt))

%total phonons
str = 'Phon/Tot(%)';
Ptp(1) = sum((S.EPr+S.EPi)*cellsize*10)/E0*100; % srim
Ptp(2) = sum(phonons)/E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Ptp(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)+ener.Phonons(:,3))/E0*100; % iradina++
else
  Ptp(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2))/E0*100; % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Ptp))

disp(' ')

%  Ionization
% ionization by ions
%disp(sprintf('Quantity\tSRIM\tiradina\tiradina++'))
str = 'Ioniz/ion';
Ii(1) = sum(S.EIi*cellsize*10); % srim
Ii(2) = NaN; %iradina
Ii(3) = sum(ener.Ionization(:,1)); %iradina++
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Ii))

% ionization by recoils
str = 'Ioniz/rec';
Ir(1) = sum(S.EIr*cellsize*10); % srim
Ir(2) = NaN; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Ir(3) = sum(ener.Ionization(:,2)+ener.Ionization(:,3)); % iradina++
else
  Ir(3) = sum(ener.Ionization(:,2)); % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Ir))

%total ionization
str = 'Ioniz/Tot';
It(1) = sum((S.EIr+S.EIi)*cellsize*10); % srim
It(2) = sum(ionization); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  It(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)+ener.Ionization(:,3)); % iradina++
else
  It(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)); % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,It))

%total ionization as %
str = 'Ion/Tot(%)';
Itp(1) = sum((S.EIr+S.EIi)*cellsize*10)/E0*100; % srim
Itp(2) = sum(ionization)/E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Itp(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)+ener.Ionization(:,3))/E0*100; % iradina++
else
  Itp(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2))/E0*100; % iradina++
endif
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Itp))

disp(' ')

%%% DE %%%
%disp(sprintf('Quantity\tSRIM\t\tiradina\t\tiradina++'))
str = 'Total ener';
de1(1) = sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10); % srim
de1(2) = sum(ionization+phonons); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  de1(3) = sum(ener.Ionization(:,1)+...
  ener.Ionization(:,2)+ener.Ionization(:,3)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)...
  +ener.Phonons(:,3)); % iradina++
else
  de1(3) = sum(ener.Ionization(:,1)...
  +ener.Ionization(:,2)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)); % iradina++
endif
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,de1))

str = 'EnerDif(%)';
de(1) = ((E0 - sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10))./E0)*100; % srim
de(2) = (E0 - sum(ionization+phonons))./E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  de(3) = (E0 - sum(ener.Ionization(:,1)+...
  ener.Ionization(:,2)+ener.Ionization(:,3)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)...
  +ener.Phonons(:,3)))./E0*100; % iradina++
else
  de(3) = (E0 - sum(ener.Ionization(:,1)...
  +ener.Ionization(:,2)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)))/E0*100; % iradina++
endif
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,de))

disp(' ')
%%% Annotation
% X = Implanted Ions
% Pi = Phonon energy by ions (eV)
% Pr = Phonon energy by recoils (eV)
% Pt = total phonon energy (eV)
% Ptp = persentage of total energy difference from E0
% Ietc... equivalent for ionization energy

%Vacancies
%Vac by ions
defects = ipp.tally.defects;
if (nb == 3 ) || (nb == 4) || (nb == 5),
  str = 'Vac U/ion';
  Vaci(1) = sum(S.Vr(:,1))*cellsize*10;% U srim
  Vaci(2) = sum(vac(:,1)); % U iradina
  Vaci(3) = sum(defects.Vacancies(:,2));% % U iradina++
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Vaci))

  % O vacancies
  str = 'Vac O/ion';
  Vaci(1) = sum(S.Vr(:,2))*cellsize*10;% O srim
  Vaci(2) = sum(vac(:,2)); % O iradina
  Vaci(3) = sum(defects.Vacancies(:,3));% % O iradina++
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Vaci))

    % sum vacancies
  str = 'Vac sum/ion';
  Vaci(1) = sum(S.Vr(:,2)+S.Vr(:,1))*cellsize*10;% sum srim
  Vaci(2) = sum(vac(:,3)); % sum iradina
  Vaci(3) = sum(defects.Vacancies(:,3)+defects.Vacancies(:,2));% % sum iradina++
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Vaci))
elseif
  str = 'Vac sum/ion';
  Vaci(1) = sum(S.Vr*cellsize*10); %srim
  Vaci(2) = sum(vac); % iradina
  Vaci(3) = sum(defects.Vacancies(:,2));% iradina ++
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Vaci))
endif

%% replacements
str = 'repl/ion';
Repl(1) = sum(S.RC*cellsize*10);% srim
Repl(2) = sum(replmnts);% iradina
if (nb == 3 )|(nb == 4) | (nb == 5),
  Repl(3) = sum(defects.Replacements(:,1)+defects.Replacements(:,2)+defects.Replacements(:,3)); %iradina++
elseif
  Repl(3) = sum(defects.Replacements(:,1)+defects.Replacements(:,2)) %iradina++
end
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,Repl))

% repl/vac ratio
for i=1:3,
ratio(i) = Repl(i)./Vaci(i);
end
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,ratio))



