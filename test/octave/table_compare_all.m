clear
clf

## !!!!! EDIT REGION - EDITABLE VALUES  !!!!
# Benchmark Number
nb = 7;
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
quantities = {'.ions.total'; '.ions.replacements'; ...
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
        %vac(:,1) = D(:,5).*5.05e6/Nions; % vac by ions
        vac(:,2) = B(:,4)./Nions; % vac by U
        vac(:,3) = C(:,4)./Nions; % vac by O
        vac(:,4) = A(:,4)/Nions; % vac_sum
        vac(:,1) = vac(:,4)-(vac(:,2)+vac(:,3)); % vac by ions
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

disp(sprintf('Quantity\tSRIM      \tiradina   \tiradina++'))

% Implanted Ions
str = 'Impl. Ions';
X(1) = sum(S.Ri.*cellsize*10e-8);  % srim
if (nb == 3 )||(nb == 4)||(nb == 5),
  X(2) = sum(ions_total); % iradina
else
  X(2) = sum(implants); % iradina
endif
X(3) = sum(ipp.tally.defects.Implantations(:,1));
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f",str,X))

disp(' ')

%  Phonons
% phonons by ions
str1 = 'Phon/ Ion\t';
str1 = [str1 num2str(sum(S.EPi*cellsize*10),4) '\t\t']; % srim
str1 = [str1 '\t\t']; %iradina
str1 = [str1 num2str(sum(ipp.tally.energy_deposition.Phonons(:,1)),4)]; %iradina++
disp(sprintf(str1))

% phonons by recoils
str2 = 'Phon/rec\t';
str2 = [str2 num2str(sum(S.EPr*cellsize*10),4) '\t\t']; % srim
str2 = [str2 '\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str2 = [str2 num2str(sum(sum(ipp.tally.energy_deposition.Phonons(:,2)+ipp.tally.energy_deposition.Phonons(:,3))),4) '\t']; % iradina++
else
  str2 = [str2 num2str(sum(ipp.tally.energy_deposition.Phonons(:,2)),4) '\t']; % iradina++
endif
disp(sprintf(str2))

%total phonons
str3 = 'Phon/TOT\t';
str3 = [str3 num2str(sum((S.EPr+S.EPi)*cellsize*10),4) '\t']; % srim
str3 = [str3 num2str(sum(phonons),4) '\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)+ipp.tally.energy_deposition.Phonons(:,3)),4)]; % iradina++
else
  str3 = [str3 ...
    num2str(sum(ipp.tally.energy_deposition.Phonons(:,1)+...
    ipp.tally.energy_deposition.Phonons(:,2)),4)]; % iradina++
endif
disp(sprintf(str3))

%total phonons
str3 = 'Phon/TOT(%%)\t';
str3 = [str3 num2str(sum((S.EPr+S.EPi)*cellsize*10)/E0*100,4) '   \t']; % srim
str3 = [str3 num2str(sum(phonons)/E0*100,4) '   \t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)+ipp.tally.energy_deposition.Phonons(:,3))/E0*100,4)]; % iradina++
else
  str3 = [str3 ...
    num2str(sum(ipp.tally.energy_deposition.Phonons(:,1)+...
    ipp.tally.energy_deposition.Phonons(:,2))/E0*100,4)]; % iradina++
endif
disp(sprintf(str3))

disp(' ')

%  Ionization
% ionization by ions
%disp(sprintf('Quantity\tSRIM\tiradina\tiradina++'))
str1 = 'Ioniz/ion\t';
str1 = [str1 num2str(sum(S.EIi*cellsize*10),4) '\t\t']; % srim
str1 = [str1 '\t']; %iradina
str1 = [str1 num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)),4)]; %iradina++
disp(sprintf(str1))

% ionization by recoils
str2 = 'Ioniz/rec\t';
str2 = [str2 num2str(sum(S.EIr*cellsize*10),4) '\t\t']; % srim
str2 = [str2 '\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str2 = [str2 num2str(sum(sum(ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3))),4) '\t']; % iradina++
else
  str2 = [str2 num2str(sum(ipp.tally.energy_deposition.Ionization(:,2)),4) '\t']; % iradina++
endif
disp(sprintf(str2))

%total ionization
str3 = 'Ioniz/TOT\t';
str3 = [str3 num2str(sum((S.EIr+S.EIi)*cellsize*10),4) '\t']; % srim
str3 = [str3 num2str(sum(ionization),4) '\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3)),4)]; % iradina++
else
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)),4)]; % iradina++
endif
disp(sprintf(str3))

%total ionization as %
str3 = 'Ioniz/TOT(%%)\t';
str3 = [str3 num2str(sum((S.EIr+S.EIi)*cellsize*10/E0*100),4) '   \t']; % srim
str3 = [str3 num2str(sum(ionization)/E0*100,4) '   \t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3))/E0*100,4) '\n']; % iradina++
else
  str3 = [str3 num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2))/E0*100,4) '\n']; % iradina++
endif
disp(sprintf(str3))

disp(' ')

%%% DE %%%
disp(sprintf('Quantity\tSRIM\t\tiradina\t\tiradina++'))
str = 'E_c\t\t';
str = [str num2str(sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10),4) '\t']; % srim
str = [str num2str(sum(ionization+phonons),4) '\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str = [str num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)+...
  ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3)+...
  ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)...
  +ipp.tally.energy_deposition.Phonons(:,3)),4) '\t']; % iradina++
else
  str = [str num2str(sum(ipp.tally.energy_deposition.Ionization(:,1)...
  +ipp.tally.energy_deposition.Ionization(:,2)+...
  ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)),4) '\t']; % iradina++
endif
disp(sprintf(str))

str = 'DE\t\t';
str = [str num2str(((E0 - sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10))./E0)*100,4) '\t\t']; % srim
str = [str num2str((E0 - sum(ionization+phonons))./E0*100,4) '\t\t']; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  str = [str num2str((E0 - sum(ipp.tally.energy_deposition.Ionization(:,1)+...
  ipp.tally.energy_deposition.Ionization(:,2)+ipp.tally.energy_deposition.Ionization(:,3)+...
  ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)...
  +ipp.tally.energy_deposition.Phonons(:,3)))./E0*100,4) '\t']; % iradina++
else
  str = [str num2str((E0 - sum(ipp.tally.energy_deposition.Ionization(:,1)...
  +ipp.tally.energy_deposition.Ionization(:,2)+...
  ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)))./E0*100,4) '\t']; % iradina++
endif
disp(sprintf(str))
