clear
clf

## !!!!! EDIT REGION - EDITABLE VALUES  !!!!
# Benchmark Number
nb = 3;
# !!!!! END OF EDIT REGION !!!!!

# Number of ions run
Nions_all = [20000 20 8000 20000 20000 20000 20000];
Nions = Nions_all(nb);

# Size of the cell
cellsize_all = [12 6 12 6 3 500 100];
cellsize = cellsize_all(3) ;

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

disp(sprintf('Quantity\tSRIM\tiradina\tiradina++'))

% Implanted Ions
str = 'Impl. Ions\t';
str = [str num2str(sum(S.Ri.*cellsize*10e-8),4) '\t']; % srim
if (nb == 3 )||(nb == 4)||(nb == 5),
  str = [str num2str(sum(ions_total),4) '\t']; % iradina
else
  str = [str num2str(sum(implants),4) '\t']; % iradina
endif
str = [str num2str(sum(ipp.tally.defects.Implantations(:,1)),4)]; % iradina++
disp(sprintf(str))

% Total Phonons


