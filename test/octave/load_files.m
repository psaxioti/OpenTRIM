function [ipp, S, ions_total, replmnts, implants, vac, ionization, phonons, pkai, ...
  titlestr, x, cellsize, Eb, S1] = load_files(nb)

## !!!!! EDIT REGION - EDITABLE VALUES  !!!!
# Benchmark Number
%nb = 1;
# !!!!! END OF EDIT REGION !!!!!

# Number of ions run
Nions_all = [20000 20 8000 20000 20000 20000 20000];
Nions = Nions_all(nb);

# multipling factor for pkai
fact = [2.99e6 1478 3.32e6 5.052e6 2.54e6 70052 30025];
factor = fact(nb);

# Size of the cell
cellsize_all = [12 6 12 6 3 500 100];
cellsize = cellsize_all(nb) ;

% binding Energy
Ebind = [3, 3, 3.3, 3.3, 3.3, 3, 3];
Eb = Ebind(nb);

# Load HDF5 Results from iradina++
ipp = load(['../iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);

#Load data from srim
S = srim2mat(['../srim/b' num2str(nb) '/']);
if (nb == 6) || (nb ==7),
  S1 = srim2mat(['../srim/ml/b' num2str(nb) '/']);,
else
  S1 = [];
end

#Load data from iradina
path = ['../iradina/b' num2str(nb) '/output/'];
quantities = {'.ions.total'; '.repl.sum'; ...
  '.vac.sum'; '.energy.electronic'; '.energy.phonons'; '.ions_vac'};

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
        vac(:,1) = B(:,4)./Nions; % vac U
        vac(:,2) = C(:,4)./Nions; % vac O
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
    case 6
      pkai = A(:,5).*factor./Nions;
      %implants2 = A(:,4)/Nions;
      %size(A)
    otherwise
      error('i case not handled');
  endswitch
end

#iradina position
xx = A(:,1)*cellsize;

# Title
titlestr = ipp.Title;
# x axis = cell centers
x = ipp.grid.cell_xyz(:,1);

%% replacement check
load('-ascii',[path 'b' num2str(nb) quantities{i}]);

end
