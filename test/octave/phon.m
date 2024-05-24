function Phonon = phon(nb, ipp, S, cellsize, phonons, E0, S1)

ener = ipp.tally.energy_deposition;

% phonons by ions
str = 'Phon/Ion';
Phonon.Pi(1) = sum(S.EPi*cellsize*10); % srim
Phonon.Pi(2) = NaN; %iradina
Phonon.Pi(3) = sum(ener.Phonons(:,1)); %iradina++
if ~isempty(S1)
    Phonon.Pi(4) = sum(S1.EPi*cellsize*10);  % srim monolayer
  endif
%disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pi))

% phonons by recoils
str = 'Phon/rec';
Phonon.Pr(1) = sum(S.EPr*cellsize*10); % srim
Phonon.Pr(2) = NaN; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Phonon.Pr(3) = sum(ener.Phonons(:,2)+ener.Phonons(:,3)); % iradina++
else
  Phonon.Pr(3) = sum(ener.Phonons(:,2)); % iradina++
endif
if ~isempty(S1)
    Phonon.Pr(4) = sum(S1.EPr*cellsize*10);  % srim monolayer
endif
%disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pr))

%total phonons
str = 'Phon Total';
Phonon.Pt(1) = sum((S.EPr+S.EPi)*cellsize*10); % srim
Phonon.Pt(2) = sum(phonons); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Phonon.Pt(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)+ener.Phonons(:,3)); % iradina++
else
  Phonon.Pt(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)); % iradina++
endif
if ~isempty(S1)
    Phonon.Pt(4) = sum((S1.EPr+S.EPi)*cellsize*10);  % srim monolayer
endif
%disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Pt))

%total phonons
str = 'Phon/Tot(%)';
Phonon.Ptp(1) = sum((S.EPr+S.EPi)*cellsize*10)/E0*100; % srim
Phonon.Ptp(2) = sum(phonons)/E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Phonon.Ptp(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2)+ener.Phonons(:,3))/E0*100; % iradina++
else
  Phonon.Ptp(3) = sum(ener.Phonons(:,1)+ener.Phonons(:,2))/E0*100; % iradina++
endif
if ~isempty(S1)
    Phonon.Ptp(4) = sum((S1.EPr+S.EPi)*cellsize*10)/E0*100;  % srim monolayer
endif

% phonons by ions
str = 'Phon/Ion';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Phonon.Pi))
% phonons by recoils
str = 'Phon/rec';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Phonon.Pr))
% phonons total
str = 'Phon total';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Phonon.Pt))
%total phonons
str = 'Phon/Tot(%)';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Phonon.Ptp))

end
