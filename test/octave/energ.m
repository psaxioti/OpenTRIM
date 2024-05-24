function energyD = energ(nb, ipp, S, cellsize, ionization,phonons, E0, S1)

  ener = ipp.tally.energy_deposition;

energyD.de1(1) = sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10); % srim
energyD.de1(2) = sum(ionization+phonons); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  energyD.de1(3) = sum(ener.Ionization(:,1)+...
  ener.Ionization(:,2)+ener.Ionization(:,3)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)...
  +ener.Phonons(:,3)); % iradina++
else
  energyD.de1(3) = sum(ener.Ionization(:,1)...
  +ener.Ionization(:,2)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)); % iradina++
endif
if ~isempty(S1)
    energyD.de1(4) = sum((S1.EIr+S1.EIi+S1.EPi+S1.EPr)*cellsize*10);  % srim monolayer
endif

energyD.de(1) = ((E0 - sum((S.EIr+S.EIi+S.EPi+S.EPr)*cellsize*10))./E0)*100; % srim
energyD.de(2) = (E0 - sum(ionization+phonons))./E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  energyD.de(3) = (E0 - sum(ener.Ionization(:,1)+...
  ener.Ionization(:,2)+ener.Ionization(:,3)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)...
  +ener.Phonons(:,3)))./E0*100; % iradina++
else
  energyD.de(3) = (E0 - sum(ener.Ionization(:,1)...
  +ener.Ionization(:,2)+...
  ener.Phonons(:,1)+ener.Phonons(:,2)))/E0*100; % iradina++
endif
if ~isempty(S1)
    energyD.de(4) = ((E0 - sum((S1.EIr+S1.EIi+S1.EPi+S1.EPr)*cellsize*10))./E0)*100;  % srim monolayer
endif
str = 'Total ener';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,energyD.de1))
str = 'EnerDif(%)';
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,energyD.de))
end
