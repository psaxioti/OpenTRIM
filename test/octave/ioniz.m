function Ionization = ioniz(nb, ipp, S, cellsize, ionization, E0, S1)

ener = ipp.tally.energy_deposition;

  % ionization by ions
%disp(sprintf('Quantity\tSRIM\tiradina\tiradina++'))
str = 'Ioniz/ion';
Ionization.Ii(1) = sum(S.EIi*cellsize*10); % srim
Ionization.Ii(2) = NaN; %iradina
Ionization.Ii(3) = sum(ener.Ionization(:,1)); %iradina++
if ~isempty(S1)
    Ionization.Ii(4) = sum(S1.EIi*cellsize*10);  % srim monolayer
endif
str = 'Ioniz/ion';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Ionization.Ii))


% ionization by recoils
str = 'Ioniz/rec';
Ionization.Ir(1) = sum(S.EIr*cellsize*10); % srim
Ionization.Ir(2) = NaN; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Ionization.Ir(3) = sum(ener.Ionization(:,2)+ener.Ionization(:,3)); % iradina++
else
  Ionization.Ir(3) = sum(ener.Ionization(:,2)); % iradina++
endif
if ~isempty(S1)
    Ionization.Ir(4) = sum(S1.EIr*cellsize*10);  % srim monolayer
endif
%disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Ir))

%total ionization
str = 'Ioniz/Tot';
Ionization.It(1) = sum((S.EIr+S.EIi)*cellsize*10); % srim
Ionization.It(2) = sum(ionization); %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Ionization.It(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)+ener.Ionization(:,3)); % iradina++
else
  Ionization.It(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)); % iradina++
endif
if ~isempty(S1)
    Ionization.It(4) = sum((S1.EIr+S1.EIi)*cellsize*10);  % srim monolayer
endif
%disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e",str,Ionization.It))

%total ionization as %
str = 'Ion/Tot(%)';
Ionization.Itp(1) = sum((S.EIr+S.EIi)*cellsize*10)/E0*100; % srim
Ionization.Itp(2) = sum(ionization)/E0*100; %iradina
if (nb == 3 )||(nb == 4)||(nb == 5),
  Ionization.Itp(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2)+ener.Ionization(:,3))/E0*100; % iradina++
else
  Ionization.Itp(3) = sum(ener.Ionization(:,1)+ener.Ionization(:,2))/E0*100; % iradina++
endif
if ~isempty(S1)
    Ionization.Itp(4) = sum((S1.EIr+S.EIi)*cellsize*10)/E0*100;  % srim monolayer
endif
str = 'Ioniz/rec';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Ionization.Ir))

str = 'Ioniz/Tot';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,Ionization.It))

str = 'Ion/Tot(%)';
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Ionization.Itp))
end
