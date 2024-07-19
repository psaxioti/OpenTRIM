function Vacan = vacan(nb, ipp, S, cellsize, vac, E0, S1)

%Vac by ions
defects = ipp.tally.defects;
if (nb == 3 ) || (nb == 4) || (nb == 5),
  str = 'Vac U/ion';
  Vacan.Vaci(1) = sum(S.Vr(:,1))*cellsize*10;% U srim
  Vacan.Vaci(2) = sum(vac(:,1)); % U iradina
  Vacan.Vaci(3) = sum(defects.Vacancies(:,2));% % U iradina++
  if ~isempty(S1)
    Vacan.Vaci(4) = sum(S1.Vr(:,1))*cellsize*10;  % srim monolayer
  endif
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Vacan.Vaci))

  % O vacancies
  str = 'Vac O/ion';
  Vacan.Vaci(1) = sum(S.Vr(:,2))*cellsize*10;% O srim
  Vacan.Vaci(2) = sum(vac(:,2)); % O iradina
  Vacan.Vaci(3) = sum(defects.Vacancies(:,3));% % O iradina++
  if ~isempty(S1)
     Vacan.Vaci(1) = sum(S1.Vr(:,2))*cellsize*10;  % O srim monolayer
  endif
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Vacan.Vaci))

    % sum vacancies
  str = 'Vac sum/ion';
  Vacan.Vaci(1) = sum(S.Vr(:,2)+S.Vr(:,1))*cellsize*10;% sum srim
  Vacan.Vaci(2) = sum(vac(:,3)); % sum iradina
  Vacan.Vaci(3) = sum(defects.Vacancies(:,3)+defects.Vacancies(:,2));% % sum iradina++
  if ~isempty(S1)
     Vacan.Vaci(4) = sum(S1.Vr(:,2)+S1.Vr(:,1))*cellsize*10; % srim monolayer
  endif
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Vacan.Vaci))
else
  str = 'Vac sum/ion';
  Vacan.Vaci(1) = sum(S.Vr*cellsize*10); %srim
  Vacan.Vaci(2) = sum(vac); % iradina
  Vacan.Vaci(3) = sum(defects.Vacancies(:,2));% iradina ++
  if ~isempty(S1)
     Vacan.Vaci(4) = sum(S1.Vr*cellsize*10); % srim monolayer
  endif
  disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Vacan.Vaci))
endif

end
