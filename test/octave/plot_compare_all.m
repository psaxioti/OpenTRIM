clear
clf

nb = 6;
[ipp, S, ions_total, replmnts, implants, vac,ionization, phonons, pkai, titlestr,...
 x, cellsize, Eb, S1] = load_files(nb);

E_init = [2e6 2e6 3e6 1e6 3e5 3e6 1e6];
E0 = E_init(nb);

Vac_plus = ipp.tally.defects.Vacancies;
Vac_all = sum(Vac_plus,2); %vacancies by U and O
# x axis = cell centers
x = ipp.grid.cell_xyz(:,1);

figure 1
plot(x,S.Ri.*cellsize*10e-8,'-dc;srim;')
hold on
plot(x,ions_total,'-^r;iradina;')
plot(x,(ipp.tally.defects.Implantations(:,1)+ipp.tally.defects.Replacements(:,1)),'-ob;iradina++;')
hold off
title([titlestr ' - Impl/I'])


figure 2
if (nb == 3 ) || (nb == 4) || (nb == 5),
  plot(x,S.Vr(:,1)*cellsize*10,'-dc;srim - U ions;')
  hold on
  plot(x,S.Vr(:,2)*cellsize*10,'-dc;srim - O ions;')
  plot(x,(S.Vi+S.Vr(:,1)+S.Vr(:,2))*cellsize*10,'-dc;srim - sum all vacancies;')
  plot(x,vac,'-^r;iradina;')
  plot(x,ipp.tally.defects.Vacancies(:,2:end),'-ob;iradina++;',x,Vac_all,'-ob;iradina++;')
elseif
  plot(x,(S.Vi+S.Vr)*cellsize*10,'-dc;srim - sum all vacancies;')
  hold on
  plot(x,vac,'-^r;iradina;')
  plot(x,(ipp.tally.defects.Vacancies(:,1)+ipp.tally.defects.Vacancies(:,2)),'-ob;iradina++;')
endif
hold off
title([titlestr ' - Vacancies'])

figure 3
plot(x,S.RC*cellsize*10,'-dc;srim;')
hold on
plot(x,replmnts,'-^r;iradina;')
if (nb == 3 )||(nb == 4) || (nb == 5),
  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)...
  +ipp.tally.defects.Replacements(:,3)),'-ob;iradina++;')
elseif
  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)),'-ob;iradina++;')
end
hold off
title([titlestr ' - Replacements/I'])

figure 4
plot(x,(S.EIi+S.EIr)*cellsize*10,'-dc;srim;')
hold on
plot(x,ionization,'-^r;iradina;')
if (nb == 3 ) || (nb == 4) || (nb == 5),
  plot(x,(ipp.tally.energy_deposition.Ionization(:,1)+ipp.tally.energy_deposition.Ionization(:,2)...
  +ipp.tally.energy_deposition.Ionization(:,3)),'-ob;Iradina++;')
elseif
  plot(x,(ipp.tally.energy_deposition.Ionization(:,1)+...
  ipp.tally.energy_deposition.Ionization(:,2)),'-ob;Iradina++;')
end
hold off
title([titlestr ' - Ionization total'])

figure 5
plot(x,(S.EPi+S.EPr)*cellsize*10,'-dc;srim;')
hold on
plot(x,phonons,'-^r;iradina;')
if (nb == 3 ) || (nb == 4) || (nb == 5),
  plot(x,(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2)...
  +ipp.tally.energy_deposition.Phonons(:,3)),'-ob;Iradina++;')
elseif
  plot(x,(ipp.tally.energy_deposition.Phonons(:,1)+ipp.tally.energy_deposition.Phonons(:,2))...
  ,'-ob;Iradina++;')
end
hold off
title([titlestr ' - Phonons total'])

% damage energy
Tdam = dam_ener(ipp, S, cellsize, E0, S1);

% PKA
figure 8
plot(x,(S.Vi)*cellsize*10,'-dc;srim;'); % srim pka/ion
hold on
%plot(x,pkai,'-^r;iradina;');
plot(x, sum(ipp.tally.defects.PKAs,2), '-ob;Iradina++;');
if ~isempty(S1)
    plot(x,S1.Vi*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off

