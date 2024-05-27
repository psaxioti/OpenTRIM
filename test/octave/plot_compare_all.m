clear

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
clf
plot(x,S.Ri.*cellsize*10e-8,'-dc;srim;')
hold on
plot(x,ions_total,'-^r;iradina;')
plot(x,(ipp.tally.defects.Implantations(:,1)+ipp.tally.defects.Replacements(:,1)),'-ob;iradina++;')
if ~isempty(S1)
    plot(x,S1.Ri.*cellsize*10e-8,'-dk;srim mnl;');  % srim monolayer
endif
hold off
title([titlestr ' - Impl/I'])
xlabel('x (nm)','interpreter','latex')
ylabel(['Implants per ion'], ...
        'interpreter','latex')
box on



figure 2
clf
if (nb == 3 ) || (nb == 4) || (nb == 5),
  plot(x,S.Vr(:,1)*cellsize*10,'-dc;srim - U ions;')
  hold on
  plot(x,S.Vr(:,2)*cellsize*10,'-dc;srim - O ions;')
  plot(x,(S.Vr(:,1)+S.Vr(:,2))*cellsize*10,'-dc;srim - sum all vacancies;')
  plot(x,vac,'-^r;iradina;')
  plot(x,ipp.tally.defects.Vacancies(:,2:end),'-ob;iradina++;',x,Vac_all,'-ob;iradina++;')
elseif
  plot(x,(S.Vr)*cellsize*10,'-dc;srim - sum all vacancies;')
  hold on
  plot(x,vac,'-^r;iradina;')
  plot(x,(ipp.tally.defects.Vacancies(:,1)+ipp.tally.defects.Vacancies(:,2)),'-ob;iradina++;')
endif
if ~isempty(S1)
    plot(x,(S1.Vr)*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off
title([titlestr ' - Vacancies'])
xlabel('x (nm)','interpreter','latex')
ylabel(['Vacancies per ion'], ...
        'interpreter','latex')
box on

figure 3
clf
plot(x,S.RC*cellsize*10,'-dc;srim;')
hold on
plot(x,replmnts,'-^r;iradina;')
if (nb == 3 )||(nb == 4) || (nb == 5),
  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)...
  +ipp.tally.defects.Replacements(:,3)),'-ob;iradina++;')
elseif
  plot(x,(ipp.tally.defects.Replacements(:,1)+ipp.tally.defects.Replacements(:,2)),'-ob;iradina++;')
end
if ~isempty(S1)
    plot(x,S1.RC*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off
title([titlestr ' - Replacements/I'])
xlabel('x (nm)','interpreter','latex')
ylabel(['replacements per ion'], ...
        'interpreter','latex')
box on

figure 4
clf
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
if ~isempty(S1)
    plot(x,(S1.EIi+S1.EIr)*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off
title([titlestr ' - Ionization total'])
xlabel('x (nm)','interpreter','latex')
ylabel(['Ionization energy (eV)'], ...
        'interpreter','latex')
box on

figure 5
clf
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
if ~isempty(S1)
    plot(x,(S1.EPi+S1.EPr)*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off
title([titlestr ' - Phonons total'])
xlabel('x (nm)','interpreter','latex')
ylabel(['Phonon energy (eV)'], ...
        'interpreter','latex')
box on

% damage energy
Tdam = dam_ener(ipp, S, cellsize, E0, S1);

% PKA
figure 8
clf
plot(x,(S.Vi)*cellsize*10,'-dc;srim;'); % srim pka/ion
hold on
%plot(x,pkai,'-^r;iradina;');
plot(x, sum(ipp.tally.defects.PKAs,2), '-ob;Iradina++;');
if ~isempty(S1)
    plot(x,S1.Vi*cellsize*10,'-dk;srim mnl;');  % srim monolayer
endif
hold off
box on
title([titlestr ' - PKAs/Ion'])
xlabel('x (nm)','interpreter','latex')
ylabel(['Number of PKAs per ion'], ...
        'interpreter','latex')
