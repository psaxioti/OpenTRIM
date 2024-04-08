clear

# Benchmark Number
nb = 3;

# Load HDF5 Results
ipp = load(['../iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);

# Title
titlestr = ipp.Title;

# x axis = cell centers
x = ipp.grid.cell_xyz(:,1);


# Figure 1 = Defects vs. x
figure 1
clf
subplot(2,2,1)
plot(x,ipp.tally.defects.Implantations)
title([titlestr ' - Impl/I'])
legend(ipp.atom.label)
xlabel('x (nm)')
subplot(2,2,2)
plot(x,ipp.tally.defects.Vacancies)
title([titlestr ' - V'])
legend(ipp.atom.label)
xlabel('x (nm)')
subplot(2,2,3)
plot(x,ipp.tally.defects.Replacements)
title([titlestr ' - R'])
legend(ipp.atom.label)
xlabel('x (nm)')
subplot(2,2,4)
plot(x,ipp.tally.defects.PKAs)
title([titlestr ' - PKAs'])
legend(ipp.atom.label)
xlabel('x (nm)')

# Figure 2 = Energy Deposition vs. x
figure 2
clf
subplot(2,2,1)
plot(x,ipp.tally.energy_deposition.Ionization)
title([titlestr ' - Ionization'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('eV')
subplot(2,2,2)
plot(x,ipp.tally.energy_deposition.Phonons)
title([titlestr ' - Phonons'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('eV')
subplot(2,2,3)
plot(x,ipp.tally.energy_deposition.PKA)
title([titlestr ' - Er'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('eV')
subplot(2,2,4)
plot(x,(ipp.tally.energy_deposition.PKA)./(ipp.tally.defects.PKAs))
title([titlestr ' - Er/PKA'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('eV')

# Figure 3 = Damage parameters vs. x
figure 3
clf
lbls = {};
for i=2:size(ipp.atom.label,1),
  lbls = { lbls{:}, [strtrim(ipp.atom.label(i,:)) ' FC-NRT']};
end
lbls = { lbls{:}, 'Total FC-NRT' };
for i=2:size(ipp.atom.label,1),
  lbls = { lbls{:}, [strtrim(ipp.atom.label(i,:)) ' FC']};
end
lbls = { lbls{:}, 'Total FC' };

subplot(2,2,1)
plot(x,ipp.tally.damage.Tdam(:,2:end), ...
  x,sum(ipp.tally.damage.Tdam(:,2:end),2), ...
  x,ipp.tally.energy_deposition.Phonons(:,2:end),...
  x,sum(ipp.tally.energy_deposition.Phonons(:,2:end),2))
title([titlestr ' - Tdam'])
legend(lbls)
xlabel('x (nm)')
ylabel('eV')
subplot(2,2,2)
plot(x,ipp.tally.damage.Vnrt(:,2:end), ...
  x,sum(ipp.tally.damage.Vnrt(:,2:end),2), ...
  x,ipp.tally.defects.Vacancies(:,2:end),...
  x,sum(ipp.tally.defects.Vacancies(:,2:end),2))
title([titlestr ' - Vac'])
legend(lbls)
xlabel('x (nm)')


lbls = {};
for i=2:size(ipp.atom.label,1),
  lbls = { lbls{:}, [strtrim(ipp.atom.label(i,:)) ' FC-NRT']};
end
lbls = { lbls{:}, 'Total FC-NRT', 'Total FC' };
pka = ipp.tally.defects.PKAs(:,2:end);

subplot(2,2,3)
plot(x,ipp.tally.damage.Tdam(:,2:end)./pka, ...
  x,sum(ipp.tally.damage.Tdam(:,2:end),2)./sum(pka,2), ...
  x,sum(ipp.tally.energy_deposition.Phonons(:,2:end),2)./sum(pka,2))
title([titlestr ' - Tdam/PKA'])
legend(lbls)
xlabel('x (nm)')
ylabel('eV')
subplot(2,2,4)
plot(x,ipp.tally.damage.Vnrt(:,2:end)./pka,...
     x,sum(ipp.tally.damage.Vnrt(:,2:end),2)./sum(pka,2),...
     x,sum(ipp.tally.defects.Vacancies(:,2:end),2)./sum(pka,2))
title([titlestr ' - Vac/PKA'])
legend(lbls)
xlabel('x (nm)')

# Figure 4 = Lost ions/energy vs. x
figure 4
clf
subplot(2,1,1)
plot(x,ipp.tally.energy_deposition.Lost)
title([titlestr ' - Lost E'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('eV')
subplot(2,1,2)
plot(x,ipp.tally.defects.Lost)
title([titlestr ' - Lost ions'])
legend(ipp.atom.label)
xlabel('x (nm)')







