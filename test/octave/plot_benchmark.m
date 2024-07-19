clear

# Benchmark Number
nb = 6;

# Load HDF5 Results
pkg load hdf5oct
ipp = h5load(['../iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);

# Title
titlestr = ipp.title;

# x axis = cell centers
x = ipp.grid.cell_xyz(:,1);


# Figure 1 = Defects vs. x
figure 1
clf
subplot(2,2,1)
plot(x,ipp.tally.defects.Implantations(:,1))
title([titlestr ' - Implanted Ions'])
legend(ipp.atom.label{1})
xlabel('x (nm)')
subplot(2,2,2)
plot(x,ipp.tally.defects.Implantations(:,2:end))
title([titlestr ' - I'])
legend(ipp.atom.label{2:end})
xlabel('x (nm)')
subplot(2,2,3)
plot(x,ipp.tally.defects.Vacancies)
title([titlestr ' - V'])
legend(ipp.atom.label)
xlabel('x (nm)')
subplot(2,2,4)
plot(x,ipp.tally.defects.Replacements)
title([titlestr ' - R'])
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
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' Tdam']};
end
lbls = { lbls{:}, 'Total Tdam' };
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' Eph']};
end
lbls = { lbls{:}, 'Total Eph' };

subplot(2,2,1)
plot(x,ipp.tally.damage.Tdam(:,2:end), ...
  x,sum(ipp.tally.damage.Tdam(:,2:end),2), ...
  x,ipp.tally.energy_deposition.Phonons(:,2:end),...
  x,sum(ipp.tally.energy_deposition.Phonons(:,2:end),2))
title([titlestr ' - Tdam'])
legend(lbls)
xlabel('x (nm)')
ylabel('eV')


lbls = {};
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' FC-NRT']};
end
lbls = { lbls{:}, 'Total FC-NRT' };
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' FC']};
end
lbls = { lbls{:}, 'Total FC' };

subplot(2,2,2)
plot(x,ipp.tally.damage.Vnrt(:,2:end), ...
  x,sum(ipp.tally.damage.Vnrt(:,2:end),2), ...
  x,ipp.tally.defects.Vacancies(:,2:end),...
  x,sum(ipp.tally.defects.Vacancies(:,2:end),2))
title([titlestr ' - Vac'])
legend(lbls)
xlabel('x (nm)')


lbls = {};
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' Tdam']};
end
lbls = { lbls{:}, 'Total Tdam', 'Total Eph' };
pka = ipp.tally.defects.PKAs(:,2:end);

subplot(2,2,3)
plot(x,ipp.tally.damage.Tdam(:,2:end)./pka, ...
  x,sum(ipp.tally.damage.Tdam(:,2:end),2)./sum(pka,2), ...
  x,sum(ipp.tally.energy_deposition.Phonons(:,2:end),2)./sum(pka,2))
title([titlestr ' - Tdam/PKA'])
legend(lbls)
xlabel('x (nm)')
ylabel('eV')

lbls = {};
for i=2:length(ipp.atom.label),
  lbls = { lbls{:}, [ipp.atom.label{i} ' FC-NRT']};
end
lbls = { lbls{:}, 'Total FC-NRT', 'Total FC' };

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
subplot(2,2,1)
plot(x,ipp.tally.ion_stat.collisions)
title([titlestr ' - Collisions'])
legend(ipp.atom.label)
xlabel('x (nm)')
subplot(2,2,2)
plot(x,ipp.tally.defects.PKAs)
title([titlestr ' - PKAs'])
legend(ipp.atom.label)
xlabel('x (nm)')

subplot(2,2,3)
plot(x,ipp.tally.ion_stat.flight_path./ipp.tally.ion_stat.collisions)
title([titlestr ' - mfp'])
legend(ipp.atom.label)
xlabel('x (nm)')
ylabel('nm')
subplot(2,2,4)
plot(x,ipp.tally.defects.Lost)
title([titlestr ' - Lost ions'])
legend(ipp.atom.label)
xlabel('x (nm)')

# PDF report name
pname = ['../iradina++/b' num2str(nb) '/b' num2str(nb) '.pdf'];

% Get the table in string form.
TString = sprintf('{\\bf\\fontsize{20}%s}\n\n',ipp.title);
TString = [TString evalc(['table_compare_all(' num2str(nb) ')'])];


% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');

figure 5
clf
set(gcf,'papertype','a4','paperorientation','landscape')
% Output the table using the annotation command.
annotation('textbox',[0.3 0 1 1],'string',TString,...
    'FontName',FixedWidth,'Units','Normalized',...
    "verticalalignment",'middle',...
    'fitboxtotext','off','edgecolor','none');

print(gcf,'-dpdf','-fillpage',pname)

for i=1:4
  figure(i)
  set(gcf,'papertype','a4','paperorientation','landscape')
  print(gcf,'-dpdf','-fillpage','-append',pname)
endfor








