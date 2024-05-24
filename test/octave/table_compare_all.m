clear

nb = 6;
[ipp, S, ions_total, replmnts, implants, vac,ionization, phonons, pkai, titlestr,...
 x, cellsize, Eb, S1] = load_files(nb);

E_init = [2e6 2e6 3e6 1e6 3e5 3e6 1e6];
E0 = E_init(nb);

ener = ipp.tally.energy_deposition;

disp(sprintf('Quantity\tSRIM      \tiradina   \tiradina++   \tmonolayer'))

% Implanted Ions
X = impl(S,ions_total, ipp, cellsize, S1);

disp(' ')

%  Phonons
Phonon = phon(nb, ipp, S, cellsize, phonons, E0, S1);

disp(' ')

%  Ionization

Ionization = ioniz(nb, ipp, S, cellsize, ionization, E0, S1);
% ionization by ions
% ionization by recoils
%total ionization
%total ionization as %

disp(' ')

%%% DE %%%
energyD = energ(nb, ipp, S, cellsize, ionization,phonons, E0, S1);


disp(' ')

%Vacancies
Vacan = vacan(nb, ipp, S, cellsize, vac, E0, S1);

%% replacements
Replac = replamnts(nb, ipp, S, cellsize, replmnts, S1);

% repl/vac ratio
str = 'repl/vac';
for i=1:length(Vacan.Vaci),
ratio(i) = Replac.Repl(i)./Vacan.Vaci(i);
end
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,ratio))

disp('')
% damage energy
Tdam = dam_ener(ipp, S, cellsize, E0, S1);

pkas = pka_f(ipp, S, cellsize, E0, pkai, S1);


