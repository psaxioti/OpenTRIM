%% Test on Lin et al. appendix A2 methods
nb = 6;
[ipp, S, ions_total, replmnts, implants, vac,ionization, phonons, titlestr,...
 x, cellsize, Eb] = load_files(nb);
%initial beam energy
E_init = [2e6 2e6 3e6 1e6 3e5 3e6 1e6];
E0 = E_init(nb);

energy = ipp.tally.energy_deposition;
Ionizpp = energy.Ionization;
Phonpp = energy.Phonons;
pka = energy.PKA;
defects = ipp.tally.defects;
% methods
% D1
% Tdam = beam energy absorbed by target atoms-target atom energy lost to ionization
TdamS(1) = (sum(sum(S.ERr(1:end))) - sum(S.EIr))*cellsize*10;
TdamPP(1) = sum(sum(pka(:,2:end))) - sum(sum(Ionizpp(:,2:end))); %%% ?????

% D2
% Tdam = Insidend beam energy - beam energy lost to ionization -
% target atom energy lost to ionization - beam energy lost to phonons
TdamS(2) = E0 - sum(S.EIi + S.EIr + S.EPi)*cellsize*10;
TdamPP(2) = E0 - sum(Ionizpp(:,1) +Phonpp(:,1)) - sum(sum(Ionizpp(:,2:end)));

% D3
% Tdam = Target atom energy lost to phonons - beam energy lost to phonons
TdamS(3) = sum(S.EPr - S.EPi)*cellsize*10;
TdamPP(3) = sum(sum(Phonpp(:,2:end))) - sum(Phonpp(:,1));

% D1*
% Tdam = Beam energy absorbed by target atoms - target atom energy lost to
% to ionization + beam energy lost to phonon production
TdamS(4) = (sum(sum(S.ERr(:,1:end))) - sum(S.EIr) + sum(S.EPi))*cellsize*10; % ?
TdamPP(4) = sum(sum(pka(:,2:end))) - sum(sum(Ionizpp(:,2:end)))...
 + sum(Phonpp(:,1)); % ????

% D2* = D3* (according to my understanding)
% Tdam = beam energy lost to lattice binding energy + target ato energy
% lost to lattice binding energy (Eb*\nu) + beam energy lost to phonon
%  production + target atom energy lost to phonons
TdamS(5) = (Eb*sum(S.Vi) + Eb*sum(sum(S.Vr(:,1:end)))+...
    sum(S.EPi+S.EPr))*cellsize*10;
TdamPP(5) = Eb*sum(sum(defects.PKAs(:,2:end)))+ ...
           Eb*sum(sum(defects.Vacancies(:,2:end)))+...
           sum(sum(Phonpp(:,1:end)));

% D2**
% Tdam = Incident ion beam energy - beam energy lost to ionization -
% target atom energy lost to ionization
TdamS(6) = E0 - sum(S.EIi + S.EIr)*cellsize*10;
TdamPP(6) = E0 - sum(Ionizpp(:,1)) - sum(sum(Ionizpp(:,2:end)));

% D3**
% Tdam = Target atom energy lost to lattice binding energy + target atom
% energy lost to phonons
TdamS(7) = (Eb*sum(sum(S.Vr(:,1:end)))+  ...
        sum(sum(S.EPr(:,1:end))))*cellsize*10;
TdamPP(7) = Eb*sum(sum(defects.Vacancies(:,2:end)))+...
           sum(sum(Phonpp(:,2:end)));

% D3***
% Tdam = beam energy lost to lattice binding energy + target atom energy
% lost to lattice binding energy + target atom energy lost to phonons
TdamS(8) = (Eb*(sum(S.Vi) + sum(sum(S.Vr(:,1:end))))+  ...
        sum(sum(S.EPr(:,1:end))))*cellsize*10;
TdamPP(8) = Eb*(sum(defects.Vacancies(:,1)) + sum(sum(defects.Vacancies(:,2:end))))+...
           sum(sum(Phonpp(:,2:end)));

% Reference value by iradina++ Tdam
Tdamppref = sum(sum(ipp.tally.damage.Tdam(:,2:end)))

disp(sprintf('Code\tD1 \t\tD2 \t\tD3 \t\tD1* \t\tD2*=D3* \tD2** \t\tD3** \t\tD3***'))
strs = 'Srim';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",strs,TdamS))
strs = 'Comp';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",strs,TdamS./Tdamppref))
strs = 'Irad++';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",strs,TdamPP))
strs = 'Comp';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",strs,TdamPP./Tdamppref))
