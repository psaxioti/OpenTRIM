clear

# Effect of correlated recombination

run_opentrim = 0;

# Load HDF5 Results
pkg load hdf5oct

cfg = jsondecode(fileread('../opentrim/b1/config.json'));

Nh = 1000;

cfg.Driver.max_no_ions = Nh;
cfg.Driver.threads = 8;
cfg.Simulation.eloss_calculation = "EnergyLossOff";
cfg.Simulation.intra_cascade_recombination = true;
cfg.Target.materials.composition.El = 0.001;

cfg.Simulation.correlated_recombination = false;

Ed = 40;
cfg.Target.materials.composition.Ed = Ed;
cfg.Target.materials.composition.Er = Ed;
cfg.Output.OutputFileBaseName = 'outEd40eVNoCorr';

fid = fopen ("cfg.json", "w");
fputs (fid, jsonencode(cfg));
fclose (fid);

if run_opentrim,
  system('opentrim -f cfg.json');
end

function [nv,nr]=makehist(Ed,bins,D)
  Td = D(5,:)';
  V = D(6,:)';
  Rep = D(7,:)';
  R = D(9,:)';

  td = Td;
  v=V;
  r=R;
  j = find(td < Ed);
  td(j)=[];
  v(j)=[];
  r(j)=[];

  nv = zeros(size(bins));
  nr = zeros(size(bins));

  for i=1:length(bins),
    j = find(td <= bins(i));
    nv(i) = mean(v(j));
    nr(i) = mean(r(j));

    td(j) = [];
    v(j) = [];
    r(j) = [];
  end

end

E = logspace(1,6,51);

D = h5load('outEd40eV.h5','/events/pka/event_data');
[V40,R40]=makehist(40,E,D);
D = h5load('outEd40eVNoCorr.h5','/events/pka/event_data');
[V40nc,R40nc]=makehist(40,E,D);

nrt = E/100;
j=find(E<100);
nrt(j)=1;
j=find(E<40);
nrt(j)=1e-6;

xi = nrt;
j=find(E>100);
c = 0.286;
b = -0.568;
xi(j) = (1-c)*nrt(j).^b + c;

figure 1
clf

loglog(E,V40+R40,'.-', ...
  E,V40,'.-',...
  E,V40nc,'.-',...
  E,nrt,E,nrt.*xi)
ylim([0.1 1e4])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T_d$ (eV)','interpreter','latex')
ylabel('$N_d$','interpreter','latex')

h = legend('OpenTRIM - no ICR, $E_d=40$ eV', ...
  'OpenTRIM - ICR, $E_d=40$ eV', ...
  'OpenTRIM - ICR, No Corr. R.', ...
  'NRT, $E_d=40$ eV', 'arc-dpa', ...
  'location','northwest');
set(h,'interpreter','latex')

title('2MeV Fe in Fe','interpreter','latex')

R_c = cfg.Target.materials.composition.Rc;
text(1e3,1,['$R_c$(nm) = ' num2str(R_c) ' '],'interpreter','latex')

sz = [4 3]*4;
set(1,'PaperUnits','centimeters')
set(1,'PaperSize',sz)
set(1,'PaperPosition',[0 0 sz(1) sz(2)])

print(1,'-dpng','icr2')

figure 2
clf

semilogx(E,(V40nc-V40)./R40)
ylim([0 1.1])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T_d$ (eV)','interpreter','latex')



title('2MeV Fe in Fe, fraction of correlated recombination',...
  'interpreter','latex')




