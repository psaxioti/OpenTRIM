clear
clf
pkg load hdf5oct

th = linspace(0,0.2,41); # degrees
mux = 0.5*(1-cos(th));
dOmega = diff(mux)*4*pi;

lbls = {'ML', 'd/100', 'd/30', ...
        'd/10', 'MH, T>1eV'};

H = zeros(length(lbls),length(mux)-1);
Nc = zeros(1,length(lbls)-3);
FP = zeros(1,length(lbls)-3);
Th = zeros(1,length(lbls)-3);
ips = zeros(1,length(lbls)-3);
tm = zeros(1,length(lbls)-3);

for i=1:length(lbls)
  fname = sprintf('../msc/HinC/b%d/out.h5',i)
  A = h5read(fname,'/exit_events/event_data');
  j = find(A(5,:)>=440 & A(2,:)==0);
  mu = 0.5*(1-A(8,j));
  Th(i) = mean(sqrt(4*mu));
  h = histc(mu,mux)/size(mu,2);
  H(i,:) = h(1:end-1)./dOmega;
  %Th(i) = (H(i,:).*(th(2:end) - 0*th(2)/2))*dOmega'/(H(i,:)*dOmega');
  Nc(i) = h5read(fname,'/tally/ion_stat/collisions')(1);
  FP(i) = h5read(fname,'/tally/ion_stat/flight_path')(1);
  ips(i) = h5read(fname,'/run_stat/ips');
  tm(i) = h5read(fname,'/run_stat/Nh')/ips(i);

end

disp(" ")
disp("MFP        ions/s(rel.) col/ion Th(rad)")
for i=1:length(lbls)-3
  printf('%10s %12.1f %7.1f %7.3f\n',lbls{i},ips(i)/ips(1),Nc(i),Th(i));
end


% load SRIM data
folders = {'../srim/msc/FC/HinC/', '../srim/msc/monolayer/HinC/' , '../srim/msc/QC/HinC/'};
for m = 1:3,
  cd(folders{m});
  addpath('../../')
  cosines = load('cosines.dat');
  cosx{m} = cosines(:,1);

range = 0:0.002:0.2;
binranges = fliplr(cos(range));
bincnts1 = fliplr(histc(cosx{m}(:,:),binranges)');
bincnts(:,m) = bincnts1(2:end);
lims = abs(diff(fliplr(binranges)));
cd('../../../../octave/')
end

% load mcnp data
MCNP = load('../mcnp/mcnp_data_h.dat');
M_theta = MCNP(:,1);
M_fnal1 = MCNP(:,2);
M_fnal2 = MCNP(:,3);
M_gaus = MCNP(:,4);

figure 1
clf
semilogy(th(2:end),H,'.-')% iradina++
hold on
semilogy(range(2:end),bincnts'/200000./(2*pi()*lims),'linestyle','--','linewidth',2)%srim
semilogy(M_theta,M_fnal1,'o-',M_theta,M_fnal2,'-',M_theta,M_gaus,'-')% mcnp
hold off
xlabel('angle (rad)')
ylabel('Prob. Density (sr^{-1})')
lbls2 = {'ML', 'd/100', 'd/30', ...
        'd/10', 'MH, T>1eV',...
        'FC','monolayer','QC',...
        'fnal1','fnal2','gauss'};
legend(lbls2)
title('270 keV H^+ in 100\mum/cm^2 C')
grid

print2png(gcf,[4 3]*4,'msc_HinC_all')



