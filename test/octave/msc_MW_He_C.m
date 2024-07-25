clear

pkg load hdf5oct

th = linspace(0,0.2,100); # degrees
mux = 0.5*(1-cos(th));
dOmega = diff(mux)*4*pi;

lbls = {'ML', 'd/100', 'd/30', ...
        'd/10', 'MH, T>1eV'};

H = zeros(length(lbls),length(mux)-1);
Nc = zeros(1,length(lbls));
FP = zeros(1,length(lbls));
Th = zeros(1,length(lbls));
ips = zeros(1,length(lbls));
tm = zeros(1,length(lbls));

for i=1:length(lbls)
  fname = sprintf('../msc/HeinC/b%d/out.h5',i)
  A = h5read(fname,'/exit_events/event_data');
  j = find(A(5,:)<=440 & A(2,:)==0);
  mu = 0.5*(1-A(8,j));
  Th(i) = mean(sqrt(4*mu));
  h = histc(mu,mux);
  h = h(1:end-1)/size(mu,2);
  H(i,:) = h./dOmega;
  %Th(i) = (H(i,:).*(th(2:end) - th(2)/2))*dOmega'/(H(i,:)*dOmega');
  Nc(i) = h5read(fname,'/tally/ion_stat/collisions')(1);
  FP(i) = h5read(fname,'/tally/ion_stat/flight_path')(1);
  ips(i) = h5read(fname,'/run_stat/ips');
  tm(i) = h5read(fname,'/run_stat/Nh')/ips(i);

end

disp(" ")
disp("MFP        ions/s(rel.) col/ion Th(rad)")
for i=1:length(lbls)
  printf('%10s %12.1f %7.1f %7.3f\n',lbls{i},ips(i)/ips(1),Nc(i),Th(i));
end


figure 1
clf
plot(th(2:end),H,'.-')
xlabel('angle (rad)')
ylabel('Prob. Density (sr^{-1})')
legend(lbls)
title('270 keV He^+ in 100\mum/cm^2 C')
grid

%print2png(gcf,[4 3]*3,'msc_HeinC')



