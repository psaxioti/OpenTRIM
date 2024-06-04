clear

th = linspace(0,0.2,41); # degrees
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
  fname = sprintf('../msc/HinC/b%d/out.exit.h5',i)
  A = load(fname);
  j = find(A.exit(5,:)>=440 & A.exit(2,:)==0);
  mu = 0.5*(1-A.exit(8,j));
  h = histc(mu,mux)/size(mu,2);
  H(i,:) = h(1:end-1)./dOmega;
  Th(i) = (H(i,:).*(th(2:end) - 0*th(2)/2))*dOmega'/(H(i,:)*dOmega');
  fname = sprintf('../msc/HinC/b%d/out.h5',i)
  A = load(fname);
  Nc(i) = A.tally.ion_stat.collisions(1);
  FP(i) = A.tally.ion_stat.flight_path(1);
  ips(i) = A.run_stat.ips;
  tm(i) = A.run_stat.Nh/ips(i);

end

num2str([ips'/ips(1) tm' Nc' FP' Th'])

figure 1
clf
plot(th(2:end),sqrt(H),'.-')
xlabel('angle (rad)')
ylabel('Prob. Density (sr^{-1})')
legend(lbls)
title('270 keV H^+ in 100\mum/cm^2 C')
grid

print2png(gcf,[12 9],'msc_HinC')
