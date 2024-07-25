% plot mcnp data
% data are stored in .dat file and are :
% row 1: theta angle(rad) -row 2-4: probab. density/sr acoerding to models fnal1, fnal2 and gauss

% load mcnp data
ions ={'h','he','he2'}
for i=1:length(ions)
  fname = sprintf('../mcnp/mcnp_data_%s.dat',ions{i});
MCNP = load(fname);
M_theta(:,i) = MCNP(:,1);
M_fnal1(:,i) = MCNP(:,2);
M_fnal2(:,i) = MCNP(:,3);
M_gaus(:,i) = MCNP(:,4);
end
titles = {'270 keV H^+ on C, 200k ions','270 keV He^+ on C, 200k ions', '270 keV He^+ on C, ions > 1 M'}
for i = 1:3,
  figure (i)
  plot(M_theta(:,i),M_fnal1(:,i),'-;fnal1;',M_theta(:,i),M_fnal2(:,i),'-;fnal2;',M_theta(:,i),M_gaus(:,i),'-;gauss;')
  xlabel('angle (rad)')
ylabel('Prob. Density (sr^{-1})')
%legend(lbls)
title(titles{i})
grid
print2png(gcf,[4 3]*3,sprintf('msc_%sinC_mcnp',ions{i}))
end
