%plot all data
clf
folders = {'../srim/msc/FC/HinC/', '../srim/msc/monolayer/HinC/' , '../srim/msc/QC/HinC/';...
 '../srim/msc/FC/HeinC', '../srim/msc/monolayer/HeinC/', '../srim/msc/QC/HeinC/'};


for i=1:2,
  figure (i)
  for m = 1:3,
  cd(folders{i,m});
  addpath('../../')
  cosines = load('cosines.dat');
  cosx = cosines(:,1);

range = 0:0.002:0.2;
binranges = fliplr(cos(range));
bincnts = fliplr(histc(cosx,binranges)');
bincnts = bincnts(2:end);
lims = abs(diff(fliplr(binranges)));
d = semilogy(range(2:end),bincnts/200000./(2*pi()*lims),'linestyle','-','linewidth',2)
hold on
cd('../../../../octave/')
end
hold off
xlabel('angle (rad)')
ylabel('Prob. Density (sr^{-1})')
legend('FC','monolayer','QC','location', 'northeast')
grid on
if i ==1,
  title('270 keV H^+ on C')
  print2png(gcf,[4 3]*3,'mscSrim_HinC')
else
  title('270 keV He^+ on C')
  print2png(gcf,[4 3]*3,'mscSrim_HeinC')
end

end
