clear

nb =6

qc = srim2mat(['../test/srim/qc/' num2str(nb)]);
fc = srim2mat(['../test/srim/fc2/' num2str(nb)]);
ml = srim2mat(['../test/srim/ml/' num2str(nb) '/1']);
for i=2:10,
  dd = srim2mat(['../test/srim/ml/' num2str(nb) '/' num2str(i)]);
  ml.Vi += dd.Vi;
  ml.Vr += dd.Vr;
  ml.Ri += dd.Ri;
end
ml.Vi /= 10;
ml.Vr /= 10;
ml.Ri /= 10;

dx = fc.x(2)-fc.x(1);
disp('QC')
disp(['pka/ion = ' num2str(sum(qc.Vi)*dx)])
disp(['vac/ion = ' num2str(sum(qc.Vr+qc.Vi)*dx)])
disp('FC')
disp(['pka/ion = ' num2str(sum(fc.Vi)*dx)])
disp(['vac/ion = ' num2str(sum(sum(fc.Vr))*dx)])
disp('ML')
disp(['pka/ion = ' num2str(sum(ml.Vi)*dx)])
disp(['vac/ion = ' num2str(sum(ml.Vr)*dx)])

ipp = load(['../test/iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);
N = double(ipp.Nions);
Ri = double(ipp.Implantations(:,1)+ipp.Replacements(:,1))/N/dx/1e-8;
V = double(sum(ipp.Vacancies(:,2:end),2))/N/dx;
disp('IRADINA++')
disp(['pka/ion = ' num2str(double(ipp.Npkas)/N)])
disp(['Sum Vt = ' num2str(sum(V)*dx)])
disp(['Sum Vnrt = ' num2str(sum(ipp.Vnrt)/N)])
disp(['Sum Vnrt-QC = ' num2str(sum(ipp.Vnrt_LSS)/N)])

x = ipp.cell_xyz(:,1);

figure 1
clf

subplot(2,1,1)
plot(x,qc.Vi+qc.Vr,...
     x,sum(fc.Vr,2),...
     x,ml.Vr,...
     x,[ipp.Vnrt_LSS'/N/dx V ipp.Vnrt'/N/dx],'--')

##subplot(2,1,1)
##plot(x,(qc.Vi+qc.Vr)/(sum(qc.Vr+qc.Vi)*dx),...
##     x,fc.Vr/(sum(fc.Vr)*dx),...
##     x,ml.Vr/(sum(ml.Vr)*dx),...
##     x,V/(sum(V)*dx))

legend('SRIM QC','SRIM FC','SRIM ML',
       'IRADINA++ QC','IRADINA++ FC','IRADINA++ FC-NRT')


subplot(2,1,2)
plot(x,qc.Ri,x,fc.Ri,x,ml.Ri,x,Ri)

legend('SRIM QC','SRIM FC','SRIM ML','IRADINA++')

figure 2
clf

subplot(2,1,1)
plot(x,[ipp.Tdam; ipp.Tdam_LSS]/N,...
     x,sum(ipp.PhononEnergy(:,2:end),2)/N)
title('IRADINA++ Tdam (eV/ion)')
legend('FC Tdam','QC Tdam (LSS)','Phonons by recoils')

subplot(2,1,2)
plot(x,qc.Vi+qc.Vr,...
  x,[ipp.Vnrt_LSS; ipp.Vnrt]/N/dx,...
  x,double(sum(ipp.Vacancies(:,2:end),2))/N/dx)
title('Vacancies/ion')
legend('SRIM QC','IRADINA++ QC',...
  'IRADINA++ NRT-FC','IRADINA++ FC')
