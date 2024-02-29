clear

nb = 7

qc = srim2mat(['../test/srim/qc/' num2str(nb)]);
fc = srim2mat(['../test/srim/fc2/' num2str(nb)]);
ml = srim2mat(['../test/srim/ml2/' num2str(nb)]);

dx = fc.x(2)-fc.x(1);

function y = NRT(T,Ed)
  y = T/Ed;
  i = find(y<1); y(i)=0;
  i = find(y>=1 & y<2.5); y(i) = 1;
  i = find(y>=2.5); y(i) = y(i)/2.5;
end

Eb = 3; Ed=40;
disp( "Vac/Ion   QC\tFC\tML")
disp(sprintf("SRIM      %.1f \t%.1f \t%.1f",...
    sum(qc.Vr+qc.Vi)*dx,...
    sum(sum(fc.Vr))*dx,...
    sum(sum(ml.Vr))*dx))

ipp = load(['../test/iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);
N = double(ipp.Nions);
Ri = double(ipp.Implantations(:,1)+ipp.Replacements(:,1))/N/dx/1e-8;
V = double(sum(ipp.Vacancies(:,2:end),2))/N/dx;
disp(sprintf("IRADINA++ %.1f \t%.1f \t%.1f",...
    sum(ipp.Vnrt_LSS)/N,...
    sum(V)*dx, ...
    sum(ipp.Vnrt)/N))

x = ipp.cell_xyz(:,1);

figure 1
clf

subplot(2,1,1)
plot(x,qc.Vi+qc.Vr,...
     x,sum(fc.Vr,2),...
     x,sum(ml.Vr,2),...
     x,ipp.Vnrt_LSS'/N/dx,'--',...
     x,double(sum(ipp.Vacancies(:,2:end),2))/N/dx,'--', ...
     x,ipp.Vnrt'/N/dx,'--')
title('Vacancy production')
xlabel('nm')
ylabel('Vacancies/ion-A')
legend('SRIM QC','SRIM FC','SRIM ML',...
       'IRADINA++ QC','IRADINA++ FC','IRADINA++ FC-NRT')
xlim([4 8]*1000)


subplot(2,1,2)
Eb = 3;
plot(x,qc.Vi*Eb + qc.Vr*Eb + qc.EPr, ... # YR Lin D3***
     x,fc.Vr*Eb + fc.EPr, ...  # YR Lin D3**
     x,ml.Vr*Eb + ml.EPr, ...  # YR Lin D3**
     x,ipp.Tdam_LSS/N/dx,'--',...
     x,sum(ipp.PhononEnergy(:,2:end),2)/N/dx,'--', ...
     x,ipp.Tdam/N/dx,'--')
title('Tdam')
xlabel('nm')
ylabel('eV/ion-A')
legend('SRIM QC','SRIM FC','SRIM ML',...
       'IRADINA++ QC','IRADINA++ FC','IRADINA++ FC-NRT')
xlim([4 8]*1000)

